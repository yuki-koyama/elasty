#include <Eigen/Geometry>
#include <algorithm>
#include <elasty/cloth-sim-object.hpp>
#include <elasty/constraint.hpp>
#include <elasty/particle.hpp>
#include <elasty/utils.hpp>
#include <iostream>
#include <tiny_obj_loader.h>

elasty::ClothSimObject::ClothSimObject(const unsigned           resolution,
                                       const double             in_plane_stiffness,
                                       const double             in_plane_compliance,
                                       const double             out_of_plane_stiffness,
                                       const double             out_of_plane_compliance,
                                       const double             dt,
                                       const Eigen::Affine3d&   transform,
                                       const InPlaneStrategy    in_plane_strategy,
                                       const OutOfPlaneStrategy out_of_plane_strategy)
{
    std::istringstream obj_data_stream(generateClothMeshObjData(2.0, 2.0, resolution, resolution));

    tinyobj::attrib_t                attrib;
    std::vector<tinyobj::shape_t>    shapes;
    std::vector<tinyobj::material_t> materials;

    std::string warn;
    std::string err;
    const bool  return_value = tinyobj::LoadObj(&attrib, &shapes, &materials, &warn, &err, &obj_data_stream);

    if (!warn.empty())
    {
        std::cerr << warn << std::endl;
    }
    if (!err.empty())
    {
        std::cerr << err << std::endl;
    }
    if (!return_value)
    {
        throw std::runtime_error("");
    }
    if (attrib.vertices.empty())
    {
        throw std::runtime_error("");
    }

    const auto& shape = shapes[0];

    assert(shapes.size() == 1);
    assert(attrib.vertices.size() % 3 == 0);
    assert(shape.mesh.indices.size() % 3 == 0);

    m_triangle_list.resize(shape.mesh.indices.size() / 3, 3);
    for (unsigned int i = 0; i < shape.mesh.indices.size() / 3; ++i)
    {
        m_triangle_list(i, 0) = shape.mesh.indices[i * 3 + 0].vertex_index;
        m_triangle_list(i, 1) = shape.mesh.indices[i * 3 + 1].vertex_index;
        m_triangle_list(i, 2) = shape.mesh.indices[i * 3 + 2].vertex_index;
    }

    if (!attrib.texcoords.empty())
    {
        m_uv_list.resize(shape.mesh.indices.size() / 3, 2 * 3);
        for (unsigned int i = 0; i < shape.mesh.indices.size() / 3; ++i)
        {
            m_uv_list(i, 2 * 0 + 0) = attrib.texcoords[2 * shape.mesh.indices[i * 3 + 0].texcoord_index + 0];
            m_uv_list(i, 2 * 0 + 1) = attrib.texcoords[2 * shape.mesh.indices[i * 3 + 0].texcoord_index + 1];
            m_uv_list(i, 2 * 1 + 0) = attrib.texcoords[2 * shape.mesh.indices[i * 3 + 1].texcoord_index + 0];
            m_uv_list(i, 2 * 1 + 1) = attrib.texcoords[2 * shape.mesh.indices[i * 3 + 1].texcoord_index + 1];
            m_uv_list(i, 2 * 2 + 0) = attrib.texcoords[2 * shape.mesh.indices[i * 3 + 2].texcoord_index + 0];
            m_uv_list(i, 2 * 2 + 1) = attrib.texcoords[2 * shape.mesh.indices[i * 3 + 2].texcoord_index + 1];
        }
    }

    std::map<unsigned int, std::shared_ptr<Particle>> map_from_obj_vertex_index_to_particle;
    for (unsigned int i = 0; i < attrib.vertices.size() / 3; ++i)
    {
        const Eigen::Vector3d position{
            attrib.vertices[3 * i + 0], attrib.vertices[3 * i + 1], attrib.vertices[3 * i + 2]};

        const Eigen::Vector3d x = transform * position;
        const Eigen::Vector3d v = Eigen::Vector3d::Zero();
        const double          m = 1.0 / double(attrib.vertices.size());

        auto particle = std::make_shared<elasty::Particle>(x, v, m);

        map_from_obj_vertex_index_to_particle[i] = particle;

        m_particles.push_back(particle);
    }

    for (unsigned int i = 0; i < shape.mesh.indices.size() / 3; ++i)
    {
        const auto p_0 = map_from_obj_vertex_index_to_particle[shape.mesh.indices[i * 3 + 0].vertex_index];
        const auto p_1 = map_from_obj_vertex_index_to_particle[shape.mesh.indices[i * 3 + 1].vertex_index];
        const auto p_2 = map_from_obj_vertex_index_to_particle[shape.mesh.indices[i * 3 + 2].vertex_index];

        if (in_plane_strategy == InPlaneStrategy::EdgeDistance || in_plane_strategy == InPlaneStrategy::Both)
        {
            m_constraints.push_back(std::make_shared<elasty::DistanceConstraint>(
                p_0, p_1, in_plane_stiffness, in_plane_compliance, dt, (p_0->x - p_1->x).norm()));
            m_constraints.push_back(std::make_shared<elasty::DistanceConstraint>(
                p_0, p_2, in_plane_stiffness, in_plane_compliance, dt, (p_0->x - p_2->x).norm()));
            m_constraints.push_back(std::make_shared<elasty::DistanceConstraint>(
                p_1, p_2, in_plane_stiffness, in_plane_compliance, dt, (p_1->x - p_2->x).norm()));
        }

        if (in_plane_strategy == InPlaneStrategy::ContinuumTriangle || in_plane_strategy == InPlaneStrategy::Both)
        {
            m_constraints.push_back(std::make_shared<elasty::ContinuumTriangleConstraint>(
                p_0, p_1, p_2, in_plane_stiffness, in_plane_compliance, dt, 1000.0, 0.10));
        }
    }

    using vertex_t   = unsigned int;
    using triangle_t = unsigned int;
    using edge_t     = std::pair<vertex_t, vertex_t>;
    std::map<edge_t, std::vector<triangle_t>> edges_and_triangles;
    for (unsigned int i = 0; i < shape.mesh.indices.size() / 3; ++i)
    {
        const vertex_t index_0 = shape.mesh.indices[i * 3 + 0].vertex_index;
        const vertex_t index_1 = shape.mesh.indices[i * 3 + 1].vertex_index;
        const vertex_t index_2 = shape.mesh.indices[i * 3 + 2].vertex_index;

        const edge_t e_01 = std::make_pair(std::min(index_0, index_1), std::max(index_0, index_1));
        const edge_t e_02 = std::make_pair(std::min(index_0, index_2), std::max(index_0, index_2));
        const edge_t e_12 = std::make_pair(std::min(index_1, index_2), std::max(index_1, index_2));

        auto register_edge = [&](const edge_t& edge) {
            if (edges_and_triangles.find(edge) == edges_and_triangles.end())
            {
                edges_and_triangles[edge] = {i};
            }
            else
            {
                edges_and_triangles[edge].push_back(i);
            }
        };

        register_edge(e_01);
        register_edge(e_02);
        register_edge(e_12);
    }

    for (const auto& key_value : edges_and_triangles)
    {
        const edge_t&                  edge      = key_value.first;
        const std::vector<triangle_t>& triangles = key_value.second;

        assert(triangles.size() == 1 || triangles.size() == 2);

        // Boundary
        if (triangles.size() == 1)
        {
            continue;
        }

        auto obtain_another_vertex = [&](const triangle_t& triangle, const edge_t& edge) {
            const vertex_t vertex_0 = shape.mesh.indices[3 * triangle + 0].vertex_index;
            const vertex_t vertex_1 = shape.mesh.indices[3 * triangle + 1].vertex_index;
            const vertex_t vertex_2 = shape.mesh.indices[3 * triangle + 2].vertex_index;

            if (vertex_0 != edge.first && vertex_0 != edge.second)
            {
                return vertex_0;
            }
            else if (vertex_1 != edge.first && vertex_1 != edge.second)
            {
                return vertex_1;
            }
            else
            {
                assert(vertex_2 != edge.first && vertex_2 != edge.second);
                return vertex_2;
            }
        };

        const triangle_t another_vertex_0 = obtain_another_vertex(triangles[0], edge);
        const triangle_t another_vertex_1 = obtain_another_vertex(triangles[1], edge);

        switch (out_of_plane_strategy)
        {
            case OutOfPlaneStrategy::Bending:
            {
                const auto p_0 = map_from_obj_vertex_index_to_particle[edge.first];
                const auto p_1 = map_from_obj_vertex_index_to_particle[edge.second];
                const auto p_2 = map_from_obj_vertex_index_to_particle[another_vertex_0];
                const auto p_3 = map_from_obj_vertex_index_to_particle[another_vertex_1];

                const Eigen::Vector3d& x_0 = p_0->x;
                const Eigen::Vector3d& x_1 = p_1->x;
                const Eigen::Vector3d& x_2 = p_2->x;
                const Eigen::Vector3d& x_3 = p_3->x;

                const Eigen::Vector3d p_10 = x_1 - x_0;
                const Eigen::Vector3d p_20 = x_2 - x_0;
                const Eigen::Vector3d p_30 = x_3 - x_0;

                const Eigen::Vector3d n_0 = p_10.cross(p_20).normalized();
                const Eigen::Vector3d n_1 = p_10.cross(p_30).normalized();

                assert(!n_0.hasNaN());
                assert(!n_1.hasNaN());

                // Typical value is 0.0 or pi
                const double dihedral_angle = std::acos(std::min(1.0, std::max(n_0.dot(n_1), -1.0)));

                assert(!std::isnan(dihedral_angle));

                m_constraints.push_back(std::make_shared<elasty::BendingConstraint>(
                    p_0, p_1, p_2, p_3, out_of_plane_stiffness, out_of_plane_compliance, dt, dihedral_angle));

                break;
            }
            case OutOfPlaneStrategy::IsometricBending:
            {
                const auto p_0 = map_from_obj_vertex_index_to_particle[edge.first];
                const auto p_1 = map_from_obj_vertex_index_to_particle[edge.second];
                const auto p_2 = map_from_obj_vertex_index_to_particle[another_vertex_0];
                const auto p_3 = map_from_obj_vertex_index_to_particle[another_vertex_1];

                m_constraints.push_back(std::make_shared<elasty::IsometricBendingConstraint>(
                    p_0, p_1, p_2, p_3, out_of_plane_stiffness, out_of_plane_compliance, dt));

                break;
            }
            case OutOfPlaneStrategy::Cross:
            {
                const auto p_2 = map_from_obj_vertex_index_to_particle[another_vertex_0];
                const auto p_3 = map_from_obj_vertex_index_to_particle[another_vertex_1];

                const Eigen::Vector3d& x_2 = p_2->x;
                const Eigen::Vector3d& x_3 = p_3->x;

                m_constraints.push_back(std::make_shared<elasty::DistanceConstraint>(
                    p_2, p_3, out_of_plane_stiffness, out_of_plane_compliance, dt, (x_2 - x_3).norm()));

                break;
            }
        }
    }

    calculateAreas();
}

void elasty::ClothSimObject::applyAerodynamicForces(const Eigen::Vector3d& global_velocity,
                                                    const double           drag_coeff,
                                                    const double           lift_coeff)
{
    constexpr double rho = 1.225; // Taken from Wikipedia: https://en.wikipedia.org/wiki/Density_of_air

    assert(drag_coeff >= lift_coeff);

    const int num_triangles = m_triangle_list.rows();

    for (int i = 0; i < num_triangles; ++i)
    {
        const auto& x_0 = m_particles[m_triangle_list.row(i)(0)]->x;
        const auto& x_1 = m_particles[m_triangle_list.row(i)(1)]->x;
        const auto& x_2 = m_particles[m_triangle_list.row(i)(2)]->x;

        const auto& v_0 = m_particles[m_triangle_list.row(i)(0)]->v;
        const auto& v_1 = m_particles[m_triangle_list.row(i)(1)]->v;
        const auto& v_2 = m_particles[m_triangle_list.row(i)(2)]->v;

        const auto& m_0 = m_particles[m_triangle_list.row(i)(0)]->m;
        const auto& m_1 = m_particles[m_triangle_list.row(i)(1)]->m;
        const auto& m_2 = m_particles[m_triangle_list.row(i)(2)]->m;

        const double m_sum = m_0 + m_1 + m_2;

        // Calculate the weighted average of the particle vecities
        const Eigen::Vector3d v_triangle = (m_0 * v_0 + m_1 * v_1 + m_2 * v_2) / m_sum;

        // Calculate the relative velocity of the triangle
        const Eigen::Vector3d v_rel         = v_triangle - global_velocity;
        const double          v_rel_squared = v_rel.squaredNorm();

        const auto            cross         = (x_1 - x_0).cross(x_2 - x_0);
        const double          area          = 0.5 * cross.norm();
        const auto            n_either_side = cross.normalized();
        const Eigen::Vector3d n             = (n_either_side.dot(v_rel) > 0.0) ? n_either_side : -n_either_side;

        const double coeff = 0.5 * rho * area;

        // Note: This wind force model was proposed by [Wilson+14]
        const Eigen::Vector3d f =
            -coeff * ((drag_coeff - lift_coeff) * v_rel.dot(n) * v_rel + lift_coeff * v_rel_squared * n);

        m_particles[m_triangle_list.row(i)(0)]->f += (m_0 / m_sum) * f;
        m_particles[m_triangle_list.row(i)(1)]->f += (m_1 / m_sum) * f;
        m_particles[m_triangle_list.row(i)(2)]->f += (m_2 / m_sum) * f;
    }
}

void elasty::ClothSimObject::calculateAreas()
{
    const int num_triangles = m_triangle_list.rows();

    m_area_list = Eigen::VectorXd(num_triangles);

    for (int i = 0; i < num_triangles; ++i)
    {
        const auto& x_0 = m_particles[m_triangle_list.row(i)[0]]->x;
        const auto& x_1 = m_particles[m_triangle_list.row(i)[1]]->x;
        const auto& x_2 = m_particles[m_triangle_list.row(i)[2]]->x;

        const double area = 0.5 * (x_1 - x_0).cross(x_2 - x_0).norm();

        m_area_list(i) = area;
    }
}
