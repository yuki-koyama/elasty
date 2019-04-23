#include <elasty/cloth-sim-object.hpp>
#include <elasty/constraint.hpp>
#include <elasty/particle.hpp>
#include <iostream>
#include <Eigen/Geometry>
#include <tiny_obj_loader.h>

elasty::ClothSimObject::ClothSimObject(const std::string& obj_path,
                                       const double distance_stiffness,
                                       const double bending_stiffness,
                                       const Eigen::Affine3d& transform,
                                       const Strategy strategy)
{
    tinyobj::attrib_t attrib;
    std::vector<tinyobj::shape_t> shapes;
    std::vector<tinyobj::material_t> materials;

    std::string warn;
    std::string err;
    const bool return_value = tinyobj::LoadObj(&attrib, &shapes, &materials, &warn, &err, obj_path.c_str());

    if (!warn.empty()) { std::cerr << warn << std::endl; }
    if (!err.empty()) { std::cerr << err << std::endl; }
    if (!return_value) { throw std::runtime_error(""); }

    if (attrib.vertices.empty() || attrib.normals.empty()) { throw std::runtime_error(""); }

    const auto& shape = shapes[0];

    assert(shapes.size() == 1);
    assert(attrib.vertices.size() % 3 == 0);
    assert(shape.mesh.indices.size() % 3 == 0);

    m_triangle_list.resize(shape.mesh.indices.size() / 3, 3);
    for (unsigned int i = 0; i < shape.mesh.indices.size() / 3; ++ i)
    {
        m_triangle_list(i, 0) = shape.mesh.indices[i * 3 + 0].vertex_index;
        m_triangle_list(i, 1) = shape.mesh.indices[i * 3 + 1].vertex_index;
        m_triangle_list(i, 2) = shape.mesh.indices[i * 3 + 2].vertex_index;
    }

    std::map<unsigned int, std::shared_ptr<Particle>> map_from_obj_vertex_index_to_particle;
    for (unsigned int i = 0; i < attrib.vertices.size() / 3; ++ i)
    {
        const Eigen::Vector3d position
        {
            attrib.vertices[3 * i + 0],
            attrib.vertices[3 * i + 1],
            attrib.vertices[3 * i + 2]
        };

        const Eigen::Vector3d x = transform * position;
        const Eigen::Vector3d v = 20.0 * Eigen::Vector3d::Random();
        const double m = 1.0 / double(attrib.vertices.size());

        auto particle = std::make_shared<elasty::Particle>(x, v, m);

        map_from_obj_vertex_index_to_particle[i] = particle;

        m_particles.push_back(particle);
    }

    for (unsigned int i = 0; i < shape.mesh.indices.size() / 3; ++ i)
    {
        const auto p_0 = map_from_obj_vertex_index_to_particle[shape.mesh.indices[i * 3 + 0].vertex_index];
        const auto p_1 = map_from_obj_vertex_index_to_particle[shape.mesh.indices[i * 3 + 1].vertex_index];
        const auto p_2 = map_from_obj_vertex_index_to_particle[shape.mesh.indices[i * 3 + 2].vertex_index];

        m_constraints.push_back(std::make_shared<elasty::DistanceConstraint>(p_0, p_1, distance_stiffness, (p_0->x - p_1->x).norm()));
        m_constraints.push_back(std::make_shared<elasty::DistanceConstraint>(p_0, p_2, distance_stiffness, (p_0->x - p_2->x).norm()));
        m_constraints.push_back(std::make_shared<elasty::DistanceConstraint>(p_1, p_2, distance_stiffness, (p_1->x - p_2->x).norm()));
    }

    using vertex_t = unsigned int;
    using triangle_t = unsigned int;
    using edge_t = std::pair<vertex_t, vertex_t>;
    std::map<edge_t, std::vector<triangle_t>> edges_and_triangles;
    for (unsigned int i = 0; i < shape.mesh.indices.size() / 3; ++ i)
    {
        const vertex_t index_0 = shape.mesh.indices[i * 3 + 0].vertex_index;
        const vertex_t index_1 = shape.mesh.indices[i * 3 + 1].vertex_index;
        const vertex_t index_2 = shape.mesh.indices[i * 3 + 2].vertex_index;

        const edge_t e_01 = std::make_pair(std::min(index_0, index_1), std::max(index_0, index_1));
        const edge_t e_02 = std::make_pair(std::min(index_0, index_2), std::max(index_0, index_2));
        const edge_t e_12 = std::make_pair(std::min(index_1, index_2), std::max(index_1, index_2));

        auto register_edge = [&](const edge_t& edge)
        {
            if (edges_and_triangles.find(edge) == edges_and_triangles.end())
            {
                edges_and_triangles[edge] =  { i };
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
        const edge_t& edge = key_value.first;
        const std::vector<triangle_t>& triangles = key_value.second;

        assert(triangles.size() == 1 || triangles.size() == 2);

        // Boundary
        if (triangles.size() == 1) { continue; }

        auto obtain_another_vertex = [&](const triangle_t& triangle, const edge_t& edge)
        {
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

        switch (strategy)
        {
            case Strategy::Bending:
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
                const double dihedral_angle = std::acos(std::max(- 1.0, std::min(+ 1.0, n_0.dot(n_1))));

                assert(!std::isnan(dihedral_angle));

                m_constraints.push_back(std::make_shared<elasty::BendingConstraint>(p_0, p_1, p_2, p_3, bending_stiffness, dihedral_angle));

                break;
            }
            case Strategy::IsometricBending:
            {
                const auto p_0 = map_from_obj_vertex_index_to_particle[edge.first];
                const auto p_1 = map_from_obj_vertex_index_to_particle[edge.second];
                const auto p_2 = map_from_obj_vertex_index_to_particle[another_vertex_0];
                const auto p_3 = map_from_obj_vertex_index_to_particle[another_vertex_1];

                m_constraints.push_back(std::make_shared<elasty::IsometricBendingConstraint>(p_0, p_1, p_2, p_3, bending_stiffness));

                break;
            }
            case Strategy::Cross:
            {
                const auto p_2 = map_from_obj_vertex_index_to_particle[another_vertex_0];
                const auto p_3 = map_from_obj_vertex_index_to_particle[another_vertex_1];

                const Eigen::Vector3d& x_2 = p_2->x;
                const Eigen::Vector3d& x_3 = p_3->x;

                m_constraints.push_back(std::make_shared<elasty::DistanceConstraint>(p_2, p_3, bending_stiffness, (x_2 - x_3).norm()));

                break;
            }
        }
    }
}
