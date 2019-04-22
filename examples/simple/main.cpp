#include <array>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <bigger/app.hpp>
#include <bigger/scene-object.hpp>
#include <bigger/materials/blinnphong-material.hpp>
#include <bigger/primitives/plane-primitive.hpp>
#include <bigger/primitives/sphere-primitive.hpp>
#include <Eigen/Geometry>
#include <elasty/constraint.hpp>
#include <elasty/engine.hpp>
#include <elasty/particle.hpp>
#include <string-util.hpp>
#include <tiny_obj_loader.h>

using ParticlePtr = std::shared_ptr<elasty::Particle>;
using ConstraintPtr = std::shared_ptr<elasty::Constraint>;
using TriangleIndices = Eigen::Matrix<uint16_t, Eigen::Dynamic, 3, Eigen::RowMajor>;

class SimObject
{
public:

    std::vector<ParticlePtr> m_particles;
    std::vector<ConstraintPtr> m_constraints;
};

class ClothSimObject : public SimObject
{
public:

    ClothSimObject()
    {
        constexpr double cloth_distance_stiffness = 0.20;
        constexpr double cloth_bending_stiffness = 0.05;

        const Eigen::Affine3d cloth_import_transform = Eigen::Translation3d(1.0, 1.0, 0.0) * Eigen::AngleAxisd(0.5 * glm::pi<double>(), Eigen::Vector3d::UnitX());

        tinyobj::attrib_t attrib;
        std::vector<tinyobj::shape_t> shapes;
        std::vector<tinyobj::material_t> materials;

        std::string warn;
        std::string err;
        const bool return_value = tinyobj::LoadObj(&attrib, &shapes, &materials, &warn, &err, "models/cloths/0.01.obj");

        if (!warn.empty()) { std::cerr << warn << std::endl; }
        if (!err.empty()) { std::cerr << err << std::endl; }
        if (!return_value) { throw std::runtime_error(""); }

        if (attrib.vertices.empty() || attrib.normals.empty()) { throw std::runtime_error(""); }

        const auto& shape = shapes[0];

        assert(shapes.size() == 1);
        assert(attrib.vertices.size() % 3 == 0);
        assert(shape.mesh.indices.size() % 3 == 0);

        m_triangle_indices.resize(shape.mesh.indices.size() / 3, 3);
        for (unsigned int i = 0; i < shape.mesh.indices.size() / 3; ++ i)
        {
            m_triangle_indices(i, 0) = shape.mesh.indices[i * 3 + 0].vertex_index;
            m_triangle_indices(i, 1) = shape.mesh.indices[i * 3 + 1].vertex_index;
            m_triangle_indices(i, 2) = shape.mesh.indices[i * 3 + 2].vertex_index;
        }

        std::map<unsigned int, ParticlePtr> map_from_obj_vertex_index_to_particle;
        for (unsigned int i = 0; i < attrib.vertices.size() / 3; ++ i)
        {
            const Eigen::Vector3d position
            {
                attrib.vertices[3 * i + 0],
                attrib.vertices[3 * i + 1],
                attrib.vertices[3 * i + 2]
            };

            const Eigen::Vector3d x = cloth_import_transform * position;
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

            m_constraints.push_back(std::make_shared<elasty::DistanceConstraint>(p_0, p_1, cloth_distance_stiffness, (p_0->x - p_1->x).norm()));
            m_constraints.push_back(std::make_shared<elasty::DistanceConstraint>(p_0, p_2, cloth_distance_stiffness, (p_0->x - p_2->x).norm()));
            m_constraints.push_back(std::make_shared<elasty::DistanceConstraint>(p_1, p_2, cloth_distance_stiffness, (p_1->x - p_2->x).norm()));
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

            enum class Strategy { Bending, IsometricBending, Cross };
            constexpr Strategy strategy = Strategy::IsometricBending;

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

                    m_constraints.push_back(std::make_shared<elasty::BendingConstraint>(p_0, p_1, p_2, p_3, cloth_bending_stiffness, dihedral_angle));

                    break;
                }
                case Strategy::IsometricBending:
                {
                    const auto p_0 = map_from_obj_vertex_index_to_particle[edge.first];
                    const auto p_1 = map_from_obj_vertex_index_to_particle[edge.second];
                    const auto p_2 = map_from_obj_vertex_index_to_particle[another_vertex_0];
                    const auto p_3 = map_from_obj_vertex_index_to_particle[another_vertex_1];

                    m_constraints.push_back(std::make_shared<elasty::IsometricBendingConstraint>(p_0, p_1, p_2, p_3, cloth_bending_stiffness));

                    break;
                }
                case Strategy::Cross:
                {
                    const auto p_2 = map_from_obj_vertex_index_to_particle[another_vertex_0];
                    const auto p_3 = map_from_obj_vertex_index_to_particle[another_vertex_1];

                    const Eigen::Vector3d& x_2 = p_2->x;
                    const Eigen::Vector3d& x_3 = p_3->x;

                    m_constraints.push_back(std::make_shared<elasty::DistanceConstraint>(p_2, p_3, cloth_bending_stiffness, (x_2 - x_3).norm()));

                    break;
                }
            }
        }

        const auto find_and_constrain_fixed_point = [&](const Eigen::Vector3d& search_position,
                                                        const Eigen::Vector3d& fixed_position,
                                                        const std::shared_ptr<elasty::Particle> particle)
        {
            if (particle->x.isApprox(search_position))
            {
                m_constraints.push_back(std::make_shared<elasty::FixedPointConstraint>(particle, 1.0, fixed_position));
            }
        };

        for (const auto& key_value : map_from_obj_vertex_index_to_particle)
        {
            const auto particle = key_value.second;

            find_and_constrain_fixed_point(Eigen::Vector3d(+ 1.0 + 1.0, + 2.0, 0.0), Eigen::Vector3d(+ 1.0 + 1.0, 3.0, 0.0), particle);
            find_and_constrain_fixed_point(Eigen::Vector3d(- 1.0 + 1.0, + 2.0, 0.0), Eigen::Vector3d(- 1.0 + 1.0, 3.0, 0.0), particle);
        }
    }

    TriangleIndices m_triangle_indices;
};

namespace
{
    inline glm::vec3 eigen2glm(const Eigen::Vector3d& eigen)
    {
        return glm::vec3(eigen.x(), eigen.y(), eigen.z());
    }
}

class SimpleEngine final : public elasty::Engine
{
public:

    void initializeScene() override
    {
        // Rod
        constexpr unsigned int num_particles = 20;
        constexpr double total_length = 2.0;
        constexpr double segment_length = total_length / double(num_particles - 1);

        std::shared_ptr<elasty::Particle> last_particle = nullptr;

        for (unsigned int i = 0; i < num_particles; ++ i)
        {
            const Eigen::Vector3d x = Eigen::Vector3d(- 1.0, 1.0 + segment_length * double(i), 0.0);
            const Eigen::Vector3d v = 50.0 * Eigen::Vector3d::Random();
            const double m = 1.0;

            auto particle = std::make_shared<elasty::Particle>(x, v, m);

            m_particles.push_back(particle);

            if (last_particle != nullptr)
            {
                addConstraint(std::make_shared<elasty::DistanceConstraint>(last_particle, particle, 0.5, segment_length));
            }

            last_particle = particle;
        }

        addConstraint(std::make_shared<elasty::FixedPointConstraint>(last_particle, 1.0, last_particle->x));

        // Cloth
        m_cloth_sim_object = std::make_shared<ClothSimObject>();
        std::copy(m_cloth_sim_object->m_particles.begin(), m_cloth_sim_object->m_particles.end(), std::back_inserter(m_particles));
        std::copy(m_cloth_sim_object->m_constraints.begin(), m_cloth_sim_object->m_constraints.end(), std::back_inserter(m_constraints));
    }

    void setExternalForces() override
    {
        const Eigen::Vector3d gravity = Eigen::Vector3d(0.0, - 9.8, 0.0);

        for (auto particle : m_particles)
        {
            particle->f = particle->m * gravity;
        }
    }

    void generateCollisionConstraints() override
    {
        for (const auto particle : m_particles)
        {
            if (particle->p.y() < 0.0)
            {
                addInstantConstraint(std::make_shared<elasty::EnvironmentalCollisionConstraint>(particle, 1.0, Eigen::Vector3d(0.0, 1.0, 0.0), 0.0));
            }
        }
    }

    void updateVelocities() override
    {
        for (auto particle : m_particles)
        {
            particle->v *= 0.999;
        }
    }

    std::shared_ptr<ClothSimObject> m_cloth_sim_object;
};

class SimpleApp final : public bigger::App
{
public:

    SimpleApp();

    void initialize(int argc, char** argv) override;
    void onReset() override;

    void updateApp() override;
    void releaseSharedResources() override;

    std::shared_ptr<SimpleEngine> m_engine;

    bool m_capture_screen;

private:

    // Shared resources
    std::shared_ptr<bigger::BlinnPhongMaterial> m_default_material;
    std::shared_ptr<bigger::BlinnPhongMaterial> m_checker_white_material;
    std::shared_ptr<bigger::BlinnPhongMaterial> m_checker_black_material;
    std::shared_ptr<bigger::SpherePrimitive> m_sphere_primitive;
    std::shared_ptr<bigger::PlanePrimitive> m_plane_primitive;
};

SimpleApp::SimpleApp()
{
    getCamera().m_position = glm::vec3(1.0f, 2.0f, - 8.0f);
    getCamera().m_target = glm::vec3(0.0f, 1.5f, 0.0f);
    getCamera().m_fov = 45.0f;

    m_capture_screen = false;
}

class ClothObject final : public bigger::SceneObject
{
public:

    ClothObject(std::shared_ptr<ClothSimObject> cloth_sim_object,
                std::shared_ptr<bigger::BlinnPhongMaterial> material) :
    bigger::SceneObject(material),
    m_cloth_sim_object(cloth_sim_object)
    {
    }

    void draw(const glm::mat4& parent_transform_matrix = glm::mat4(1.0f)) override
    {
        // TODO
    }

    std::shared_ptr<ClothSimObject> m_cloth_sim_object;
};

class ParticlesObject final : public bigger::SceneObject
{
public:

    ParticlesObject(std::shared_ptr<SimpleEngine> engine,
                    std::shared_ptr<bigger::SpherePrimitive> sphere_primitive,
                    std::shared_ptr<bigger::BlinnPhongMaterial> material) :
    bigger::SceneObject(material),
    m_engine(engine),
    m_sphere_primitive(sphere_primitive)
    {
    }

    void draw(const glm::mat4& parent_transform_matrix = glm::mat4(1.0f)) override
    {
        m_material->submitUniforms();

        for (auto& particle : m_engine->m_particles)
        {
            const glm::mat4 translate_matrix = glm::translate(eigen2glm(particle->x));
            const glm::mat4 scale_matrix = glm::scale(glm::vec3(0.05f));

            const glm::mat4 transform = parent_transform_matrix * translate_matrix * scale_matrix;

            bgfx::setTransform(glm::value_ptr(transform));

            m_sphere_primitive->submitPrimitive(m_material->m_program, true);
        }
    }

private:

    std::shared_ptr<SimpleEngine> m_engine;
    std::shared_ptr<bigger::SpherePrimitive> m_sphere_primitive;
};

class CheckerBoardObject final : public bigger::SceneObject
{
public:

    CheckerBoardObject(std::shared_ptr<bigger::PlanePrimitive> plane_primitive,
                       std::shared_ptr<bigger::BlinnPhongMaterial> checker_white_material,
                       std::shared_ptr<bigger::BlinnPhongMaterial> checker_black_material) :
    bigger::SceneObject(nullptr),
    m_plane_primitive(plane_primitive),
    m_checker_white_material(checker_white_material),
    m_checker_black_material(checker_black_material)
    {
    }

    void draw(const glm::mat4& parent_transform_matrix = glm::mat4(1.0f))
    {
        const glm::mat4 transform = parent_transform_matrix * getTransform();

        for (int i = - m_resolution; i <= m_resolution; ++ i)
        {
            for (int j = - m_resolution; j <= m_resolution; ++ j)
            {
                const glm::mat4 local_transform = glm::translate(transform, glm::vec3(float(i), 0.0f, float(j)));

                bgfx::setTransform(glm::value_ptr(local_transform));

                auto target_material = ((i + j) % 2 == 0) ? m_checker_white_material : m_checker_black_material;
                target_material->submitUniforms();
                m_plane_primitive->submitPrimitive(target_material->m_program);
            }
        }
    }

private:

    const int m_resolution = 4;

    std::shared_ptr<bigger::PlanePrimitive> m_plane_primitive;
    std::shared_ptr<bigger::BlinnPhongMaterial> m_checker_white_material;
    std::shared_ptr<bigger::BlinnPhongMaterial> m_checker_black_material;
};

void SimpleApp::initialize(int argc, char** argv)
{
    // Register and apply BGFX configuration settings
    reset(BGFX_RESET_VSYNC | BGFX_RESET_MSAA_X8);

    m_default_material = std::make_shared<bigger::BlinnPhongMaterial>();
    m_checker_white_material = std::make_shared<bigger::BlinnPhongMaterial>();
    m_checker_white_material->u_diffuse = glm::vec3(0.8);
    m_checker_white_material->u_specular = glm::vec3(0.0);
    m_checker_black_material = std::make_shared<bigger::BlinnPhongMaterial>();
    m_checker_black_material->u_diffuse = glm::vec3(0.3);
    m_checker_black_material->u_specular = glm::vec3(0.0);
    m_sphere_primitive = std::make_shared<bigger::SpherePrimitive>();
    m_sphere_primitive->initializePrimitive();
    m_plane_primitive = std::make_shared<bigger::PlanePrimitive>();
    m_plane_primitive->initializePrimitive();

    m_engine = std::make_unique<SimpleEngine>();
    m_engine->initializeScene();

    addSceneObject(std::make_shared<ParticlesObject>(m_engine, m_sphere_primitive, m_default_material));
    addSceneObject(std::make_shared<CheckerBoardObject>(m_plane_primitive, m_checker_white_material, m_checker_black_material));
    addSceneObject(std::make_shared<ClothObject>(m_engine->m_cloth_sim_object, m_default_material));
}

void SimpleApp::onReset()
{
    constexpr uint32_t bg_color = 0x303030ff;

    bgfx::setViewClear(0, BGFX_CLEAR_COLOR | BGFX_CLEAR_DEPTH, bg_color, 1.0f, 0);
}

void SimpleApp::updateApp()
{
    // Display ImGui components
    ImGui::Begin("Config");
    {
        ImGui::Text("frame: %d", m_frame);
        ImGui::Text("time: %.2f", m_time);
        ImGui::Text("fps: %.2f", 1.0f / m_last_dt);
        ImGui::Separator();
        ImGui::SliderFloat3("camera.position", glm::value_ptr(getCamera().m_position), - 10.0f, 10.0f);
        ImGui::SliderFloat("camera.fov", &(getCamera().m_fov), 10.0f, 120.0f);
        ImGui::Separator();
        m_default_material->drawImgui();
    }
    ImGui::End();

    ImGui::Begin("Control");
    {
        if (ImGui::Button("Reset"))
        {
            m_engine->clearScene();
            m_engine->initializeScene();
        }
    }
    ImGui::End();

    // Request screen capture
    if (m_capture_screen)
    {
        const std::string dir_path = ".";
        const std::string file_name = stringutil::ConvertWithZeroPadding(m_frame, 5);
        const std::string file_path = dir_path + "/" + file_name;
        bgfx::requestScreenShot(BGFX_INVALID_HANDLE, file_path.c_str());
    }

    // Physics
    m_engine->stepTime();
}

void SimpleApp::releaseSharedResources()
{
    m_sphere_primitive = nullptr;
    m_default_material = nullptr;
    m_checker_white_material = nullptr;
    m_checker_black_material = nullptr;
}

int main(int argc, char** argv)
{
    SimpleApp app;
#if 0
    app.m_capture_screen = true;
    return app.runApp(argc, argv, bgfx::RendererType::OpenGL);
#else
    return app.runApp(argc, argv);
#endif
}
