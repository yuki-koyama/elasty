#include <array>
#include <bigger/app.hpp>
#include <bigger/materials/blinnphong-material.hpp>
#include <bigger/primitives/dynamic-mesh-primitive.hpp>
#include <bigger/primitives/plane-primitive.hpp>
#include <bigger/primitives/sphere-primitive.hpp>
#include <bigger/scene-object.hpp>
#include <elasty/cloth-sim-object.hpp>
#include <elasty/constraint.hpp>
#include <elasty/engine.hpp>
#include <elasty/particle.hpp>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <string-util.hpp>

// #define CAPTURE_SCREEN
// #define EXPORT_ALEMBIC

#ifdef EXPORT_ALEMBIC
#include <elasty/alembic-manager.hpp>
#endif

namespace
{
    inline glm::vec3 eigen2glm(const Eigen::Vector3d& eigen)
    {
        return glm::vec3(eigen.x(), eigen.y(), eigen.z());
    }
} // namespace

class SimpleEngine final : public elasty::AbstractEngine
{
public:
    void initializeScene() override
    {
        // Rod
        constexpr unsigned int num_particles  = 30;
        constexpr double       total_length   = 2.0;
        constexpr double       segment_length = total_length / double(num_particles - 1);

        std::shared_ptr<elasty::Particle> last_particle = nullptr;

        for (unsigned int i = 0; i < num_particles; ++i)
        {
            const Eigen::Vector3d x = Eigen::Vector3d(-1.0, 1.0 + segment_length * double(i), 0.0);
            const Eigen::Vector3d v = 50.0 * Eigen::Vector3d::Random();
            const double          m = 1.0;

            auto particle = std::make_shared<elasty::Particle>(x, v, m);

            m_particles.push_back(particle);

            if (last_particle != nullptr)
            {
                addConstraint(std::make_shared<elasty::DistanceConstraint>(
                    last_particle, particle, 0.5, 0.0, getDeltaPhysicsTime(), segment_length));
            }

            last_particle = particle;
        }

        // Pin the tip of the rod
        addConstraint(std::make_shared<elasty::FixedPointConstraint>(
            last_particle, 1.0, 0.0, getDeltaPhysicsTime(), last_particle->x));

        // Register the cloth object
        std::copy(m_cloth_sim_object->m_particles.begin(),
                  m_cloth_sim_object->m_particles.end(),
                  std::back_inserter(m_particles));
        std::copy(m_cloth_sim_object->m_constraints.begin(),
                  m_cloth_sim_object->m_constraints.end(),
                  std::back_inserter(m_constraints));

        // Pin two of the corners of the cloth
        constexpr double range_radius = 0.1;
        for (const auto& particle : m_particles)
        {
            const Eigen::Vector3d translation{0.0, 0.0, 0.0};

            if ((particle->x - Eigen::Vector3d(+1.0 + 1.0, 3.0, 0.0)).norm() < range_radius)
            {
                m_constraints.push_back(std::make_shared<elasty::FixedPointConstraint>(
                    particle, 1.0, 0.0, getDeltaPhysicsTime(), particle->x + translation));
            }
            if ((particle->x - Eigen::Vector3d(-1.0 + 1.0, 3.0, 0.0)).norm() < range_radius)
            {
                m_constraints.push_back(std::make_shared<elasty::FixedPointConstraint>(
                    particle, 1.0, 0.0, getDeltaPhysicsTime(), particle->x + translation));
            }
        }

        // Create a solid blob object
        constexpr double shape_matching_stiffness     = 0.2;
        constexpr double shape_matching_total_mass    = 5.0;
        constexpr int    shape_matching_num_particles = 1000;

        std::vector<std::shared_ptr<elasty::Particle>> shape_matching_particles;
        for (int i = 0; i < shape_matching_num_particles; ++i)
        {
            const double          size     = 0.4;
            const Eigen::Vector3d rand_vec = Eigen::Vector3d::Random();
            const Eigen::Vector3d offset   = Eigen::Vector3d(0.0, 5.0, -2.0);

            const Eigen::Vector3d x = size * (rand_vec / rand_vec.lpNorm<Eigen::Infinity>()) + offset;
            const Eigen::Vector3d v = Eigen::Vector3d::Random();
            const double          m = shape_matching_total_mass / shape_matching_num_particles;

            auto particle = std::make_shared<elasty::Particle>(x, v, m);

            shape_matching_particles.push_back(particle);
            m_particles.push_back(particle);
        }
        auto shape_matching_constaint = std::make_shared<elasty::ShapeMatchingConstraint>(
            shape_matching_particles, shape_matching_stiffness, 0.0, getDeltaPhysicsTime());
        m_constraints.push_back(shape_matching_constaint);
    }

    void setExternalForces() override
    {
        const Eigen::Vector3d gravity = Eigen::Vector3d(0.0, -9.8, 0.0);

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
                addInstantConstraint(std::make_shared<elasty::EnvironmentalCollisionConstraint>(
                    particle, 1.0, 0.0, getDeltaPhysicsTime(), Eigen::Vector3d(0.0, 1.0, 0.0), 0.0));
            }
        }
    }

    void updateVelocities() override
    {
        const double decay_rate = std::exp(std::log(0.9) * getDeltaPhysicsTime());

        for (auto particle : m_particles)
        {
            particle->v *= decay_rate;
        }
    }

    std::shared_ptr<elasty::ClothSimObject> m_cloth_sim_object;
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
    const double          m_cloth_in_plane_stiffness      = 0.95;
    const double          m_cloth_in_plane_compliance     = 5e-2;
    const double          m_cloth_out_of_plane_stiffness  = 0.05;
    const double          m_cloth_out_of_plane_compliance = 1e+4;
    const unsigned        m_cloth_resolution              = 20;
    const Eigen::Affine3d m_cloth_import_transform        = Eigen::Affine3d(Eigen::Translation3d(1.0, 3.0, 1.0));

    // Shared resources
    std::shared_ptr<bigger::BlinnPhongMaterial> m_default_material;
    std::shared_ptr<bigger::BlinnPhongMaterial> m_checker_white_material;
    std::shared_ptr<bigger::BlinnPhongMaterial> m_checker_black_material;
    std::shared_ptr<bigger::SpherePrimitive>    m_sphere_primitive;
    std::shared_ptr<bigger::PlanePrimitive>     m_plane_primitive;
};

SimpleApp::SimpleApp()
{
    getCamera().m_position = glm::vec3(1.0f, 2.0f, -8.0f);
    getCamera().m_target   = glm::vec3(0.0f, 1.5f, 0.0f);
    getCamera().m_fov      = 45.0f;

    m_capture_screen = false;
}

class ClothObject final : public bigger::SceneObject
{
public:
    ClothObject(std::shared_ptr<elasty::ClothSimObject>     cloth_sim_object,
                std::shared_ptr<bigger::BlinnPhongMaterial> material)
        : bigger::SceneObject(material), m_cloth_sim_object(cloth_sim_object)
    {
        std::vector<bigger::PositionNormalUvVertex> vertex_data = generateVertexData();

        std::vector<uint16_t> triangle_list;
        for (unsigned int i = 0; i < cloth_sim_object->m_triangle_list.rows(); ++i)
        {
            triangle_list.push_back(cloth_sim_object->m_triangle_list(i, 0));
            triangle_list.push_back(cloth_sim_object->m_triangle_list(i, 1));
            triangle_list.push_back(cloth_sim_object->m_triangle_list(i, 2));
        }

        m_dynamic_mesh_primitive = std::make_unique<bigger::DynamicMeshPrimitive>(vertex_data, triangle_list);

#ifdef EXPORT_ALEMBIC
        m_alembic_manager = elasty::createClothAlembicManager("./cloth.abc", m_cloth_sim_object, 1.0 / 60.0);
        m_alembic_manager->submitCurrentStatus();
#endif
    }

    void update() override
    {
#ifdef EXPORT_ALEMBIC
        assert(m_alembic_manager != nullptr);
        m_alembic_manager->submitCurrentStatus();
#endif
    }

    void draw(const glm::mat4& parent_transform_matrix = glm::mat4(1.0f)) override
    {
        m_material->submitUniforms();

        const std::vector<bigger::PositionNormalUvVertex> vertex_data = generateVertexData();
        m_dynamic_mesh_primitive->updateVertexData(vertex_data);

        // Do not cull back-facing triangles
        bgfx::setState(BGFX_STATE_DEFAULT & (~BGFX_STATE_CULL_CW));

        m_dynamic_mesh_primitive->submitPrimitive(m_material->m_program);
    }

    std::shared_ptr<elasty::ClothSimObject> m_cloth_sim_object;

private:
#ifdef EXPORT_ALEMBIC
    std::shared_ptr<elasty::AlembicManager> m_alembic_manager;
#endif

    std::unique_ptr<bigger::DynamicMeshPrimitive> m_dynamic_mesh_primitive;

    std::vector<bigger::PositionNormalUvVertex> generateVertexData() const
    {
        std::vector<bigger::PositionNormalUvVertex> vertex_data(m_cloth_sim_object->m_particles.size());
        for (unsigned int i = 0; i < m_cloth_sim_object->m_particles.size(); ++i)
        {
            vertex_data[i] = {{m_cloth_sim_object->m_particles[i]->x(0),
                               m_cloth_sim_object->m_particles[i]->x(1),
                               m_cloth_sim_object->m_particles[i]->x(2)},
                              glm::vec3(0.0f)};
        }

        for (unsigned int i = 0; i < m_cloth_sim_object->m_triangle_list.rows(); ++i)
        {
            const auto& i_0 = m_cloth_sim_object->m_triangle_list(i, 0);
            const auto& i_1 = m_cloth_sim_object->m_triangle_list(i, 1);
            const auto& i_2 = m_cloth_sim_object->m_triangle_list(i, 2);

            const glm::vec3& x_0 = vertex_data[i_0].position;
            const glm::vec3& x_1 = vertex_data[i_1].position;
            const glm::vec3& x_2 = vertex_data[i_2].position;

            const glm::vec3 area_scaled_face_normal = glm::cross(x_1 - x_0, x_2 - x_0);

            vertex_data[i_0].normal += area_scaled_face_normal;
            vertex_data[i_1].normal += area_scaled_face_normal;
            vertex_data[i_2].normal += area_scaled_face_normal;
        }

        for (unsigned int i = 0; i < m_cloth_sim_object->m_particles.size(); ++i)
        {
            vertex_data[i].normal = glm::normalize(vertex_data[i].normal);
        }

        return vertex_data;
    }
};

class ParticlesObject final : public bigger::SceneObject
{
public:
    ParticlesObject(std::shared_ptr<SimpleEngine>               engine,
                    std::shared_ptr<bigger::SpherePrimitive>    sphere_primitive,
                    std::shared_ptr<bigger::BlinnPhongMaterial> material)
        : bigger::SceneObject(material), m_engine(engine), m_sphere_primitive(sphere_primitive)
    {
    }

    void draw(const glm::mat4& parent_transform_matrix = glm::mat4(1.0f)) override
    {
        m_material->submitUniforms();

        for (const auto& particle : m_engine->getParticles())
        {
            const glm::mat4 translate_matrix = glm::translate(glm::mat4(1.0f), eigen2glm(particle->x));
            const glm::mat4 scale_matrix     = glm::scale(glm::mat4(1.0f), glm::vec3(scale));

            const glm::mat4 transform = parent_transform_matrix * translate_matrix * scale_matrix;

            bgfx::setTransform(glm::value_ptr(transform));

            m_sphere_primitive->submitPrimitive(m_material->m_program, true);
        }
    }

private:
    const float scale = 0.02f;

    std::shared_ptr<SimpleEngine>            m_engine;
    std::shared_ptr<bigger::SpherePrimitive> m_sphere_primitive;
};

class CheckerBoardObject final : public bigger::SceneObject
{
public:
    CheckerBoardObject(std::shared_ptr<bigger::PlanePrimitive>     plane_primitive,
                       std::shared_ptr<bigger::BlinnPhongMaterial> checker_white_material,
                       std::shared_ptr<bigger::BlinnPhongMaterial> checker_black_material)
        : bigger::SceneObject(nullptr),
          m_plane_primitive(plane_primitive),
          m_checker_white_material(checker_white_material),
          m_checker_black_material(checker_black_material)
    {
    }

    void draw(const glm::mat4& parent_transform_matrix = glm::mat4(1.0f))
    {
        const glm::mat4 transform = parent_transform_matrix * getTransform();

        for (int i = -m_resolution; i <= m_resolution; ++i)
        {
            for (int j = -m_resolution; j <= m_resolution; ++j)
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

    std::shared_ptr<bigger::PlanePrimitive>     m_plane_primitive;
    std::shared_ptr<bigger::BlinnPhongMaterial> m_checker_white_material;
    std::shared_ptr<bigger::BlinnPhongMaterial> m_checker_black_material;
};

void SimpleApp::initialize(int argc, char** argv)
{
    // Register and apply BGFX configuration settings
    reset(BGFX_RESET_VSYNC | BGFX_RESET_MSAA_X8);

    m_default_material                   = std::make_shared<bigger::BlinnPhongMaterial>();
    m_checker_white_material             = std::make_shared<bigger::BlinnPhongMaterial>();
    m_checker_white_material->m_diffuse  = glm::vec3(0.9);
    m_checker_white_material->m_specular = glm::vec3(0.0);
    m_checker_black_material             = std::make_shared<bigger::BlinnPhongMaterial>();
    m_checker_black_material->m_diffuse  = glm::vec3(0.6);
    m_checker_black_material->m_specular = glm::vec3(0.0);
    m_sphere_primitive                   = std::make_shared<bigger::SpherePrimitive>();
    m_plane_primitive                    = std::make_shared<bigger::PlanePrimitive>();

    m_engine                     = std::make_unique<SimpleEngine>();
    m_engine->m_cloth_sim_object = std::make_shared<elasty::ClothSimObject>(m_cloth_resolution,
                                                                            m_cloth_in_plane_stiffness,
                                                                            m_cloth_in_plane_compliance,
                                                                            m_cloth_out_of_plane_stiffness,
                                                                            m_cloth_out_of_plane_compliance,
                                                                            m_engine->getDeltaPhysicsTime(),
                                                                            m_cloth_import_transform);
    m_engine->initializeScene();

    addSceneObject(std::make_shared<ParticlesObject>(m_engine, m_sphere_primitive, m_default_material));
    addSceneObject(
        std::make_shared<CheckerBoardObject>(m_plane_primitive, m_checker_white_material, m_checker_black_material));
    addSceneObject(std::make_shared<ClothObject>(m_engine->m_cloth_sim_object, m_default_material), "cloth");
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
        getCamera().drawImgui();
        ImGui::Separator();
        m_default_material->drawImgui();
    }
    ImGui::End();

    ImGui::Begin("Control");
    {
        if (ImGui::Button("Reset"))
        {
            // Clear
            m_engine->m_cloth_sim_object = nullptr;
            m_engine->clearScene();

            // Init
            m_engine->m_cloth_sim_object = std::make_shared<elasty::ClothSimObject>(m_cloth_resolution,
                                                                                    m_cloth_in_plane_stiffness,
                                                                                    m_cloth_in_plane_compliance,
                                                                                    m_cloth_out_of_plane_stiffness,
                                                                                    m_cloth_out_of_plane_compliance,
                                                                                    m_engine->getDeltaPhysicsTime(),
                                                                                    m_cloth_import_transform);
            m_engine->initializeScene();

            // Re-register
            m_scene_objects["cloth"] = std::make_shared<ClothObject>(m_engine->m_cloth_sim_object, m_default_material);
        }
    }
    ImGui::End();

    // Request screen capture
    if (m_capture_screen)
    {
        const std::string dir_path  = ".";
        const std::string file_name = stringutil::ConvertWithZeroPadding(m_frame, 5);
        const std::string file_path = dir_path + "/" + file_name;
        bgfx::requestScreenShot(BGFX_INVALID_HANDLE, file_path.c_str());
    }

    // Physics
    m_engine->proceedFrame();
}

void SimpleApp::releaseSharedResources()
{
    m_plane_primitive        = nullptr;
    m_sphere_primitive       = nullptr;
    m_default_material       = nullptr;
    m_checker_white_material = nullptr;
    m_checker_black_material = nullptr;
}

int main(int argc, char** argv)
{
    SimpleApp app;
#ifdef CAPTURE_SCREEN
    app.m_capture_screen = true;
    return app.runApp(argc, argv, bgfx::RendererType::OpenGL);
#else
    return app.runApp(argc, argv);
#endif
}
