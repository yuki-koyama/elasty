#include <array>
#include <memory>
#include <bigger/app.hpp>
#include <bigger/scene-object.hpp>
#include <bigger/materials/blinnphong-material.hpp>
#include <bigger/primitives/plane-primitive.hpp>
#include <bigger/primitives/sphere-primitive.hpp>
#include <elasty/constraint.hpp>
#include <elasty/engine.hpp>

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
        constexpr unsigned int num_particles = 20;
        constexpr double total_length = 3.0;
        constexpr double segment_length = total_length / double(num_particles - 1);

        m_particles.resize(num_particles);

        for (unsigned int i = 0; i < num_particles; ++ i)
        {
            m_particles[i].x = Eigen::Vector3d(0.0, 1.0 + segment_length * double(i), 0.0);
            m_particles[i].v = 50.0 * Eigen::Vector3d::Random();
            m_particles[i].m = 1.0;
            m_particles[i].i = i;
        }

        addConstraint(std::make_shared<elasty::FixedPointConstraint>(this, std::vector<unsigned int>{ num_particles - 1 }, 1.0, m_particles[num_particles - 1].x));

        for (unsigned int i = 0; i < num_particles - 1; ++ i)
        {
            addConstraint(std::make_shared<elasty::DistanceConstraint>(this, std::vector<unsigned int>{ i, i + 1 }, 0.5, segment_length));
        }
    }

    void setExternalForces() override
    {
        const Eigen::Vector3d gravity = Eigen::Vector3d(0.0, - 9.8, 0.0);

        for (auto& particle : m_particles)
        {
            particle.f = particle.m * gravity;
        }
    }

    void generateCollisionConstraints() override
    {
        for (auto& particle : m_particles)
        {
            if (particle.p.y() < 0.0)
            {
                addInstantConstraint(std::make_shared<elasty::EnvironmentalCollisionConstraint>(this, std::vector<unsigned int>{ particle.i }, 1.0, Eigen::Vector3d(0.0, 1.0, 0.0), 0.0));
            }
        }
    }

    void updateVelocities() override
    {
        for (auto& particle : m_particles)
        {
            particle.v *= 0.999;
        }
    }
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
    getCamera().m_position = glm::vec3(1.0f, 2.0f, - 10.0f);
    getCamera().m_target = glm::vec3(0.0f, 1.0f, 0.0f);
}

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

    void draw(const glm::mat4& parent_transform_matrix = glm::mat4(1.0f))
    {
        for (auto& particle : m_engine->m_particles)
        {
            const glm::mat4 translate_matrix = glm::translate(eigen2glm(particle.x));
            const glm::mat4 scale_matrix = glm::scale(glm::vec3(0.1f));

            const glm::mat4 transform = parent_transform_matrix * translate_matrix * scale_matrix;

            bgfx::setTransform(glm::value_ptr(transform));

            m_material->submitUniforms();
            m_sphere_primitive->submitPrimitive(m_material->m_program);
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
        ImGui::Text("time: %.2f", m_time);
        ImGui::Text("fps: %.2f", 1.0f / m_last_dt);
        ImGui::Separator();
        ImGui::SliderFloat3("camera.position", glm::value_ptr(getCamera().m_position), - 10.0f, 10.0f);
        ImGui::SliderFloat("camera.fov", &(getCamera().m_fov), 10.0f, 120.0f);
        ImGui::Separator();
        m_default_material->drawImgui();
    }
    ImGui::End();

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
    return app.run(argc, argv);
}
