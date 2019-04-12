#include <array>
#include <memory>
#include <bigger/app.hpp>
#include <bigger/scene-object.hpp>
#include <bigger/materials/blinnphong-material.hpp>
#include <bigger/primitives/sphere-primitive.hpp>
#include <elasty/elasty.hpp>

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
        m_particles.resize(1);

        m_particles[0].x = Eigen::Vector3d::Zero();
        m_particles[0].v = Eigen::Vector3d(0.0, 10.0, 0.0);
        m_particles[0].m = 1.0;
    }

    void stepTime() override
    {
        elasty::Particle& particle = m_particles[0];

        constexpr double dt = 1.0 / 60.0;
        const Eigen::Vector3d gravity = Eigen::Vector3d(0.0, - 9.8, 0.0);

        // Apply external forces
        const Eigen::Vector3d external_forces = particle.m * gravity;
        particle.v = particle.v + dt * (1.0 / particle.m) * external_forces;

        // Calculate predicted positions
        particle.p = particle.x + dt * particle.v;

        // Generate collision constraints
        if (particle.p.y() < 0.0)
        {
            addInstantConstraint(std::make_shared<elasty::EnvironmentalCollisionConstraint>(this, std::vector<unsigned int>{ 0 }, Eigen::Vector3d(0.0, 1.0, 0.0), 0.0));
        }

        // Solve constraints
        constexpr unsigned int num_iterations = 10;
        for (unsigned int i = 0; i < num_iterations; ++ i)
        {
            for (auto constraint : m_constraints)
            {
                projectConstraint(constraint);
            }

            for (auto constraint : m_instant_constraints)
            {
                projectConstraint(constraint);
            }
        }

        // Apply the results
        particle.v = (particle.p - particle.x) * (1.0 / dt);
        particle.x = particle.p;

        // Update velocities
        particle.v *= 0.999;

        // Clear instant constraints
        m_instant_constraints.clear();
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
    std::shared_ptr<bigger::SpherePrimitive> m_sphere_primitive;
};

SimpleApp::SimpleApp()
{
    getCamera().m_position = glm::vec3(1.0f, 1.0f, - 10.0f);
}

class ParticleObject final : public bigger::SceneObject
{
public:
    ParticleObject(std::shared_ptr<SimpleEngine> engine,
                   std::shared_ptr<bigger::SpherePrimitive> sphere_primitive,
                   std::shared_ptr<bigger::BlinnPhongMaterial> material) :
    bigger::SceneObject(material),
    m_engine(engine),
    m_sphere_primitive(sphere_primitive)
    {
    }

    void draw(const glm::mat4& parent_transform_matrix = glm::mat4(1.0f))
    {
        const glm::mat4 translate_matrix = glm::translate(eigen2glm(m_engine->m_particles[0].x));
        const glm::mat4 scale_matrix = glm::scale(glm::vec3(0.1f));

        const glm::mat4 transform = parent_transform_matrix * translate_matrix * scale_matrix;

        bgfx::setTransform(glm::value_ptr(transform));

        m_sphere_primitive->submitPrimitive(m_material->m_program);
    }

private:

    std::shared_ptr<SimpleEngine> m_engine;
    std::shared_ptr<bigger::SpherePrimitive> m_sphere_primitive;
};

void SimpleApp::initialize(int argc, char** argv)
{
    // Register and apply BGFX configuration settings
    reset(BGFX_RESET_VSYNC | BGFX_RESET_MSAA_X8);

    m_default_material = std::make_shared<bigger::BlinnPhongMaterial>();
    m_sphere_primitive = std::make_shared<bigger::SpherePrimitive>();
    m_sphere_primitive->initializePrimitive();

    m_engine = std::make_unique<SimpleEngine>();
    m_engine->initializeScene();

    addSceneObject(std::make_shared<ParticleObject>(m_engine, m_sphere_primitive, m_default_material));
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
}

int main(int argc, char** argv)
{
    SimpleApp app;
    return app.run(argc, argv);
}
