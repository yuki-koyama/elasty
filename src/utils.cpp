#include <Alembic/AbcCoreOgawa/All.h>
#include <Alembic/AbcGeom/All.h>
#include <elasty/cloth-sim-object.hpp>
#include <elasty/constraint.hpp>
#include <elasty/particle.hpp>
#include <elasty/utils.hpp>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <utility>

std::vector<float> elasty::packParticlePositions(const std::vector<std::shared_ptr<Particle>>& particles)
{
    std::vector<float> verts;
    verts.reserve(3 * particles.size());
    for (const auto& particle : particles)
    {
        verts.push_back(particle->x(0));
        verts.push_back(particle->x(1));
        verts.push_back(particle->x(2));
    }
    return verts;
}

void elasty::setRandomVelocities(const std::vector<std::shared_ptr<Particle>>& particles, const double scale)
{
    for (const auto& particle : particles)
    {
        particle->v = scale * Eigen::Vector3d::Random();
    }
}

void elasty::generateFixedPointConstraints(const Eigen::Vector3d&                            search_position,
                                           const Eigen::Vector3d&                            fixed_position,
                                           const std::vector<std::shared_ptr<Particle>>&     particles,
                                           std::vector<std::shared_ptr<AbstractConstraint>>& constraints)
{
    constexpr double dummy_delta_time = 1.0 / 60.0;

    for (const auto& particle : particles)
    {
        if (particle->x.isApprox(search_position))
        {
            constraints.push_back(
                std::make_shared<elasty::FixedPointConstraint>(particle, 1.0, 0.0, dummy_delta_time, fixed_position));
        }
    }
}

class elasty::AlembicManager
{
public:
    AlembicManager(const std::string&                    output_file_path,
                   const std::shared_ptr<ClothSimObject> cloth_sim_object,
                   const double                          delta_time = 1.0 / 60.0)
        : m_cloth_sim_object(cloth_sim_object),
          m_archive(Alembic::AbcCoreOgawa::WriteArchive(), output_file_path.c_str())
    {
        using namespace Alembic::Abc;
        using namespace Alembic::AbcGeom;

        const TimeSampling time_sampling(delta_time, 0);
        const uint32_t     time_sampling_index = m_archive.addTimeSampling(time_sampling);

        m_mesh_obj = OPolyMesh(OObject(m_archive, kTop), "cloth");
        m_mesh_obj.getSchema().setTimeSampling(time_sampling_index);
    }

    void submitCurrentStatus()
    {
        using namespace Alembic::Abc;
        using namespace Alembic::AbcGeom;

        const size_t             num_verts = m_cloth_sim_object->m_particles.size();
        const std::vector<float> verts     = packParticlePositions(m_cloth_sim_object->m_particles);

        // If this is the first call, set a sample with full properties including vertex positions, indices, counts, and
        // UVs (if exists); otherwise, set a sample with only vertex positions.
        if (m_is_first)
        {
            const size_t               num_indices = m_cloth_sim_object->m_triangle_list.size();
            const std::vector<int32_t> indices     = [&]() {
                std::vector<int32_t> indices;
                indices.reserve(num_indices);
                for (unsigned int i = 0; i < m_cloth_sim_object->m_triangle_list.rows(); ++i)
                {
                    indices.push_back(m_cloth_sim_object->m_triangle_list(i, 0));
                    indices.push_back(m_cloth_sim_object->m_triangle_list(i, 1));
                    indices.push_back(m_cloth_sim_object->m_triangle_list(i, 2));
                }
                return indices;
            }();
            const size_t               num_counts = m_cloth_sim_object->m_triangle_list.rows();
            const std::vector<int32_t> counts(num_counts, 3);

            // If the model has UV info, include it into the geometry sample
            const OV2fGeomParam::Sample geom_param_sample =
                m_cloth_sim_object->hasUv()
                    ? OV2fGeomParam::Sample(
                          V2fArraySample((const V2f*) m_cloth_sim_object->m_uv_list.data(), num_indices),
                          kFacevaryingScope)
                    : OV2fGeomParam::Sample();

            const OPolyMeshSchema::Sample sample(V3fArraySample((const V3f*) verts.data(), num_verts),
                                                 Int32ArraySample(indices.data(), num_indices),
                                                 Int32ArraySample(counts.data(), num_counts),
                                                 geom_param_sample);

            m_mesh_obj.getSchema().set(sample);

            m_is_first = false;
        }
        else
        {
            const OPolyMeshSchema::Sample sample(V3fArraySample((const V3f*) verts.data(), num_verts));
            m_mesh_obj.getSchema().set(sample);
        }
    }

private:
    bool                                  m_is_first = true;
    const std::shared_ptr<ClothSimObject> m_cloth_sim_object;
    Alembic::Abc::OArchive                m_archive;
    Alembic::AbcGeom::OPolyMesh           m_mesh_obj;
};

std::shared_ptr<elasty::AlembicManager>
elasty::createAlembicManager(const std::string&                    file_path,
                             const std::shared_ptr<ClothSimObject> cloth_sim_object,
                             const double                          delta_time)
{
    return std::make_shared<AlembicManager>(file_path, cloth_sim_object, delta_time);
}

void elasty::submitCurrentStatus(const std::shared_ptr<AlembicManager> alembic_manager)
{
    alembic_manager->submitCurrentStatus();
}

void elasty::exportCurrentClothStateAsObj(const std::string&                    output_file_path,
                                          const std::shared_ptr<ClothSimObject> cloth_sim_object)
{
    const auto& particles = cloth_sim_object->m_particles;
    const auto& triangles = cloth_sim_object->m_triangle_list;
    const auto& uvs       = cloth_sim_object->m_uv_list;

    // Convert to the OBJ format
    std::ofstream file(output_file_path);
    for (const auto& particle : particles)
    {
        file << "v " << particle->x(0) << " " << particle->x(1) << " " << particle->x(2) << std::endl;
    }
    for (unsigned row = 0; row < uvs.rows(); ++row)
    {
        file << "vt " << uvs.row(row)(0) << " " << uvs.row(row)(1) << std::endl;
    }
    for (unsigned row = 0; row < triangles.rows(); ++row)
    {
        file << "f ";
        file << triangles.row(row)(0) + 1 << "/" << triangles.row(row)(0) + 1 << " ";
        file << triangles.row(row)(1) + 1 << "/" << triangles.row(row)(1) + 1 << " ";
        file << triangles.row(row)(2) + 1 << "/" << triangles.row(row)(2) + 1 << std::endl;
    }
    file.close();
}

std::string elasty::generateClothMeshObjData(const double   width,
                                             const double   height,
                                             const unsigned horizontal_resolution,
                                             const unsigned vertical_resolution)
{
    std::vector<Eigen::Vector3d> vertices;
    std::vector<Eigen::Vector3i> triangles;
    std::vector<Eigen::Vector2d> uvs;

    // Vertices
    for (unsigned v_index = 0; v_index <= vertical_resolution; ++v_index)
    {
        for (unsigned h_index = 0; h_index <= horizontal_resolution; ++h_index)
        {
            const double u = [&]() {
                if (v_index % 2 == 0 || h_index == 0)
                {
                    return static_cast<double>(h_index) / static_cast<double>(horizontal_resolution);
                }
                else
                {
                    return (static_cast<double>(h_index) - 0.5) / static_cast<double>(horizontal_resolution);
                }
            }();
            const double v = static_cast<double>(v_index) / static_cast<double>(vertical_resolution);

            const double x = (u - 0.5) * width;
            const double y = (v - 0.5) * height;

            vertices.push_back(Eigen::Vector3d(x, 0.0, y));
            uvs.push_back(Eigen::Vector2d(u, v));

            // Additional vetex at the even-indexed row
            if (v_index % 2 == 1 && h_index == horizontal_resolution)
            {
                vertices.push_back(Eigen::Vector3d(0.5 * width, 0.0, y));
                uvs.push_back(Eigen::Vector2d(1.0, v));
            }
        }
    }

    // Triangles
    for (unsigned v_index = 0; v_index < vertical_resolution; ++v_index)
    {
        if (v_index % 2 == 0)
        {
            const unsigned top_row_begin    = (2 * (horizontal_resolution + 1) + 1) * (v_index / 2);
            const unsigned bottom_row_begin = top_row_begin + horizontal_resolution + 1;

            for (unsigned h_index = 0; h_index <= horizontal_resolution; ++h_index)
            {
                if (h_index == 0)
                {
                    triangles.push_back(
                        Eigen::Vector3i(top_row_begin + h_index, bottom_row_begin + 0, bottom_row_begin + 1));
                }
                else
                {
                    triangles.push_back(Eigen::Vector3i(
                        top_row_begin + h_index, top_row_begin + h_index - 1, bottom_row_begin + h_index));
                    triangles.push_back(Eigen::Vector3i(
                        top_row_begin + h_index, bottom_row_begin + h_index, bottom_row_begin + h_index + 1));
                }
            }
        }
        else
        {
            const unsigned top_row_begin =
                (2 * (horizontal_resolution + 1) + 1) * ((v_index - 1) / 2) + horizontal_resolution + 1;
            const unsigned bottom_row_begin = top_row_begin + horizontal_resolution + 2;

            for (unsigned h_index = 0; h_index <= horizontal_resolution; ++h_index)
            {
                if (h_index == 0)
                {
                    triangles.push_back(Eigen::Vector3i(
                        top_row_begin + h_index, bottom_row_begin + h_index, top_row_begin + h_index + 1));
                }
                else
                {
                    triangles.push_back(Eigen::Vector3i(
                        top_row_begin + h_index, bottom_row_begin + h_index - 1, bottom_row_begin + h_index));
                    triangles.push_back(Eigen::Vector3i(
                        top_row_begin + h_index, bottom_row_begin + h_index, top_row_begin + h_index + 1));
                }
            }
        }
    }

    // Convert to the OBJ format
    std::stringstream sstream;
    for (const auto& vertex : vertices)
    {
        sstream << "v " << vertex(0) << " " << vertex(1) << " " << vertex(2) << std::endl;
    }
    for (const auto& uv : uvs)
    {
        sstream << "vt " << uv(0) << " " << uv(1) << std::endl;
    }
    for (const auto& triangle : triangles)
    {
        sstream << "f ";
        sstream << triangle(0) + 1 << "/" << triangle(0) + 1 << " ";
        sstream << triangle(1) + 1 << "/" << triangle(1) + 1 << " ";
        sstream << triangle(2) + 1 << "/" << triangle(2) + 1 << std::endl;
    }

    return sstream.str();
}
