#pragma once

#include <Integrators/Integrator.h>
#include <Integrators/DirectIntegrator.h>
#include <Math/Vectorf.h>
#include <mutex>
#include <tbb/mutex.h>

namespace Integrator
{
	class PhotonMapper : public DirectIntegrator
	{
	protected:

		using Vector3i = Math::Vector<int, 3>;
		using Vector3f = Math::Vector3f;
		using Vector2f = Math::Vector2f;

		template <class T>
		using Collection = std::vector<T>;

		using PhotonFloat = float;

		template <class Float>
		class Photon
		{
		public:
			const Primitive* m_primitive;
			Math::Vector<Float, 3> m_point;
			uint8_t m_depth;
			Math::Vector<Float, 3> m_dir;
			Float m_beta[3];

			Photon(Hit const& hit, uint8_t depth, RGBColor const& beta):
				m_primitive(hit.primitve),
				m_point(hit.point),
				m_depth(depth),
				m_dir(hit.to_view)
			{
				m_beta[0] = beta[0];
				m_beta[1] = beta[1];
				m_beta[2] = beta[2];
			}

			void fillHit(Hit& res)const
			{
				res.primitve = m_primitive;
				res.geometry = m_primitive->geometry();
				res.point = m_point;
				res.to_view = m_dir;

				res.primitive_uv = m_primitive->uv(m_point);
				res.tex_uv = m_primitive->tuv(res.primitive_uv);

				res.primitive_normal = m_primitive->normal(m_point, res.primitive_uv);
				res.normal = m_primitive->shading_normal(m_point, res.primitive_uv);

				res.z = -1;
				res.facing = res.normal * res.to_view > 0;
			}

			RGBColor beta()const
			{
				return RGBColor(m_beta[0], m_beta[1], m_beta[2]);
			}

		};

		using Photonf = Photon<PhotonFloat>;
		using PhotonCollection = Collection<Photonf>;

	public:

		static int photon_size()
		{
			return sizeof(Photonf);
		}

	protected:

		double m_relative_radius;
		double m_radius, m_radius2;
		Vector3f m_pixel_size;
		
		BoundingBox m_bb;
		Vector3i m_size;

		std::vector<PhotonCollection> m_map;

		int index_int(Vector3i const& ijk)const
		{
			return ijk[2] + m_size[2] * (ijk[1] + ijk[0] * m_size[1]);
		}

		Vector3i cell_index(Vector3f const& xyz)const
		{
			return m_size.simdMin(Vector3i(0, 0, 0).simdMax((Vector3i((xyz - m_bb[0]).simdDiv(m_pixel_size)))));
		}
		
		int index(Vector3f const& xyz)const
		{
			return index_int(Vector3i((xyz - m_bb[0]).simdDiv(m_pixel_size)));
		}

		std::mutex mtx;

		void addPhotonToMap(Photonf const& photon)
		{
			int id = index(photon.m_point);
			mtx.lock();
			PhotonCollection& list = m_map[id];
			list.push_back(photon);
			mtx.unlock();
		}

		bool insideMap(Vector3i const& cell_index)const
		{
			return cell_index[0] >= 0 && cell_index[0] < m_size[0] &&
				   cell_index[1] >= 0 && cell_index[1] < m_size[1] &&
				   cell_index[2] >= 0 && cell_index[2] < m_size[2];
		}

		int m_number_of_photons;

	public:

		PhotonMapper(unsigned int sample_per_pixel, unsigned int width, unsigned int height):
			DirectIntegrator(sample_per_pixel, width, height)
		{}

		void buildMap(Scene const& scene, double relative_radius, int pcount)
		{
			m_relative_radius = relative_radius;
			m_bb = scene.m_sceneBoundingBox;
			m_bb[0] -= Vector3f(0.001, 0.001, 0.001);
			m_bb[1] += Vector3f(0.001, 0.001, 0.001);
			Vector3f dim = m_bb.diag();
			double max_dir = dim.simdAbs().max();
			
			m_radius = max_dir * m_relative_radius;
			m_radius2 = m_radius * m_radius;
			m_pixel_size = m_radius * 2.0; 
			Vector3f sizef = dim.simdDiv(m_pixel_size);
			m_size = sizef.ceil();
			int n_cells = m_size.prod();
			assert(n_cells > 0);
			std::cout << "Photon map size: " << m_size << std::endl;
			m_map.resize(n_cells);

			m_number_of_photons = pcount;
			OMP_PARALLEL_FOR
				for (int sample = 0; sample < pcount; ++sample)
				{
					Math::Sampler sampler(sample);
					SurfaceLightSample sls;
					sampleOneLight(scene, sampler, sls, sample);
					
					DirectionSample dirSample = sls.geo->getMaterial()->sampleLightDirection(sls, sampler);
					Hit hit;
					RGBColor beta = dirSample.bsdf / (sls.pdf * dirSample.pdf) * std::abs(dirSample.direction * sls.normal);
					Ray ray(sls.vector, dirSample.direction);
					for (int len = 2; len <= m_max_len-1; ++len)
					{
						if (scene.full_intersection(ray, hit))
						{
							if (!hit.geometry->getMaterial()->delta())
							{
								Photonf photon = { hit, (uint8_t)(len - 2), beta / ((double)m_number_of_photons) };
								addPhotonToMap(photon);
							}

							hit.geometry->getMaterial()->sampleBSDF(hit, 1, 1, dirSample, sampler, true);
							beta *= dirSample.bsdf / dirSample.pdf * std::abs(dirSample.direction * hit.primitive_normal);
							ray = { hit.point, dirSample.direction };
						}
					}
				}
		}

		
		virtual RGBColor sendRay(Scene const& scene, Ray const& pray, Math::Sampler& sampler)const final override
		{
			Ray ray = pray;
			RGBColor beta = scene.m_camera.We<true>(ray.direction()) / scene.m_camera.pdfWeSolidAngle<true>(ray.direction());
			RGBColor res = 0;
			for (int len = 2; len <= m_max_len; ++len)
			{
				Hit hit;
				if (scene.full_intersection(ray, hit))
				{
					const Material& material = *hit.geometry->getMaterial();

					
					res += beta * material.Le(hit.facing, hit.tex_uv);


					if (material.delta())
					{
						DirectionSample next_dir;
						material.sampleBSDF(hit, 1, 1, next_dir, sampler, false);
						ray = { hit.point, next_dir.direction };
						beta *= next_dir.bsdf / next_dir.pdf * std::abs(next_dir.direction * hit.primitive_normal);
					}
					else
					{
						// Add the photons near the hit
						const Vector3i ijk = cell_index(hit.point);
						RGBColor photons_contrib = 0;
						for (int i = -1; i <= 1; ++i)
						{
							for (int j = -1; j <= 1; ++j)
							{
								for (int k = -1; k <= 1; ++k)
								{
									const Vector3i cell = ijk + Vector3i(i, j, k);
									if (insideMap(cell))
									{
										const PhotonCollection& photons = m_map[index_int(cell)];
										for (Photonf const& photon : photons)
										{
											int path_len = photon.m_depth + 2 + len - 1;
											if (path_len <= m_max_len)
											{
												const Vector3f d = hit.point - photon.m_point;
												const double dist2 = d.norm2();
												if (dist2 < m_radius2)
												{
													const Geometry::Material* mat = photon.m_primitive->geometry()->getMaterial();
													if (mat == hit.geometry->getMaterial())
													{
														double kernel = 1 / (Math::pi * m_radius2);
														double attenuation = 1;// -std::sqrt(dist2 / m_radius2);
														photons_contrib += beta * photon.beta() * mat->BSDF(hit, photon.m_dir, -ray.direction()) * attenuation * kernel;
													}
												}
											}
										}
									}
								}
							}
						}
						res += photons_contrib;
						break;
					}


					
				}
				else
				{
					res += beta * scene.getBackgroundColor(ray.direction());
					break;
				}
			}
			return res;
		}

	};
}