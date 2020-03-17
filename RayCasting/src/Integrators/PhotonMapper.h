#pragma once

#include <Integrators/Integrator.h>
#include <Integrators/DirectIntegrator.h>
#include <Math/Vectorf.h>
#include <mutex>
#include <tbb/mutex.h>
#include <Integrators/PhotonMap.h>

namespace Integrator
{
	class PhotonMapper : public DirectIntegrator
	{
	protected:

		using Vector3i = Math::Vector<int, 3>;
		using Vector3f = Math::Vector3f;
		using Vector2f = Math::Vector2f;

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

			Math::Vector<Float, 3> point()const
			{
				return m_point;
			}

		};

		using Photonf = Photon<PhotonFloat>;

	public:

		static int photon_size()
		{
			return sizeof(Photonf);
		}

	protected:

		PhotonMap<Photonf> m_map;

		double m_relative_radius;
		double m_radius, m_radius2;

		int m_number_of_photons;

		double m_alpha, m_beta;

		double kernel(double dist2)const
		{
			return 1.0 / (Math::pi * m_radius2);
			double dist = std::sqrt(dist2);
			return m_alpha - dist * m_beta;
		}

	public:

		PhotonMapper(unsigned int sample_per_pixel, unsigned int width, unsigned int height):
			DirectIntegrator(sample_per_pixel, width, height)
		{}

		void buildMap(Scene const& scene, double relative_radius, int pcount)
		{
			tic();
			m_relative_radius = relative_radius;
			BoundingBox m_bb = scene.m_sceneBoundingBox;
			m_bb[0] -= Vector3f(0.001, 0.001, 0.001);
			m_bb[1] += Vector3f(0.001, 0.001, 0.001);
			Vector3f dim = m_bb.diag();
			double max_dir = dim.simdAbs().max();
			
			m_radius = max_dir * m_relative_radius;
			m_radius2 = m_radius * m_radius;
			Math::Vector3f m_pixel_size = m_radius * 2.0; 
			Vector3f sizef = dim.simdDiv(m_pixel_size);
			Vector3i m_size = sizef.ceil();
			
			m_map.init(m_bb, m_size);

			m_alpha = 3.0 / (Math::pi * m_radius2);
			m_beta = 3.0 / (Math::pi * m_radius2 * m_radius);

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
								m_map.addPhoton(photon);
							}

							hit.geometry->getMaterial()->sampleBSDF(hit, 1, 1, dirSample, sampler, true);
							beta *= dirSample.bsdf / dirSample.pdf * std::abs(dirSample.direction * hit.primitive_normal);
							ray = { hit.point, dirSample.direction };
						}
					}
				}
			toc();
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
						RGBColor photons_contrib;
						m_map.loopThroughPhotons<std::function<void(Photonf const&)>>([&](Photonf const& photon) {
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
										const double k = kernel(dist2);
										photons_contrib += beta * photon.beta() * mat->BSDF(hit, photon.m_dir, -ray.direction()) * k;
									}
								}
							}
							}, hit.point);
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