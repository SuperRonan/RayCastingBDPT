#pragma once

#include <Integrators/Integrator.h>

namespace Integrator
{
	class PhotonMapper : public Integrator
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
			Math::Vector<Float, 2> m_dir;
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
		double m_pixel_size;
		
		BoundingBox m_bb;
		Vector3i m_size;

		std::vector<Photonf> m_map;

		int index(Vector3i const& ijk)
		{
			return ijk[2] + m_size[2] * (ijk[1] + ijk[0] * m_size[1]);
		}
		
		int index(Vector3f const& xyz)
		{
			return index(Vector3i(xyz / m_pixel_size));
		}

		void addPhotonToMap(Photonf const& photon)
		{

		}

		int m_number_of_photons;

	public:

		PhotonMapper(double relative_radius):
			Integrator(),
			m_relative_radius(relative_radius)
		{}

		int buildMap(Scene const& scene, int pcount, Math::Sampler & sampler)
		{
			m_bb = scene.m_sceneBoundingBox;
			Vector3f dim = m_bb.diag();
			double max_dir = dim.simdAbs().max();
			
			m_radius = max_dir * m_relative_radius;
			m_radius2 = m_radius * m_radius;
			m_pixel_size = m_radius * 2.0; 
			Vector3f sizef = dim / m_pixel_size;
			m_size = sizef.ceil();
			int n_cells = m_size.prod();
			m_map.resize(n_cells);

			m_number_of_photons = pcount;
			OMP_PARALLEL_FOR
				for (int sample = 0; sample < pcount; ++sample)
				{
					SurfaceLightSample sls;
					sampleOneLight(scene, sampler, sls, sample);
					
					DirectionSample dirSample = sls.geo->getMaterial()->sampleLightDirection(sls, sampler);
					Hit hit;
					RGBColor beta = dirSample.bsdf / (sls.pdf * dirSample.pdf);
					Ray ray(sls.vector, dirSample.direction);
					for (int len = 2; len <= m_max_len; ++len)
					{
						if (scene.full_intersection(ray, hit))
						{
							if (!hit.geometry->getMaterial()->delta())
							{
								Photonf photon = { hit, len - 2, beta };

							}
						}
					}
				}
		}


	};
}