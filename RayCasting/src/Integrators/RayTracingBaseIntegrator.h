#pragma once

#include <Integrators/DirectIntegrator.h>
#include <Geometry/GeometryBase.h>
#include <Geometry/BoundedStack.h>


namespace Integrator
{
	class RayTracingBaseIntegrator : public DirectIntegrator
	{
	public:

		unsigned int m_diffuse_samples = 1;
		unsigned int m_specular_samples = 1;




	protected:

		__forceinline static bool samePoint(Hit const& hit, double dist, const GeometryBase* geo)
		{
#ifdef TRICK_DIRECT
			return hit.geometry == geo;
#else
			return samePoint(hit, dist);
#endif
		}

		__forceinline static bool samePoint(Hit const& hit, double dist)
		{
			return abs(dist - hit.z) < 0.00000001;
		}

	public:

		RayTracingBaseIntegrator(unsigned int sample_per_pixel, unsigned int width, unsigned int height) :
			DirectIntegrator(sample_per_pixel, width, height),
			m_diffuse_samples(1),
			m_specular_samples(1)
		{
		}


		RGBColor addDirectIllumination(Scene const& scene, Hit const& hit, Math::Sampler& sampler)const
		{
			LightSampleStack ps;

			RGBColor res;
			for (const GeometryBase* const sl : scene.m_surface_lights)
			{
				if (sl != hit.geometry)
				{
					sl->sampleLights(ps, sampler, m_direct_samples); //old
					//sl->sampleLights(ps, hit, sampler, m_direct_samples); //new clever but biased (need to correct it)
				}
			}

			//switching from point to direction from the hit
			StackN<double> distStack;
			for (SurfaceLightSample& sls : ps)
			{
				Math::Vector3f dir = (sls.vector - hit.point);
				double dist = dir.norm();
				distStack.push(dist);
				sls.vector = dir / dist;
			}

			ColorStack bsdfs;
			hit.geometry->getMaterial()->BSDF(hit, bsdfs, ps);

			for (size_t i = 0; i < bsdfs.size(); ++i)
			{
				RGBColor const& bsdf = bsdfs[i];
				const double cosi = std::abs(hit.primitive_normal * ps[i].vector);
				RGBColor prod = bsdf * cosi;
				if (!prod.isBlack()) //or almost black ???
				{
					Hit light_hit;
					Ray ray(hit.point, ps[i].vector);
					if (scene.full_intersection(ray, light_hit) && samePoint(light_hit, distStack[i], ps[i].geo))
					{
						RGBColor color = light_hit.geometry->getMaterial()->Le(light_hit.facing, light_hit.tex_uv);
						double dist2 = light_hit.z * light_hit.z;
						double cosl = light_hit.primitive_normal * (-ray.direction());

						res += prod * cosl * color / (ps[i].pdf * m_direct_samples * dist2);
					}
				}
			}
			if (res.anythingWrong())
			{
				return 0;
			}
			return res;
		}

		//clever version
		bool sampleOneLight(Scene const& scene, Hit const& hit, Math::Sampler & sampler, SurfaceLightSample& res)const
		{
			//select one light
			const GeometryBase* light;
			double pdf;

			
			if (!scene.sampleOneLight(sampler, pdf, light))
			{
				return false;
			}

			/*
			if (light == hit.geometry)
			{
				return false;
			}
			*/

			//select a point on this light
			//light->sampleLight(res, hit, sampler);
			light->sampleLight(res, sampler);
			res.pdf *= pdf;

			return true;
		}


		//old version
		bool sampleOneLight(Scene const& scene, Math::Sampler& sampler, SurfaceLightSample& res)const
		{
			if (scene.m_surface_lights.empty())
			{
				return false;
			}

			//select one light
			const GeometryBase* light;
			double pdf;
			scene.sampleOneLight(sampler, pdf, light);

			//select a point on this light
			light->sampleLight(res, sampler);
			res.pdf *= pdf;

			return true;
		}



	};
}