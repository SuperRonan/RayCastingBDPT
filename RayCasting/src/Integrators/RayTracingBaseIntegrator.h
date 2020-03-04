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
						double cosl = std::abs(light_hit.primitive_normal * (-ray.direction()));

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

		RGBColor addOneDirectIllumination(Scene const& scene, Hit const& hit, Math::Sampler& sampler)const
		{
			Geometry::SurfaceLightSample sample;
			sampleOneLight(scene, hit, sampler, sample); // Clever
			//sampleOneLight(scene, sampler, sample); // Uniform
			Math::Vector3f to_light = sample.vector - hit.point;
			const double dist2 = to_light.norm2();
			const double dist = std::sqrt(dist2);
			to_light /= dist;
			const RGBColor bsdf = hit.geometry->getMaterial()->BSDF(hit, to_light, hit.to_view);
			const double cos_on_light = -(sample.normal * to_light);
			const double cos_theta = std::abs(to_light * hit.normal);
			const RGBColor Le = sample.geo->getMaterial()->Le(cos_on_light > 0, sample.uv);
			const RGBColor contrib = bsdf * Le * std::abs(cos_on_light) * cos_theta / dist2;
			if (contrib.isBlack() || contrib.anythingWrong())
				return 0;
			Hit light_hit;
			Ray ray(hit.point, to_light);
			if (scene.full_intersection(ray, light_hit) && samePoint(light_hit, dist, sample.geo))
			{
				return contrib / sample.pdf;
			}
			return 0;
		}





	};
}