#pragma once

#include <Integrators/DirectIntegrator.h>
#include <Geometry/GeometryBase.h>


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


		RGBColor addOneDirectIllumination(Scene const& scene, Hit const& hit, Math::Sampler& sampler)const
		{
			Geometry::SurfaceSample sample;
			sampleOneLight(scene, hit, sampler, sample); // Clever
			//sampleOneLight(scene, sampler, sample); // Uniform
			Math::Vector3f to_light = sample.vector - hit.point;
			const double dist2 = to_light.norm2();
			const double dist = std::sqrt(dist2);
			to_light /= dist;
			const RGBColor bsdf = hit.geometry->getMaterial()->BSDF(hit, to_light, hit.to_view);
			const double cos_on_light = -(sample.normal * to_light);
			const double cos_theta = std::abs(to_light * hit.normal);
			const RGBColor Le = sample.geo->getMaterial()->Le(sample.normal, sample.uv, -to_light);
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