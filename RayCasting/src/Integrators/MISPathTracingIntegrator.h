#pragma once

#include <Integrators/RayTracingBaseIntegrator.h>
#include <cassert>


namespace Integrator
{
	class MISPathTracingIntegrator: public RayTracingBaseIntegrator
	{
	protected:

	public:

		MISPathTracingIntegrator(unsigned int sample_per_pixel, unsigned int width, unsigned int height) :
			RayTracingBaseIntegrator(sample_per_pixel, width, height)
		{}


		static __forceinline void balance(__out double& w0, __out double& w1, const double p00, const double p01, const double p10, const double p11)
		{
			w0 = p00 / (p00 + p10);
			w1 = p11 / (p01 + p11);			
		}

		static __forceinline void power(__out double& w0, __out double& w1, const double p00, const double p01, const double p10, const double p11, double beta = 2)
		{
			w0 = pow(p00, beta) / (pow(p00, beta)+ pow(p10, beta));
			w1 = pow(p11, beta) / (pow(p01, beta) + pow(p11, beta));
		}


		RGBColor MISAddDirectIllumination(Scene const& scene, Hit const& hit, Math::Sampler& sampler)const
		{
			RGBColor res = 0;

			//sample the surface
			SurfaceSample light_sample;
			//sampleOneLight(scene, sampler, light_sample);
			sampleOneLight(scene, hit, sampler, light_sample);
			Math::Vector3f to_light = (light_sample.vector - hit.point);
			const double dist2 = to_light.norm2();
			const double dist = sqrt(dist2);
			to_light /= dist;
			RGBColor bsdf = hit.geometry->getMaterial()->BSDF(hit, to_light);
			if (!bsdf.isBlack())
			{
				Ray ray(hit.point, to_light);
				Hit light_hit;
				if (scene.full_intersection(ray, light_hit) && samePoint(light_hit, dist))
				{
					RGBColor contribution = bsdf * std::abs(to_light * hit.primitive_normal) * light_hit.geometry->getMaterial()->Le(light_hit.primitive_normal, light_hit.tex_uv, light_hit.to_view);
					double surface_pdf = light_sample.pdf * dist2 / std::abs(light_hit.primitive_normal * to_light);
					res += contribution / (surface_pdf + hit.geometry->getMaterial()->pdf(hit, to_light));
				}
			}
			return res;
		}

		

		RGBColor sendRay(Scene const& scene, Ray const& pray, Math::Sampler& sampler)const final override
		{
			Ray ray(pray);
			RGBColor T = scene.m_camera.We(ray.direction()) / scene.m_camera.pdfWeSolidAngle(ray.direction());
			RGBColor res = 0;
			int len = 1;
			Hit hit, prev_hit;
			bool prev_delta = true;
			double dir_pdf;
			while (len < m_max_len)
			{
				prev_hit = hit;
				if (scene.full_intersection(ray, hit))
				{
					++len;
					if (prev_delta)
					{
						res += T * hit.geometry->getMaterial()->Le(hit.primitive_normal, hit.tex_uv, hit.to_view);
					}
					else
					{
						RGBColor contribution = hit.geometry->getMaterial()->Le(hit.primitive_normal, hit.tex_uv, hit.to_view);

						double surface_pdf = scene.pdfSamplingLight(hit.geometry, prev_hit, hit.point) * hit.z * hit.z / (std::abs(hit.primitive_normal * ray.direction()));

						res += T * contribution / (dir_pdf + surface_pdf) * dir_pdf;
						//res += T * contribution;
					}

					prev_delta = hit.geometry->getMaterial()->delta();

					if (!prev_delta && len < m_max_len)
					{
						res += T * MISAddDirectIllumination(scene, hit, sampler);
					}

					double xi = sampler.generateContinuous<double>();
					if (xi < m_alpha)
					{
						DirectionSample next_dir;
						hit.geometry->getMaterial()->sampleBSDF(hit, next_dir, sampler);
						T *= next_dir.bsdf * std::abs(hit.primitive_normal * next_dir.direction) / next_dir.pdf / m_alpha;
						ray = Ray(hit.point, next_dir.direction);
						dir_pdf = next_dir.pdf;
					}
					else
					{
						break;
					}
					if (T.isBlack())
						break;
				}
				else
				{
					res += T * scene.getBackgroundColor(ray.direction());
					break;
				}
			}
			return res;
		}


		
	};
}