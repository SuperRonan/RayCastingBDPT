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
			const double p_bsdf = 0.5;
			const double p_surface = 1.0 - p_bsdf;
			double xi = sampler.generateContinuous<double>();
			if (xi < p_bsdf)
			{
				//sample the bsdf
				DirectionSample next_dir;
				hit.geometry->getMaterial()->sampleBSDF(hit, 1, 1, next_dir, sampler);
				RGBColor contribution = next_dir.bsdf * std::abs(next_dir.direction * hit.primitive_normal);
				if (!contribution.isBlack())
				{
					Ray ray(hit.point, next_dir.direction);
					Hit light_hit;
					if (scene.full_intersection(ray, light_hit))
					{
						contribution *= light_hit.geometry->getMaterial()->Le(light_hit.facing, light_hit.tex_uv);
						
						double surface_pdf = scene.pdfSamplingLight(light_hit.geometry) * light_hit.z * light_hit.z / (std::abs(light_hit.primitive_normal * ray.direction()));

						res += contribution / (next_dir.pdf +surface_pdf) / p_bsdf;
					}
				}
			}
			else
			{
				//sample the surface
				SurfaceLightSample light_sample;
				sampleOneLight(scene, sampler, light_sample);
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
						RGBColor contribution = bsdf * std::abs(to_light * hit.primitive_normal) * light_hit.geometry->getMaterial()->Le(light_hit.facing, light_hit.tex_uv);
						double surface_pdf = light_sample.pdf * dist2 / std::abs(light_hit.primitive_normal * to_light);
						res += contribution / (surface_pdf + hit.geometry->getMaterial()->pdf(hit, to_light)) / p_surface;
					}
				}
			}
			return res;
		}

		

		RGBColor sendRay(Scene const& scene, Ray const& pray, Math::Sampler& sampler)const final override
		{
			Ray ray(pray);
			RGBColor T = scene.m_camera.We(ray.direction()) / scene.m_camera.pdfWeSolidAngle(ray.direction());
			RGBColor res = 0;
			int depth = 0;
			Hit hit;
			bool use_direct = true;
			while (depth <= m_max_depth)
			{
				if (scene.full_intersection(ray, hit))
				{
					if (use_direct)
					{
						res += T * hit.geometry->getMaterial()->Le(hit.facing, hit.tex_uv);
					}

					use_direct = hit.geometry->getMaterial()->delta();

					if (!use_direct)
					{
						res += T * MISAddDirectIllumination(scene, hit, sampler);
					}

					double xi = sampler.generateContinuous<double>();
					if (xi < m_alpha)
					{
						DirectionSample next_dir;
						hit.geometry->getMaterial()->sampleBSDF(hit, 1, 1, next_dir, sampler);
						T *= next_dir.bsdf * std::abs(hit.primitive_normal * next_dir.direction) / next_dir.pdf / m_alpha;
						ray = Ray(hit.point, next_dir.direction);
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
					res += scene.getBackgroundColor(ray.direction());
					break;
				}
				++depth;
			}
			return res;
		}


		
	};
}