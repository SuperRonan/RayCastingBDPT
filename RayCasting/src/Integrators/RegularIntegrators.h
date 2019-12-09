#pragma once

#include <Integrators/DirectIntegrator.h>


namespace Integrator
{
	class NormalIntegrator final: public DirectIntegrator
	{
	protected:

	public:

		NormalIntegrator(unsigned int sample_per_pixel, unsigned int width, unsigned int height):
			DirectIntegrator(sample_per_pixel, width, height)
		{}

		RGBColor sendRay(Scene const& scene, Ray const& ray, Math::Sampler& sampler)const final override
		{
			Hit hit;

			if (scene.full_intersection(ray, hit))
			{
				RGBColor res;
				res[0] = 0.5 + hit.normal[0] * 0.5;
				res[1] = 0.5 + hit.normal[1] * 0.5;
				res[2] = 0.5 + hit.normal[2] * 0.5;
				return res;
			}
			return 0;
			RGBColor res;
			res[0] = 0.5 + -ray.direction()[0] * 0.5;
			res[1] = 0.5 + -ray.direction()[1] * 0.5;
			res[2] = 0.5 + -ray.direction()[2] * 0.5;
			return res;
			
		}
	};

	class UVIntegrator final: public DirectIntegrator
	{
	protected:

	public:

		UVIntegrator(unsigned int sample_per_pixel, unsigned int width, unsigned int height) :
			DirectIntegrator(sample_per_pixel, width, height)
		{}

		RGBColor sendRay(Scene const& scene, Ray const& ray, Math::Sampler& sampler)const final override
		{
			Hit hit;

			if (scene.full_intersection(ray, hit))
			{
				RGBColor res;
				Math::Vector2f uv = hit.tex_uv;
				uv = uv.simdMul(uv);
				res[0] = uv[0];
				res[1] = 0.02;
				res[2] = uv[1];
				return res;
			}
			return 0;
		}
	};

	class AmbientOcclusionIntegrator final : public DirectIntegrator
	{
	protected:

	public:

		AmbientOcclusionIntegrator(unsigned int sample_per_pixel, unsigned int width, unsigned int height) :
			DirectIntegrator(sample_per_pixel, width, height)
		{}

		RGBColor sendRay(Scene const& scene, Ray const& ray, Math::Sampler& sampler)const final override
		{
			Hit hit;
			RGBColor sky = 1;
			if (scene.full_intersection(ray, hit))
			{
				Math::RandomDirection dir_sampler(&sampler, hit.primitive_normal);
				Math::Vector3f dir = dir_sampler.generate();
				double cos_dir = dir * hit.primitive_normal;
				if (cos_dir < 0)
				{
					cos_dir = -cos_dir;
					dir = -dir;
				}
				double pdf = cos_dir / Math::pi;

				Ray ray(hit.point, dir);
				if (scene.full_intersection(ray, hit))
				{
					return 0;
				}
				else
				{
					return 1.0;
				}
			}
			return 0;
		}
	};

	class MIDIntegrator final: public DirectIntegrator
	{
	protected:

	public:

		MIDIntegrator(unsigned int sample_per_pixel, unsigned int width, unsigned int height) :
			DirectIntegrator(sample_per_pixel, width, height)
		{}

		RGBColor sendRay(Scene const& scene, Ray const& ray, Math::Sampler& sampler)const final override
		{
			Hit hit;

			if (scene.full_intersection(ray, hit))
			{
				return hit.geometry->getMaterial()->ID_COLOR() * std::abs(hit.normal * ray.direction());
			}
			return 0;
		}
	};

	class BOXIntegrator final: public DirectIntegrator 
	{
	protected:

	public:

		BOXIntegrator(unsigned int sample_per_pixel, unsigned int width, unsigned int height) :
			DirectIntegrator(sample_per_pixel, width, height)
		{}

		RGBColor sendRay(Scene const& scene, Ray const& ray, Math::Sampler& sampler)const final override
		{
			return scene.count_box(ray);
		}
	};
}