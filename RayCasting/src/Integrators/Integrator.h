#pragma once


#include <Geometry/Scene.h>
#include <Geometry/RGBColor.h>
#include <Visualizer/Visualizer.h>
#include <string>
#include <Image/Image.h>
#include <settings.h>
#include <Auto/RenderResult.h>

namespace Integrator
{
	using namespace Geometry;

	class Integrator
	{
		
	public:
		unsigned int m_sample_per_pixel;

	protected:
		unsigned int m_max_depth;

		
		

	public:
		double m_alpha = 0.5;

	protected:

		static __forceinline size_t pixelSeed(size_t x, size_t y, size_t width, size_t height, size_t pass)
		{
#ifdef SAMPLER_BIAS
			return (pass) * (width * height);
#else
			return (y * width + x) + (pass) * (width * height);
#endif
		}

		
		//__forceinline double subPixel(size_t sub_p)const
		//{
		//	return double(2*sub_p + 1) / double(2*m_sub_pixel_samples);
		//}

	public:

		virtual void setDepth(unsigned int d)
		{
			m_max_depth = d;
		}

		virtual void render(Scene const&, Visualizer::Visualizer &) = 0;

		virtual void fastRender(Scene const&, Visualizer::Visualizer &) = 0;

		virtual void debug(Scene const&, Visualizer::Visualizer &) = 0;

		virtual void render(Scene const&, size_t width, size_t height, Auto::RenderResult & res) = 0;




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

		//clever version
		bool sampleOneLight(Scene const& scene, Hit const& hit, Math::Sampler& sampler, SurfaceLightSample& res)const
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
		bool sampleOneLight(Scene const& scene, Math::Sampler& sampler, SurfaceLightSample& res, int index=0)const
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
			light->sampleLight(res, sampler, index);
			res.pdf *= pdf;

			return true;
		}

	};
}