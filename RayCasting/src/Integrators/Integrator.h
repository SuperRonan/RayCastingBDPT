#pragma once


#include <Geometry/Scene.h>
#include <Geometry/RGBColor.h>
#include <Visualizer/Visualizer.h>
#include <string>
#include <Image/Image.h>
#include <settings.h>
#include <Auto/RenderResult.h>
#include <System/ProgressReporter.h>

namespace Integrator
{
	using namespace Geometry;

	class Integrator
	{
		
	public:
		unsigned int m_sample_per_pixel;

	protected:
		unsigned int m_max_len;

		
		

	public:
		double m_alpha = 1;

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

		virtual void setLen(unsigned int len)
		{
			m_max_len = len;
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

	};
}