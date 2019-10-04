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



		//void setPass(unsigned int pass)
		//{
		//	m_maximum_pass = pass;
		//}

		//void setPass(unsigned int sub_pixel_division, unsigned int pass_per_pixel)
		//{
		//	m_maximum_pass = sub_pixel_division * sub_pixel_division * pass_per_pixel;
		//}


	};
}