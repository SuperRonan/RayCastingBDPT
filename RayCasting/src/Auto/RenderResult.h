#pragma once

#include <Image/Image.h>
#include <Geometry/RGBColor.h>

namespace Auto
{
	class RenderResult
	{
	public:
		Image::Image<Geometry::RGBColor> image;

		double time;

		RenderResult() :
			image(),
			time(0)
		{
			
		}
	};
}