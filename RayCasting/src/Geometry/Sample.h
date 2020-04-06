#pragma once

#include <Math/Vectorf.h>
#include <Geometry/RGBColor.h>

namespace Geometry
{
	class GeometryBase;
	class Primitive;
	struct DirectionSample
	{
		double pdf;
		RGBColor bsdf;

		Math::Vector3f direction;
	};

	//describes a point sampled on a light
	struct SurfaceSample
	{
		double pdf;
		const GeometryBase* geo;
		const Primitive* primitive;
		Math::Vector2f uv;
		Math::Vector3f normal;
		Math::Vector3f vector;
	};
}