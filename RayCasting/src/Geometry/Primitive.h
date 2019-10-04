#pragma once


#include <Math/Vectorf.h>
#include <Geometry/BoundingBox.h>

namespace Geometry
{
	//Just a guideline for BVH
	class Primitive
	{
	protected:

	public:

		virtual Math::Vector3f center()const = 0;

		virtual BoundingBox box()const = 0;
	};
}