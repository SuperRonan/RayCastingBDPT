#pragma once


#include <Math/Vectorf.h>
#include <Geometry/BoundingBox.h>

namespace Geometry
{
	class Primitive
	{
	protected:

	public:

		virtual Math::Vector3f center()const = 0;

		virtual BoundingBox box()const = 0;

		virtual const GeometryBase* geometry()const = 0;

		virtual Math::Vector2f uv(Math::Vector3f const& point)const = 0;
		virtual Math::Vector2f tuv(Math::Vector2f const& uv)const = 0;

		virtual Math::Vector3f point(Math::Vector3f const& uv)const = 0;
		
		virtual Math::Vector3f normal(Math::Vector3f const& point, Math::Vector2f const& uv)const = 0;
		virtual Math::Vector3f shading_normal(Math::Vector3f const& point, Math::Vector2f const& uv)const = 0;
	};
}