#pragma once
#include <Math/Vector.h>
#include <Math/Vectorf.h>

namespace Geometry
{


	class GeometryBase;
	class Primitive;

	class Hit
	{
	

	public:
		double z;

		const GeometryBase * geometry;
		const Primitive * primitve;

		Math::Vector3f point;
		
		Math::Vector3f to_view;
		
		bool facing;

		Math::Vector3f normal;
		Math::Vector3f primitive_normal;
		
		Math::Vector2f primitive_uv;
		Math::Vector2f tex_uv;

		//Math::Vector3f reflected;

		Hit():
			z(-1),
			geometry(nullptr),
			primitve(nullptr)
		{}

		Hit(Hit const& other) = default;


		Hit & operator=(Hit const& other)noexcept
		{
			std::memcpy(this, &other, sizeof(Hit));
			return *this;
			z = other.z;
			geometry = other.geometry;
			primitve = other.primitve;
			point = other.point;
			to_view = other.to_view;
			facing = other.facing;
			normal = other.normal;
			primitive_normal = other.primitive_normal;
			primitive_uv = other.primitive_uv;
			tex_uv = other.tex_uv;
			//reflected = other.reflected;
			return *this;
		}

		void set_normal(Math::Vector3f const& vec)
		{
			normal = vec;
			//reflected = normal * (2 * (normal* to_view)) - to_view;
		}

		Math::Vector3f primitive_reflected()const
		{
			return primitive_normal.reflect(to_view);
		}

		Math::Vector3f shading_refelcted()const
		{
			return normal.reflect(to_view);
		}

		bool valid()const
		{
			return z != -1;
		}
	};
}