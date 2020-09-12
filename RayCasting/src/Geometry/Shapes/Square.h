#ifndef _Geometry_Square_H
#define _Geometry_Square_H

#include <Geometry/Shapes/Geometry.h>
#include <algorithm>

namespace Geometry
{
	////////////////////////////////////////////////////////////////////////////////////////////////////
	/// \class	Square
	///
	/// \brief	A square on plane (X,Y), centered in (0,0,0).
	///
	/// \author	F. Lamarche, Université de Rennes 1
	/// \date	04/12/2013
	////////////////////////////////////////////////////////////////////////////////////////////////////
	class Square : public GeometryCollection
	{

	protected:

		unsigned int m_u_div=1, m_v_div=1;


	public:

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	Square::Square(Material * material)
		///
		/// \brief	Constructor.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	04/12/2013
		///
		/// \param [in,out]	material	If non-null, the material.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		Square(Material * material, Math::Vector3f const& origin = { 0, 0, 0 }, Math::Vector3f const& u = { 1, 0, 0 }, Math::Vector3f const& v = { 0, 1, 0 })
			: GeometryCollection(material)
		{
			int p0 = addVertex(origin + u * 0.5 + v * 0.5) ;
			int p1 = addVertex(origin + u * 0.5 - v * 0.5) ;
			int p2 = addVertex(origin - u * 0.5 + v * 0.5) ;
			int p3 = addVertex(origin - u * 0.5 - v * 0.5) ;
			addTriangle(p0,p1,p2) ;
			addTriangle(p3,p2,p1) ;
			m_surface = u.norm() * v.norm();
		}

		virtual bool build_lights()
		{
			Math::Vector3f u_axis = (m_vertices[1] - m_vertices[0]), v_axis = (m_vertices[2] - m_vertices[0]);
			m_current_sum = sqrt(u_axis.norm2() * v_axis.norm2());
			m_surface = m_current_sum;
			return m_material->is_emissive();
		}

		virtual void divide(unsigned int div=1)
		{
			unsigned int sqrt_div = std::max(1.0, std::sqrt(div));
			if (sqrt_div * sqrt_div == div)
			{
				m_u_div = sqrt_div;
				m_v_div = sqrt_div;
			}
			else
			{
				m_u_div = div;
				m_v_div = 1;
			}

			m_divisions = m_u_div * m_v_div;

			check_capacity(div);
		}

		virtual void sampleLight(SurfaceSample& res, Math::Sampler & sampler, unsigned int i=0)const override
		{
			
			unsigned int offset = (m_offset + i) % m_divisions;
			double sub_u = sampler.generateContinuous<double>();
			double sub_v = sampler.generateContinuous<double>();
			

			unsigned int index_u = offset / m_v_div;
			unsigned int index_v = offset % m_v_div;

			double u = (double(index_u) + sub_u) / double(m_u_div);
			double v = (double(index_v) + sub_v) / double(m_v_div);


			const Triangle* tri = u + v > 1 ? &m_triangles[1] : &m_triangles[0];
			Math::Vector3f u_axis = (m_vertices[1] - m_vertices[0]), v_axis = (m_vertices[2] - m_vertices[0]);
			Math::Vector3f point = m_vertices[0] + u_axis * u + v_axis * v;
			res = { 1.0 / surface(), this, tri, {u, v}, GeometryCollection::m_triangles[0].normal(), point };
			
		}

		virtual void sampleLight(SurfaceSample& res, double u, double v)const override
		{
			const Triangle* tri = u + v > 1 ? &m_triangles[1] : &m_triangles[0];
			Math::Vector3f u_axis = (m_vertices[1] - m_vertices[0]), v_axis = (m_vertices[2] - m_vertices[0]);
			Math::Vector3f point = m_vertices[0] + u_axis * u + v_axis * v;
			res = { 1.0 / surface(), this, tri, {u, v}, GeometryCollection::m_triangles[0].normal(), point };

		}

	} ;





	class OneTriangle : public GeometryCollection
	{
	public:

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	Square::Square(Material * material)
		///
		/// \brief	Constructor.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	04/12/2013
		///
		/// \param [in,out]	material	If non-null, the material.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		OneTriangle(Material * material)
			: GeometryCollection(material)
		{
			int p0 = addVertex(Math::makeVector(0.5, 0.5, 0.0));
			int p1 = addVertex(Math::makeVector(0.5, -0.5, 0.0));
			int p2 = addVertex(Math::makeVector(-0.5, 0.5, 0.0));
			addTriangle(p0, p1, p2);
		}
	};
}

#endif
