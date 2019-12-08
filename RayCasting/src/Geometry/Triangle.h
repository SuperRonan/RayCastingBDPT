#ifndef _Geometry_Triangle
#define _Geometry_Triangle

#include <Math/Vectorf.h>
#include <Geometry/Ray.h>
#include <cassert>
//#include <GeometryCollection/generic_texture.hpp>
//#include <GeometryCollection/TriangleTetxture.h>
#include <Geometry/TriangleMultipleTexture.hpp>
#include <valarray>
#include <Math/Sampler.h>
#include <Geometry/Primitive.h>




namespace Geometry
{

	class GeometryBase;

	class Triangle : public Primitive
	{
	protected:
		/// \brief	Pointers to the three vertices
		Math::Vector3f * m_vertex[3] ; 
		/// \brief Pointers to the texture coordinates
		Math::Vector2f* m_textureCoordinates[3] = { nullptr, nullptr, nullptr };
		/// \brief	The vertex 0 (added to enhance cache consistency)
		Math::Vector3f m_vertex0 ;
		/// \brief	The u axis.
		Math::Vector3f m_uAxis ; 
		/// \brief	The v axis.
		Math::Vector3f m_vAxis ;
		/// \brief	The normal.
		Math::Vector3f m_normal ;
		
		/// \brief	Per vertex normal.
		Math::Vector3f m_vertexNormal[3];

		const GeometryBase * m_geometry;

		

	public:

		

		Math::Vector3f center()const final override
		{
			return (vertex(0) + vertex(1) + vertex(2)) / 3;
		}

		BoundingBox box()const final override
		{
			BoundingBox res(m_vertex0, m_vertex0);
			res.update(vertex(1));
			res.update(vertex(2));
			return res;
		}



		/// <summary>
		/// Sets a vertex normal.
		/// </summary>
		/// <param name="index">The index of the vertex.</param>
		/// <param name="normal">The normal.</param>
		void setVertexNormal(unsigned int index, Math::Vector3f const & normal)
		{
			m_vertexNormal[index] = normal;
		}


		/// <summary>
		/// Gets a vertex normal.
		/// </summary>
		/// <param name="index">The index of the vertex.</param>
		/// <returns></returns>
		const Math::Vector3f & getVertexNormal(unsigned int index) const
		{
			return m_vertexNormal[index];
		}

		/// <summary>
		/// Gets the vertex normal oriented toward a given point.
		/// </summary>
		/// <param name="index">The index of the vertex.</param>
		/// <param name="toward">The point.</param>
		/// <returns></returns>
		Math::Vector3f getVertexNormal(unsigned int index, const Math::Vector3f & toward) const
		{
			const Math::Vector3f & normal = m_vertexNormal[index];
			if ((toward - vertex(index))*normal < 0.0)
			{
				return -normal;
			}
			return normal;
		}


		Math::Vector3f getVertexNormal(unsigned int index, bool facing) const
		{
			const Math::Vector3f & normal = m_vertexNormal[index];
			if (!facing)
			{
				return -normal;
			}
			return normal;
		}

		/// <summary>
		/// Gets a pointer to the vertex normals array (size = 3, one normal per vertex).
		/// </summary>
		/// <returns></returns>
		const Math::Vector3f * getVertexNormals() const
		{
			return m_vertexNormal;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	void Triangle::update()
		///
		/// \brief	Updates precomputed data. This method should be called if vertices are externally 
		/// 		modified.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	04/12/2013
		////////////////////////////////////////////////////////////////////////////////////////////////////
		void update()
		{
			m_vertex0 = *m_vertex[0] ;
			m_uAxis = (*m_vertex[1])-(*m_vertex[0]) ;
			m_vAxis = (*m_vertex[2])-(*m_vertex[0]) ;
			m_normal = m_uAxis^m_vAxis ;
			m_normal = m_normal*(1.0f/m_normal.norm()) ;
			m_vertexNormal[0] = m_normal;
			m_vertexNormal[1] = m_normal;
			m_vertexNormal[2] = m_normal;
		}

		/// <summary>
		/// Initializes a new instance of the <see cref="Triangle"/> class.
		/// </summary>
		/// <param name="a">A pointer to the first vertex.</param>
		/// <param name="b">A pointer to the second vertex.</param>
		/// <param name="c">A pointer to the third vertex.</param>
		/// <param name="ta">The texture coordinates of the first vertex.</param>
		/// <param name="tb">The texture coordinates of the second vertex.</param>
		/// <param name="tc">The texture coordinates of the third vertex.</param>
		/// <param name="material">The material.</param>
		Triangle(Math::Vector3f * a, Math::Vector3f * b, Math::Vector3f * c,
			Math::Vector2f * ta, Math::Vector2f * tb, Math::Vector2f * tc,
			const GeometryBase * geo,
			const Math::Vector3f * normals = nullptr) :
			m_geometry(geo)
		{
			m_vertex[0] = a;
			m_vertex[1] = b;
			m_vertex[2] = c;
			m_textureCoordinates[0] = ta;
			m_textureCoordinates[1] = tb;
			m_textureCoordinates[2] = tc;
			//m_material = material;
			update();
			if (normals != nullptr)
			{
				::std::copy(normals, normals+3, m_vertexNormal);
			}
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	Triangle::Triangle(Math::Vector3 * a, Math::Vector3 * b, Math::Vector3 * c,
		/// 	Material * material)
		///
		/// \brief	Constructor.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	04/12/2013
		///
		/// \param	a					A pointer to the first vertex.
		/// \param	b					A pointer to the second vertex.
		/// \param	c					A pointer to the third vertex.
		/// \param [in,out]	material	If non-null, the material.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		Triangle(Math::Vector3f * a, Math::Vector3f * b, Math::Vector3f * c, const GeometryBase * geo, const Math::Vector3f * normals = nullptr)
			:
			m_geometry(geo)
		{
			m_vertex[0] = a ;
			m_vertex[1] = b ;
			m_vertex[2] = c ;
			std::fill(m_textureCoordinates, m_textureCoordinates + 3, nullptr);
			update() ;
			if (normals != nullptr)
			{
				::std::copy(normals, normals + 3, m_vertexNormal);
			}
		}


		~Triangle()
		{
			
		}


		Triangle(Triangle const& other) :
			m_geometry(other.m_geometry)
		{
			//std::cout << "copy constructor!" << std::endl;
			//TODO
			m_vertex[0] = other.m_vertex[0];
			m_vertex[1] = other.m_vertex[1];
			m_vertex[2] = other.m_vertex[2];
			std::copy(other.m_textureCoordinates, other.m_textureCoordinates + 3, m_textureCoordinates);
			//std::fill(m_textureCoordinates, m_textureCoordinates + 3, nullptr);
			update();
			::std::copy(other.m_vertexNormal, other.m_vertexNormal + 3, m_vertexNormal);

			
		}


		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	Triangle::Triangle()
		///
		/// \brief	Default constructor.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	04/12/2013
		////////////////////////////////////////////////////////////////////////////////////////////////////
		Triangle()
			:
			m_geometry(nullptr)
		{
			std::fill(m_vertex, m_vertex + 3, nullptr);
			std::fill(m_textureCoordinates, m_textureCoordinates + 3, nullptr);
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	Material * Triangle::material() const
		///
		/// \brief	Gets the material.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	04/12/2013
		///
		/// \return	null if it fails, else.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		


		const GeometryBase * geometry()const
		{
			return m_geometry;
		}

		const GeometryBase * & geometry()
		{
			return m_geometry;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	Math::Vector3 const & Triangle::vertex(int i) const
		///
		/// \brief	Gets the ith vertex
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	04/12/2013
		///
		/// \param	i	Vertex index in [0;2].
		///
		/// \return	the vertex.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		Math::Vector3f const & vertex(int i) const
		{ 
			assert(i >= 0 && i < 3);
			return *(m_vertex[i]) ; 
		}

		Math::Vector3f const& vertex0()const
		{
			return m_vertex0;
		}


		/// <summary>
		/// Gets the textures coordinates of a vertex.
		/// </summary>
		/// <param name="i">The i.</param>
		/// <returns></returns>
		Math::Vector2f const & textureCoordinate(int i) const
		{
			assert(i >= 0 && i < 3);
			//assert(m_textureCoordinates[i] != nullptr);
			return *(m_textureCoordinates[i]);
		}

		/// <summary>
		/// Interpolates the texture coordinate given the u,v coordinates of an intersection.
		/// </summary>
		/// <param name="u">The u.</param>
		/// <param name="v">The v.</param>
		/// <returns></returns>
		Math::Vector2f interpolateTextureCoordinate(double u, double v) const
		{
			if (has_texture_coordinates())
				return textureCoordinate(0)* (1 - u - v) + textureCoordinate(1) * u + textureCoordinate(2) * v;
			return { u, v };
		}


		bool has_texture_coordinates()const
		{
			return m_textureCoordinates[0] != nullptr;
		}

		/// <summary>
		/// Samples the texture given the u,v coordinates of an intersection.
		/// </summary>
		/// <param name="u">The u.</param>
		/// <param name="v">The v.</param>
		/// <returns>The color of the texture at the given u,v coordinates or (1.0,1.0,1.0) if no texture is associated with the material.</returns>
		/*
		RGBColor sampleTexture(double u, double v) const
		{
			if(m_material->hasTexture() && m_textureCoordinates[0]!=NULL)
			{
				RGBColor texel = m_material->getTexture().pixel(interpolateTextureCoordinate(u, v));
				return texel;
			}
			return RGBColor(1.0f, 1.0f, 1.0f);
		}
		*/

		/// <summary>
		/// Samples the triangle given the u,v coordinates of an intersection
		/// </summary>
		/// <param name="u"></param>
		/// <param name="v"></param>
		/// <returns></returns>
		Math::Vector3f samplePoint(double u, double v) const
		{
			return m_uAxis*u + m_vAxis*v + m_vertex0;
		}

		/// <summary>
		/// Samples the normal given the u,v coordinates of an intersection. The normal is oriented toward the 'toward' point.
		/// </summary>
		/// <param name="u">The u coordinate.</param>
		/// <param name="v">The v coordinate.</param>
		/// <param name="toward">The toward point.</param>
		/// <returns></returns>
		Math::Vector3f sampleNormal(double u, double v, const Math::Vector3f & toward) const
		{
			Math::Vector3f result = (m_vertexNormal[0]*(1 - u - v) + m_vertexNormal[1]*u + m_vertexNormal[2]*v).normalized();
			if ((toward - (m_vertex0+m_uAxis*u+m_vAxis*v))*result < 0.0)
			{
				return (result*(-1.0)).normalized();
			}
			return result.normalized();
		}

		Math::Vector3f sampleNormal(double u, double v, bool facing=true) const
		{
			Math::Vector3f result = (m_vertexNormal[0] * (1 - u - v) + m_vertexNormal[1] * u + m_vertexNormal[2] * v).normalized();
			if (!facing)
			{
				result = -result;
			}
			return result.normalized();
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	const Math::Vector3 & Triangle::uAxis() const
		///
		/// \brief	Gets the u axis.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	04/12/2013
		///
		/// \return	.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		const Math::Vector3f & uAxis() const		
		{ 
			return m_uAxis ;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	const Math::Vector3 & Triangle::vAxis() const
		///
		/// \brief	Gets the v axis.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	04/12/2013
		///
		/// \return	.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		const Math::Vector3f & vAxis() const
		{ 
			return m_vAxis ;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	const Math::Vector3 & Triangle::normal() const
		///
		/// \brief	Gets the normal.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	04/12/2013
		///
		/// \return	.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		const Math::Vector3f & normal() const
		{ return m_normal ; }

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	Math::Vector3 Triangle::normal(Math::Vector3 const & point) const
		///
		/// \brief	Gets the normal directed toward the half space containing the provided point.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	04/12/2013
		///
		/// \param	point	The point.
		///
		/// \return	.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		Math::Vector3f normal(Math::Vector3f const & point) const
		{
			if((point-m_vertex0)*m_normal<0.0)
			{ return m_normal*(-1.0) ; }
			return m_normal ; 
		}

		Math::Vector3f normal(bool facing) const
		{
			if (!facing)
			{
				return m_normal * (-1.0);
			}
			return m_normal;
		}

		/// <summary>
		/// Returns the direction of a reflected ray, from a surface normal and the direction of the incident ray.
		/// </summary>
		/// <param name="n">The n.</param>
		/// <param name="dir">The dir.</param>
		/// <returns></returns>
		static Math::Vector3f reflectionDirection(Math::Vector3f const & n, Math::Vector3f const & dir)
		{
			Math::Vector3f reflected(dir - n*(2.0f*(dir*n)));
			return reflected;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	Math::Vector3 Triangle::reflectionDirection(Math::Vector3 const & dir) const
		///
		/// \brief Returns the direction of a reflected ray, from the direction of the incident ray.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	04/12/2013
		///
		/// \param	dir	The direction of the incident ray.
		///
		/// \return	The direction of the reflected ray.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		Math::Vector3f reflectionDirection(Math::Vector3f const & dir) const
		{
			Math::Vector3f n = normal();
			Math::Vector3f reflected(dir-n*(2.0f*(dir*n))) ; 
			return reflected ;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	Math::Vector3 Triangle::reflectionDirection(Ray const & ray) const
		///
		/// \brief	Returns the direction of the reflected ray from a ray description.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	04/12/2013
		///
		/// \param	ray	The incident ray.
		///
		/// \return	.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		Math::Vector3f reflectionDirection(Ray const & ray) const
		{
			Math::Vector3f n = normal() ;
			//if(n*(ray.source()-vertex(0))<=0.0)
			if(n*(ray.source()-m_vertex0)<=0.0)
			{ n = n*(-1.0) ; }
			Math::Vector3f reflected(ray.direction()-n*(2.0f*(ray.direction()*n))) ; 
			return reflected ;
		}



		bool facing(const Math::Vector3f & to_view)const
		{
			return to_view * m_normal >= 0.0;
		}



		/*
		const * get_light_texture(const Math::Vector3f to_view)const
		{
			return get_light_texture(facing(to_view));
		}
		//*/

		size_t get_light_tex_index(size_t light_index, bool facing)const
		{
			assert(facing == 0 || facing == 1);
			return (light_index << 1) | facing;
		}




		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	bool Triangle::intersection(Ray const & r, double & t, double & u, double & v) const
		///
		/// \brief	Computes the intersection between a ray and this triangle.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	04/12/2013
		///
		/// \param	r		 	The tested ray.
		/// \param 	t	The distance between the ray source and the intersection point (t>=0).
		/// \param	u		 	The u coordinate of the intersection (useful for texture mapping).
		/// \param	v		 	The v coordinate of the intersection (useful for texture mapping).
		///
		/// \return	True if an intersection has been found, false otherwise.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		bool intersection(Ray const & r, double & t, double & u, double & v) const
		{
			/*
			double denum = m_normal * r.direction();
			double num = m_normal * (m_vertex0 - r.source());
			t = num / denum;
			if (t < 0.00000001)
			{
				return false;
			}
			Math::Vector3f inter_point = r.source() + r.direction() * t;
			Math::Vector3f point_to_a = inter_point - m_vertex0;
			v = (point_to_a * m_vAxis) / (m_vAxis.norm2());

			if (v < 0 || v > 1)
			{
				return false;
			}
			
			u = (point_to_a * m_uAxis) / (m_uAxis.norm2());
			
			if (u < 0 || u+v > 1)
			{
				return false;
			}
			
			return true;
			//*/
			

			/* find vectors for two edges sharing vert0 */
			const Math::Vector3f & edge1(uAxis()) ;
			const Math::Vector3f & edge2(vAxis()) ;

			/* begin calculating determinant - also used to calculate U parameter */
			Math::Vector3f pvec(r.direction() ^ edge2);

			/* if determinant is near zero, ray lies in plane of triangle */
			double det = edge1 * pvec ;
		
			if (fabs(det)<0.000000001) // if(det > -0.000001 && det < 0.000001) 
			{
				return false ; 
			}

			double inv_det = 1.0f / det;

			/* calculate distance from vert0 to ray origin */
			//Math::Vector3 tvec(r.source() - vertex(0));
			Math::Vector3f tvec(r.source() - m_vertex0);

			/* calculate U parameter and test bounds */
			u = (tvec * pvec) * inv_det;

			//std::cout<<"u = "<<u<<std::endl ;

			if (fabs(u-0.5)>0.5) //u < 0.0 || u > 1.0) //
			{
				return  false ;
			}

			/* prepare to test V parameter */
			Math::Vector3f qvec(tvec ^ edge1) ;

			/* calculate V parameter and test bounds */
			v = (r.direction() * qvec) * inv_det;
			if (v < 0.0 || u + v > 1.000000001)
			{
				return false ;
			}

			/* calculate t, ray intersects triangle */
			t = (edge2 * qvec) * inv_det;

			return t>=0.0001 ;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	bool Triangle::generalIntersection(Ray const & r, double & t, double & u, double & v) const
		///
		/// \brief	Computes the intersection between the ray and the plane supporting the triangle.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	04/12/2013
		///
		/// \param	r	The ray.
		/// \param	t	The distance between the ray source and the intersection.
		/// \param	u	The u coordinate of the intersection.
		/// \param	v	The v coordinate of the intersection.
		///
		/// \return	True if the ray is not parallel to the plane.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		/*
		bool generalIntersection(Ray const & r, double & t, double & u, double & v) const
		{
			// find vectors for two edges sharing vert0 
			const Math::Vector3f & edge1(uAxis()) ;
			const Math::Vector3f & edge2(vAxis()) ;
			double det,inv_det;

			// begin calculating determinant - also used to calculate U parameter 
			Math::Vector3f pvec(r.direction() ^ edge2);

			// if determinant is near zero, ray lies in plane of triangle 
			det = edge1 * pvec ;
		
			if (det > -0.000001 && det < 0.000001) 
			{
				return false ; 
			}

			inv_det = 1.0f / det;

			// calculate distance from vert0 to ray origin 
			Math::Vector3f tvec(r.source() - m_vertex0);

			// calculate U parameter and test bounds 
			u = (tvec * pvec) * inv_det;

			// prepare to test V parameter 
			Math::Vector3f qvec(tvec ^ edge1) ;

			// calculate V parameter and test bounds 
			v = (r.direction() * qvec) * inv_det;

			// calculate t, ray intersects triangle 
			t = (edge2 * qvec) * inv_det;

			return true ;
		}
		*/

		/// <summary>
		///  Returns the surface of the triangle
		/// </summary>
		/// <returns></returns>
		double surface() const
		{
			return (m_uAxis^m_vAxis).norm() / 2.0;
		}

		/// <summary>
		/// Computes random barycentric coordinates
		/// </summary>
		/// <returns></returns>
		Math::Vector3f randomBarycentric(Math::Sampler & sampler) const
		{
			double r = sampler.generateContinuous<double>();
			double s = sampler.generateContinuous<double>();
			double a = double(1.0) - sqrt(s);
			double b = (double)((1.0 - r)*sqrt(s));
			double c = r*sqrt(s);
			return Math::makeVector(a, b, c);
		}

		/// <summary>
		/// Computes a point on a triangle from the barycentric coordinates
		/// </summary>
		/// <param name="barycentric"> The barycentric coordinates of the point</param>
		/// <returns></returns>
		Math::Vector3f pointFromBraycentric(const Math::Vector3f & barycentric) const
		{
			const Math::Vector3f & tmp = barycentric;
			return ((*m_vertex[0]) * tmp[0] + (*m_vertex[1]) * tmp[1] + (*m_vertex[2]) * tmp[2]);
		}

		/// <summary>
		/// Samples the texture from the provided barycentric coordinates
		/// </summary>
		/// <param name="barycentic"> The barycentric coordinates of the point</param>
		/// <returns> The color of the texture at the given barycentric coordinates</returns>
		/*
		RGBColor sampleTexture(const Math::Vector3f & barycentic) const
		{
			if (m_material->hasTexture())
			{
				Math::Vector2f textureCoord = textureCoordinate(0)*barycentic[0] + textureCoordinate(1)*barycentic[1] + textureCoordinate(2)*barycentic[2];
				return m_material->getTexture().pixel(textureCoord);
			}
			return RGBColor(1.0, 1.0, 1.0);
		}
		*/

		/// <summary>
		///  Computes a random point on the triangle
		/// </summary>
		/// <returns></returns>
		Math::Vector3f randomPoint(Math::Sampler & sampler) const
		{
			/*
			double u = (double)rand() / RAND_MAX;
			double v = (double)rand() / RAND_MAX;
			
			if (u + v > 1)
			{
				u = 1 - u;
				v = 1 - v;
			}

			return samplePoint(u, v);
			*/

			

			double r = sampler.generateContinuous<double>();
			double s = sampler.generateContinuous<double>();
			double a = double(1.0) - sqrt(s);
			double b = (double)((1.0 - r)*sqrt(s));
			double c = r*sqrt(s);
			return ((*m_vertex[0]) * a + (*m_vertex[1]) * b + (*m_vertex[2]) * c);

			
		}



		Math::Vector2f randomUV(Math::Sampler & sampler)const
		{
			double u = sampler.generateContinuous<double>();
			double v = sampler.generateContinuous<double>();

			if (u + v > 1)
			{
				u = 1 - u;
				v = 1 - v;
			}
			return { u, v };
		}

		/// <summary>
		/// Computes the distance between the point and the plane on which the triangle lies.
		/// </summary>
		/// <param name="point"></param>
		/// <returns></returns>
		double planeDistance(const Math::Vector3f & point) const
		{
			return fabs((point - m_vertex0)*m_normal);
		}
		
	} ;
}




#endif
 
