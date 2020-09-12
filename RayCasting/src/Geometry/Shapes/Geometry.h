#ifndef _Geometry_Geometry_H
#define _Geometry_Geometry_H

//#include <Geometry/CastedRay.h>
#include <Geometry/Triangle.h>
#include <Geometry/ComputeVertexNormals.h>
#include <Geometry/Materials/Material.h>
#include <Math/Quaternion.h>
#include <Math/Vectorf.h>
#include <vector>
#include <deque>
#include <map>
#include <System/aligned_allocator.h>
#include <Math/Constant.h>
#include <Geometry/GeometryBase.h>
#include <Geometry/BoundingBox.h>

namespace Geometry
{

	class GeometryCollection : public GeometryBase
	{
	public:

		using TriangleCollection = ::std::vector<Triangle>;

	protected:
	    /// \brief	The vertices.
	    std::deque<Math::Vector3f> m_vertices ;
		/// \brief The texture coordinates
		std::deque<Math::Vector2f> m_textureCoordinates;
		/// \brief	The triangles.
		TriangleCollection m_triangles ;

		
		::std::vector<double> m_surface_sum;
		double m_current_sum;

		BoundingBox m_box;


		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	void GeometryCollection::updateTriangles()
		///
		/// \brief	Updates all the triangles of the geometry (normals, u and v vectors). This method should
		/// 		be called if some transformations arer applied on the vertices of the geometry.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	04/12/2013
		////////////////////////////////////////////////////////////////////////////////////////////////////
		void updateTriangles()
		{
			for(int cpt=0 ; cpt<(int)m_triangles.size() ; cpt++)
			{
				m_triangles[cpt].update() ;
			}
		}

	public:
		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	const std::deque<Math::Vector3> & GeometryCollection::getVertices() const
		///
		/// \brief	Gets the vertices.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	04/12/2013
		///
		/// \return	The vertices.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		const std::deque<Math::Vector3f> & getVertices() const
		{ return m_vertices ; }

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	const std::vector<Triangle> & GeometryCollection::getTriangles() const
		///
		/// \brief	Gets the triangles.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	04/12/2013
		///
		/// \return	The triangles.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		const TriangleCollection& getTriangles() const
		{ return m_triangles ; }


		TriangleCollection & getTriangles()
		{
			return m_triangles;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	GeometryCollection::GeometryCollection()
		///
		/// \brief	Default constructor.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	04/12/2013
		////////////////////////////////////////////////////////////////////////////////////////////////////
		GeometryCollection(Material * material=nullptr):
			GeometryBase(material),
			m_current_sum(0)
		{}

		/// <summary>
		/// Copy constructor.
		/// </summary>
		/// <param name="geometry">The geometry.</param>
		GeometryCollection(const GeometryCollection & geometry):
			GeometryBase(geometry),
			m_current_sum(geometry.m_current_sum)
		{
			merge(geometry);
		}




		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	unsigned int GeometryCollection::addVertex(const Math::Vector3 & vertex)
		///
		/// \brief	Adds a vertex.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	04/12/2013
		///
		/// \param	vertex	The vertex.
		///
		/// \return	The index of the added vertex.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		unsigned int addVertex(const Math::Vector3f & vertex)
		{ 
			m_vertices.push_back(vertex) ; 
			m_box.update(vertex);
			return m_vertices.size()-1 ;
		}

		/// <summary>
		/// Adds the texture coordinates.
		/// </summary>
		/// <param name="coord">The coordinates.</param>
		/// <returns></returns>
		unsigned int addTextureCoordinates(const Math::Vector2f & coord)
		{
			m_textureCoordinates.push_back(coord);
			return m_textureCoordinates.size() - 1;
		}

		/// <summary>
		/// Adds a triangle.
		/// </summary>
		/// <param name="i1">Zero-based index of the first vertex.</param>
		/// <param name="i2">Zero-based index of the second vertex.</param>
		/// <param name="i3">Zero-based index of the third vertex.</param>
		/// <param name="material">The material.</param>
		/// <param name="normals">The normals.</param>
		void addTriangle(int i1, int i2, int i3, const Math::Vector3f * normals = nullptr)
		{
			if(m_textureCoordinates.empty())
			{
				m_triangles.push_back(Triangle(&m_vertices[i1], &m_vertices[i2], &m_vertices[i3], this, normals));
			}
			else
			{
				m_triangles.push_back(Triangle(&m_vertices[i1], &m_vertices[i2], &m_vertices[i3], &m_textureCoordinates[i1], &m_textureCoordinates[i2], &m_textureCoordinates[i3], this, normals));
			}
		}

		/// <summary>
		/// Adds the triangle.
		/// </summary>
		/// <param name="i1">The index of the first vertex.</param>
		/// <param name="i2">The index of the second vertex.</param>
		/// <param name="i3">The index of the third vertex.</param>
		/// <param name="m_begin">The index of the first texture coordinates.</param>
		/// <param name="m_current">The index of the second texture coordinates.</param>
		/// <param name="t3">The index of the third texture coordinates.</param>
		/// <param name="material">The material.</param>
		/// <param name="normals">The vertices normals if any, nullptr means that normals must be  set to the face normal.</param>
		void addTriangle(int i1, int i2, int i3, int t1, int t2, int t3, const Math::Vector3f * normals = nullptr)
		{
			m_triangles.push_back(Triangle(&m_vertices[i1], &m_vertices[i2], &m_vertices[i3], &m_textureCoordinates[t1], &m_textureCoordinates[t2], &m_textureCoordinates[t3], this, normals));
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	void GeometryCollection::addTriangle(const Triangle & triangle)
		///
		/// \brief	Adds a triangle.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	04/12/2013
		///
		/// \param	triangle	The triangle.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		void addTriangle(const Triangle & triangle)
		{
			int i1 = addVertex(triangle.vertex(0)) ;
			int i2 = addVertex(triangle.vertex(1)) ;
			int i3 = addVertex(triangle.vertex(2)) ;
			addTriangle(i1,i2,i3, triangle.getVertexNormals()) ;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	void GeometryCollection::addTriangle(const Math::Vector3 & p0, const Math::Vector3 & p1,
		/// 	const Math::Vector3 & p2, Material * material)
		///
		/// \brief	Adds a triangle.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	04/12/2013
		///
		/// \param	p0					The first vertex.
		/// \param	p1					The second vertex.
		/// \param	p2					The third vertex.
		/// \param [in,out]	material	If non-null, the material.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		void addTriangle(const Math::Vector3f & p0, const Math::Vector3f & p1, const Math::Vector3f & p2, const Math::Vector3f * normals = nullptr)
		{
			int i1 = addVertex(p0) ;
			int i2 = addVertex(p1) ;
			int i3 = addVertex(p2) ;
			addTriangle(i1,i2,i3, normals) ;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	void GeometryCollection::merge(const GeometryCollection & geometry)
		///
		/// \brief	Merges the provided geometry with this one.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	04/12/2013
		///
		/// \param	geometry The geometry that should be merged
		////////////////////////////////////////////////////////////////////////////////////////////////////
		void merge(const GeometryCollection & geometry)
		{
			m_box.update(geometry.m_box);
			::std::map<const Math::Vector3f*, int> vectorToIndex ;
			for(size_t cpt=0 ; cpt<geometry.getVertices().size() ; cpt++)
			{
				const Math::Vector3f * vertex = &geometry.getVertices()[cpt];
				auto it = vectorToIndex.find(vertex) ;
				if(it==vectorToIndex.end())
				{
					int index = addVertex(geometry.getVertices()[cpt]) ;
					vectorToIndex.insert(::std::make_pair(vertex, index)) ;
				}
			}
			::std::map<const Math::Vector2f*, int> textureToIndex;
			for (size_t cpt = 0; cpt < geometry.m_textureCoordinates.size(); ++cpt)
			{
				const Math::Vector2f * vertex = &geometry.m_textureCoordinates[cpt];
				auto it = textureToIndex.find(vertex);
				if (it == textureToIndex.end())
				{
					int index = addTextureCoordinates(geometry.m_textureCoordinates[cpt]);
					textureToIndex.insert(::std::make_pair(vertex, index));
				}
			}

			for(int cpt=0 ; cpt<(int)geometry.getTriangles().size() ; cpt++)
			{
				const Triangle & current = geometry.getTriangles()[cpt] ;
				auto it1 = vectorToIndex.find(&current.vertex(0)) ;
				auto it2 = vectorToIndex.find(&current.vertex(1)) ;
				auto it3 = vectorToIndex.find(&current.vertex(2)) ;
				auto tex1 = textureToIndex.find(&current.textureCoordinate(0));
				auto tex2 = textureToIndex.find(&current.textureCoordinate(1));
				auto tex3 = textureToIndex.find(&current.textureCoordinate(2));
				if(tex1==textureToIndex.end())
				{
					addTriangle(it1->second, it2->second, it3->second, current.getVertexNormals());
				}
				else
				{
					addTriangle(it1->second, it2->second, it3->second, tex1->second, tex2->second, tex3->second);
				}
			}
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	bool GeometryCollection::intersection(CastedRay & ray)
		///
		/// \brief	Computes the intersection between this geometry and the provided Ray. This algorithm
		///			is straightforward... 
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	04/12/2013
		///
		/// \param [in,out]	ray	The ray.
		///
		/// \return	true if it succeeds, false if it fails.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		//bool intersection(CastedRay<Triangle> & ray)
		//{
		//	for(auto it=m_triangles.begin(), end=m_triangles.end() ; it!=end ; ++it)
		//	{
		//		ray.intersect((*it)) ;
		//	}
		//	return ray.validIntersectionFound() ;
		//}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	void GeometryCollection::translate(Math::Vector3 const & t)
		///
		/// \brief	Translate this geometry.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	04/12/2013
		///
		/// \param	t	The translation vector.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		void translate(Math::Vector3f const & t)
		{
			for(auto it=m_vertices.begin(), end=m_vertices.end() ; it!=end ; ++it)
			{
				(*it) = (*it)+t ; 
			}
			updateTriangles();
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	void GeometryCollection::scale(double v)
		///
		/// \brief	Applies a scale factor on this geometry.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	04/12/2013
		///
		/// \param	v	The scale factor.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		void scale(double v)
		{
			for(auto it=m_vertices.begin(), end=m_vertices.end() ; it!=end ; ++it)
			{
				(*it) = (*it)*v ;
			}
			updateTriangles() ;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	void GeometryCollection::scaleX(double v)
		///
		/// \brief	Scales geometry on X axis.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	04/12/2013
		///
		/// \param	v	The scale factor on X axis.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		void scaleX(double v)
		{
			for(auto it=m_vertices.begin(), end=m_vertices.end() ; it!=end ; ++it)
			{
				(*it)[0] = (*it)[0]*v ;
			}
			updateTriangles() ;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	void GeometryCollection::scaleY(double v)
		///
		/// \brief	Scales geometry on Y axis.
		/// 		
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	04/12/2013
		///
		/// \param	v	The scale factor.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		void scaleY(double v)
		{
			for(auto it=m_vertices.begin(), end=m_vertices.end() ; it!=end ; ++it)
			{
				(*it)[1] = (*it)[1]*v ;
			}
			updateTriangles() ;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	void GeometryCollection::scaleZ(double v)
		///
		/// \brief	Scales the geometry on Z axis.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	04/12/2013
		///
		/// \param	v	The scale factor.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		void scaleZ(double v)
		{
			for(auto it=m_vertices.begin(), end=m_vertices.end() ; it!=end ; ++it)
			{
				(*it)[2] = (*it)[2]*v ;
			}
			updateTriangles() ;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	void GeometryCollection::rotate(Math::Quaternion const & q)
		///
		/// \brief	Rotates this geometry.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	04/12/2013
		///
		/// \param	q	Quaternion describing the rotation.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		void rotate(Math::Quaternion<double> const & q)
		{
			for(auto it=m_vertices.begin(), end=m_vertices.end() ; it!=end ; ++it)
			{
				(*it) = q.rotate(*it).v() ;
			}
			updateTriangles() ;
		}

		/// <summary>
		/// Computes the per vertex normals.
		/// </summary>
		/// <param name="angle">The angle limit for smoothing surface.</param>
		void computeVertexNormals(double angle)
		{
			double cosAngleLimit = cos(angle);
			::std::vector<Triangle *> triangles;
			for (auto it = m_triangles.begin(), end = m_triangles.end(); it != end; ++it)
			{
				triangles.push_back(&(*it));
			}
			ComputeVertexNormals normalsComputation(triangles);
			normalsComputation.compute(cosAngleLimit);
		}


		BoundingBox box()const override
		{
			return m_box;
		}

		void recomputeBox()
		{
			m_box = BoundingBox();
			for (Triangle const& t : m_triangles)
			{
				m_box.update(t.box());
			}
		}

		virtual bool build_lights()
		{
			if (!m_material->is_emissive())
			{
				//no need to build the lights
				//maybe use an assert
				return false;
			}
			else
			{
				m_current_sum = 0;
				m_surface_sum.clear();

				for (Triangle const& tri : m_triangles)
				{
					m_current_sum += tri.surface();
					m_surface_sum.push_back(m_current_sum);
				}
				
				m_surface = m_current_sum;
				return true;
			}
		}



	public:

		virtual void sampleLight(SurfaceSample& res, Math::Sampler& sampler, unsigned int i=0)const override
		{
			unsigned int offset = (m_offset + i) % m_divisions;
			double rnd = sampler.generateContinuous<double>(offset, offset + 1) / double(m_divisions);
			rnd *= m_current_sum;
			auto found = std::lower_bound(m_surface_sum.cbegin(), m_surface_sum.cend(), rnd);
			size_t index = found - m_surface_sum.cbegin();
			const Triangle& tri = m_triangles[index];
			Math::Vector2f uv = tri.randomUV(sampler);
			Math::Vector3f point = tri.samplePoint(uv[0], uv[1]);
			
			double pdf = 1.0 / surface();
			
			//TODO maybe use the interpolated normal
			res = { pdf, this, &tri, tri.interpolateTextureCoordinate(uv[0], uv[1]), tri.normal(), point };
		}

		virtual void sampleLight(SurfaceSample& res, double xi1, double xi2)const override
		{
			xi1 *= m_current_sum;
			auto found = std::lower_bound(m_surface_sum.cbegin(), m_surface_sum.cend(), xi1);
			size_t index = found - m_surface_sum.cbegin();
			xi1 = (xi1 - index) / m_surface_sum.size();
			const Triangle& tri = m_triangles[index];
			Math::Vector2f uv = tri.sampleUV(xi1, xi2);
			Math::Vector3f point = tri.samplePoint(uv[0], uv[1]);

			double pdf = 1.0 / surface();

			//TODO maybe use the interpolated normal
			res = { pdf, this, &tri, tri.interpolateTextureCoordinate(uv[0], uv[1]), tri.normal(), point };
		}


		virtual void divide(unsigned int div = 1)
		{
			m_divisions = div;
		}


	} ;

}


#endif
