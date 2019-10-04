#ifndef _Geometry_ComputeVertexNormals_H
#define _Geometry_ComputeVertexNormals_H

#include <Geometry/Triangle.h>
#include <Math/Vectorf.h>
#include <vector>
#include <map>
#include <iterator>

namespace Geometry
{
	/// <summary>
	/// This class compute per vertex normals for a set of triangles.
	/// </summary>
	class ComputeVertexNormals
	{
	protected:
		/// \brief An edge composed of two vector3f
		typedef std::pair<const Math::Vector3f *, const Math::Vector3f *> Edge;

		/// \brief Creates an edge (the smallest pointer is the first, the greater is the second)
		Edge makeEdge(const Math::Vector3f* v1, const Math::Vector3f * v2)
		{
			return ::std::make_pair(::std::min(v1, v2), ::std::max(v1, v2));
		}

		/// \brief The set of triangles
		::std::vector<Triangle*> m_triangles;
		/// \brief A map between an edge and its associated triangles
		::std::map<Edge, ::std::vector<size_t>> m_edgeToTriangle;
		/// \brief The groups of faces
		::std::vector<::std::vector<size_t>> m_groups;

	public:
		/// <summary>
		/// Initializes a new instance of the <see cref="ComputeVertexNormals"/> class.
		/// </summary>
		/// <param name="triangles">The set of triangles for which per vertex normals should be computed.</param>
		ComputeVertexNormals(const ::std::vector<Triangle*> triangles)
			: m_triangles(triangles)
		{}

		/// <summary>
		/// Computes per vertex normals for the set of triangles.
		/// </summary>
		/// <param name="cosLimit">The cosine of the mximum angle allowing normal interpolation.</param>
		void compute(double cosLimit)
		{
			computeEdgeToTriangle();
			filterEdges(cosLimit);
			computeGroups();
			computeNormalsForGroups();
			ensure(cosLimit);
		}

	protected:
		/// <summary>
		/// Computes the per vertex normals for computed groups of triangles.
		/// </summary>
		void computeNormalsForGroups()
		{
			for (auto it = m_groups.begin(), end = m_groups.end(); it != end; ++it)
			{
				computeNormalsForGroup(*it);
			}
		}

		/// <summary>
		/// Computes the per vertex normals for a given groups of triangles.
		/// </summary>
		/// <param name="group">The group.</param>
		void computeNormalsForGroup(const std::vector<size_t> & group)
		{
			// We create a map accumulating the sum of the normals
			::std::map<const Math::Vector3f *, Math::Vector3f> m_normals;
			for (auto it = group.begin(), end = group.end(); it != end; ++it)
			{
				size_t current = *it;
				const Triangle & triangle = *m_triangles[current];
				for (int i = 0; i < 3; ++i)
				{
					m_normals[&triangle.vertex(i)] += triangle.normal();
				}
			}
			// We set the normals of the triangles belonging to the group
			for(auto it=group.begin(), end=group.end() ; it!=end ; ++it)
			{
				Triangle & triangle = *m_triangles[*it];//**it;
				for (size_t i = 0; i < 3; ++i)
				{
					auto found = m_normals.find(&triangle.vertex(i));
					if (found != m_normals.end())
					{
						triangle.setVertexNormal(i, found->second.normalized());
					}
				}
			}
		}

		/// <summary>
		/// Computes the groups of triangles based on the edge to triangle map.
		/// </summary>
		void computeGroups()
		{
			// 1- We compute the triangle adjency graph.
			::std::vector<::std::vector<size_t>> graph(m_triangles.size());
			for (auto it = m_edgeToTriangle.begin(), end = m_edgeToTriangle.end(); it != end; ++it)
			{
				const ::std::vector<size_t> & triangles = it->second;
				graph[triangles[0]].push_back(triangles[1]);
				graph[triangles[1]].push_back(triangles[0]);
			}

			// 2 - We compute the connected components.
			::std::vector<bool> explored(m_triangles.size(), false);
			::std::vector<size_t> toExplore;
			for (size_t cpt = 0; cpt < graph.size() ; ++cpt)
			{
				if (explored[cpt]) { continue; }
				//::std::cout<<"Computing group "<<m_groups.size()<<"/"<<m_triangles.size()<<"..."<<::std::flush;
				m_groups.push_back(std::vector<size_t>());
				toExplore.push_back(cpt);
				while (!toExplore.empty())
				{
					size_t current = toExplore.back();
					toExplore.pop_back();
					if (explored[current]) { continue; }
					explored[current] = true;
					m_groups.back().push_back(current);
					std::copy(graph[current].begin(), graph[current].end(), ::std::back_inserter(toExplore));
				}
				//::std::cout << "OK" << ::std::endl;
			}
		}

		/// <summary>
		/// Computes the edge to triangle map.
		/// </summary>
		void computeEdgeToTriangle()
		{
			for (size_t cpt=0 ; cpt<m_triangles.size() ; ++cpt)
			{
				const Triangle & triangle = *m_triangles[cpt];
				for(size_t i=0; i<3 ; ++i)
				{
					m_edgeToTriangle[makeEdge(&triangle.vertex(i), &triangle.vertex((i+1)%3))].push_back(cpt);
				}
			}
		}

		/// <summary>
		/// Filters the edge to triangle map by removing edges shared by triangles for which the dot product 
		/// of normals is lower than the given threshold.
		/// </summary>
		/// <param name="cosLimit">The cosine limit.</param>
		void filterEdges(double cosLimit)
		{
			::std::vector<::std::map<Edge, ::std::vector<size_t>>::iterator> toRemove;
			::std::vector<size_t> trianglesToRemove;
			for (auto it = m_edgeToTriangle.begin(), end = m_edgeToTriangle.end(); it != end; ++it)
			{
				const ::std::vector<size_t> & triangles = it->second;
				if (triangles.size() < 2 )
				{
					toRemove.push_back(it);
				}
				else if (triangles.size() == 2)
				{
					if (m_triangles[triangles[0]]->normal()*m_triangles[triangles[1]]->normal() < cosLimit)
					{
						toRemove.push_back(it);
					}
				}
				else
				{
					//::std::cerr << "3D object topology is strange... an edge is shared between " << triangles.size() << " faces. We skip this edge." << ::std::endl;
					toRemove.push_back(it);
					::std::copy(it->second.begin(), it->second.end(), ::std::back_inserter(trianglesToRemove));
				}
			}
			for (auto it = toRemove.begin(), end = toRemove.end(); it != end; ++it)
			{
				m_edgeToTriangle.erase(*it);
			}
			toRemove.clear();
			::std::sort(trianglesToRemove.begin(), trianglesToRemove.end());
			trianglesToRemove.erase(::std::unique(trianglesToRemove.begin(), trianglesToRemove.end()), trianglesToRemove.end());
			for (auto it = m_edgeToTriangle.begin(), end = m_edgeToTriangle.end(); it != end; ++it)
			{
				if (std::binary_search(trianglesToRemove.begin(), trianglesToRemove.end(), it->second[0]))
				{
					toRemove.push_back(it);
				}
				else if (std::binary_search(trianglesToRemove.begin(), trianglesToRemove.end(), it->second[1]))
				{
					toRemove.push_back(it);
				}
			}
			for (auto it = toRemove.begin(), end = toRemove.end(); it != end; ++it)
			{
				m_edgeToTriangle.erase(*it);
			}
		}

		/// <summary>
		/// Ensures that the cosine of the angle between normals of each triangle is greater than the given threshold.
		/// </summary>
		/// <param name="cosLimit">The cosine limit.</param>
		void ensure(double cosLimit)
		{
			for (auto it = m_triangles.begin(), end = m_triangles.end(); it != end; ++it)
			{
				Triangle & triangle = (**it);
				// If the angle between normals is greater that the threshold, we reset the per vertex normals to the triangle normal
				if (triangle.getVertexNormal(0)*triangle.getVertexNormal(1) < cosLimit ||
					triangle.getVertexNormal(1)*triangle.getVertexNormal(2) < cosLimit ||
					triangle.getVertexNormal(0)*triangle.getVertexNormal(2) < cosLimit)
				{
					triangle.setVertexNormal(0, triangle.normal());
					triangle.setVertexNormal(1, triangle.normal());
					triangle.setVertexNormal(2, triangle.normal());
				}
				// If the norm of a normal is 0 (it should not be...), we reset it to the triangle normal.
				if (triangle.getVertexNormal(0).norm() == 0.0) { triangle.setVertexNormal(0, triangle.normal()); }
				if (triangle.getVertexNormal(1).norm() == 0.0) { triangle.setVertexNormal(1, triangle.normal()); }
				if (triangle.getVertexNormal(2).norm() == 0.0) { triangle.setVertexNormal(2, triangle.normal()); }
			}
		}
	};
}

#endif