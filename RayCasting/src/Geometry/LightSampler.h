#ifndef _Geometry_LightSampler_H
#define _Geometry_LightSampler_H

#include <random>
#include <Geometry/Triangle.h>
#include <Geometry/Geometry.h>
#include <Geometry/PointLight.h>
#include "Phong.h"

namespace Geometry
{
	/// <summary>
	/// A light sampler. To initialize this sampler, just add triangles or geometries. The sampler
	/// will keep triangles with a material that emits light.
	/// </summary>
	/*
	class LightSampler
	{
	protected:
		/// <summary> A random number generator. </summary>
		mutable std::mt19937_64 randomGenerator;
		::std::vector<double> surfaceSum;
		::std::vector<const Triangle *> allTriangles;
		double currentSum ;

	public:
		/// <summary>
		/// Constructor
		/// </summary>
		LightSampler()
			: currentSum(0.0)
		{}

		/// <summary>
		/// Adds a triangle to the light sampler.
		/// </summary>
		/// <param name="triangle"></param>
		void add(const Triangle & triangle)
		{
			if (!triangle.material()->getEmissive().isBlack())
			{
				currentSum += triangle.surface();
				surfaceSum.push_back(currentSum);
				allTriangles.push_back(&triangle);
			}
		}

		/// <summary>
		/// Adds a geometry to the light sampler.
		/// </summary>
		/// <param name="geometry"></param>
		void add(const GeometryCollection & geometry)
		{
			auto & triangles = geometry.getTriangles();
			add(triangles.begin(), triangles.end());
		}

		/// <summary>
		/// Adds a range of triangles or geometries to the light sampler.
		/// </summary>
		template <class iterator>
		void add(iterator begin, iterator end)
		{
			for (iterator it = begin; it != end; ++it)
			{
				add(*it);
			}
		}



		///
		/// params are out
		///
		///
		void generate(const Triangle * & tri, PointLight & l)const
		{
			double random = double(randomGenerator()) / double(randomGenerator.max());
			random = random * currentSum;
			auto found = std::lower_bound(surfaceSum.begin(), surfaceSum.end(), random);
			size_t index = found - surfaceSum.begin();
			tri = allTriangles[index];
			Math::Vector3f barycentric = tri->randomBarycentric();
			l = PointLight(tri->pointFromBraycentric(barycentric), tri->sampleTexture(barycentric) * (static_cast<const Phong *>(tri->material())->getEmissive())); //* emisive?
		}

		/// <summary>
		/// Generates nb point lights by sampling the kept triangles.
		/// </summary>
		template <class iterator>
		iterator generate(iterator output, size_t nb)
		{
			for (size_t cpt = 0; cpt < nb; ++cpt)
			{
				(*output) = generate();
				++output;
			}
			return output;
		}

		/// <summary>
		/// Generates end-begin point lights by sampling the kept triangles.
		/// </summary>
		template <class iterator>
		iterator generate(iterator begin, iterator end)
		{
			for (; begin != end; ++begin)
			{
				(*begin) = generate();
				++begin;
			}
			return begin;
		}

		/// <summary>
		/// Returns true if the light sampler can sample lights.
		/// </summary>
		/// <returns></returns>
		bool hasLights() const
		{
			return allTriangles.size() > 0;
		}
	};
	*/
}

#endif