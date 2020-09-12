#pragma once

#include <random>
#include <Geometry/Materials/Material.h>
#include <Math/Sampler.h>
#include <Geometry/BoundingBox.h>
#include <Geometry/Sample.h>

namespace Geometry
{


	class GeometryBase
	{
	protected:

		unsigned int m_divisions=1;
		mutable unsigned int m_offset=0;

		Material * m_material;

		double m_surface;
		

	public:


		GeometryBase(GeometryBase const& other) = default;

		GeometryBase(Material* mat, unsigned int div = 1, unsigned int offset = 0) :
			m_divisions(div),
			m_material(mat),
			m_offset(offset)
		{}

		virtual bool build_lights() = 0;

		

		virtual void divide(unsigned int div=1) = 0;

		virtual void sampleLight(SurfaceSample& res, Math::Sampler& sampler, unsigned int i=0)const = 0;
		
		virtual void sampleLight(SurfaceSample& res, double xi1, double xi2)const = 0;

	public:

		virtual void sampleLight(SurfaceSample& res, const Hit& hit, Math::Sampler& sampler, unsigned int i = 0)const
		{
			sampleLight(res, sampler, i);
		}

		virtual void sampleLight(SurfaceSample& res, const Hit& hit, double xi1, double xi2)const
		{
			sampleLight(res, xi1, xi2);
		}


		double surface()const
		{
			return m_surface;
		}

		void set_material(Material* mat)
		{
			m_material = mat;
		}

		const Material* getMaterial()const
		{
			return m_material;
		}

		Material* getMaterial()
		{
			return m_material;
		}

		virtual void increment_offset(unsigned int i)const
		{
			m_offset += i;
			m_offset %= m_divisions;
		}
		
		////////////////////////
		//check if the reauested capacity div 
		////////////////////////
		void check_capacity(unsigned int div)
		{
			if (m_divisions != div)
			{
				std::cout << "Warning: actual number of divisions: " << m_divisions << ", instead of " << div << std::endl;
			}
		}

		void reset()const
		{
			m_offset = 0;
			
		}

		//the pdf of sampling one point uniformaly on the entire geometry
		double pdfSamplingPoint()const
		{
			return 1.0 / surface();
		}

		virtual double pdfSamplingPoint(Hit const& hit, Math::Vector3f const& point)const
		{
			return pdfSamplingPoint();
		}

		unsigned int offset()const
		{
			return m_offset;
		}

		unsigned int divisions()const
		{
			return m_divisions;
		}

		virtual BoundingBox box()const = 0;
	};

}