#ifndef _Geometry_Disk_H
#define _Geometry_Disk_H

#include <Geometry/Geometry.h>
#include <Geometry/Material.h>
#include <Geometry/Primitive.h>

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

namespace Geometry
{


	class Disk : public GeometryBase, public Primitive
	{
	protected:
		
		Math::Vector3f m_center;
		//these vectors are always normalized
		Math::Vector3f m_normal;
		//aligned with the x axis
		Math::Vector3f m_tg;
		//aligned with the y axis
		Math::Vector3f m_ctg;
		
		double m_radius, m_radius2;

		unsigned int m_angle_div=1, m_radius_div=1;

	public:

		Disk(Material* mat, Math::Vector3f center = 0.0, double radius = 1.0, Math::Vector3f const& normal = { 0, 0, 1 }) :
			GeometryBase(mat),
			m_center(center),
			m_normal(normal.normalized()),
			m_radius(radius),
			m_radius2(radius*radius)
		{
			m_surface = Math::pi * m_radius2;
			m_tg = Math::Vector3f(1, 0, 0);
			m_tg = (m_tg - normal) * (normal * m_tg);
			std::cout << m_tg << std::endl;
			if (m_tg.norm() < std::numeric_limits<double>::epsilon() * 10)
			{
				m_tg = { 0, 1, 0 };
				m_tg = (m_tg - normal) * (normal * m_tg);
				std::cout << m_tg << std::endl; 
				if (m_tg.norm() < std::numeric_limits<double>::epsilon() * 10)
				{
					m_tg = { 0, 0, 1 };
					m_tg = (m_tg - normal) * (normal * m_tg);
					std::cout << m_tg << std::endl;
					assert(m_tg.norm() > std::numeric_limits<double>::epsilon() * 10);
				}
			}
			m_tg = m_tg.normalized();
			m_ctg = normal ^ m_tg;
			m_ctg = m_ctg.normalized();
			m_tg = m_ctg ^ m_normal;
			m_tg = m_tg.normalized();
		}

		BoundingBox box()const final override
		{
			assert(m_tg[1] == 0);
			assert(m_ctg[0] == 0);
			BoundingBox res = { -m_normal*m_radius, m_normal*m_radius };
			res.update(-m_tg * m_radius);
			res.update(m_tg * m_radius);
			res.update(-m_ctg * m_radius);
			res.update(m_ctg * m_radius);
			return BoundingBox(m_center + res.min(), m_center + res.max());
		}

		Math::Vector3f center()const final override
		{
			return m_center;
		}

		const Math::Vector3f& normal()const
		{
			return m_normal;
		}

		Math::Vector3f normal(bool facing)const
		{
			if (facing)	return m_normal;
			return -m_normal;
		}

		const Math::Vector3f& tg()const
		{
			return m_tg;
		}

		const Math::Vector3f& ctg()const
		{
			return m_ctg;
		}

		double radius()const
		{
			return m_radius;
		}

		double radius2()const
		{
			return m_radius2;
		}

		bool facing(Math::Vector3f const& dir)const
		{
			return dir * m_normal > 0;
		}

		virtual bool build_lights()final override
		{
			return true;
		}

		virtual void divide(unsigned int div = 1.0)final override
		{
			//todo
		}

		virtual void sampleLight(SurfaceLightSample& res, Math::Sampler& sampler, unsigned int i = 0)const
		{

		}

		virtual void sampleLights(LightSampleStack& res, Math::Sampler& sampler, unsigned int n = 0)const
		{

		}


		template <class out_t>
		out_t& printAxis(out_t& out)const
		{
			out << "normal: "<<m_normal << std::endl;
			out << "tg: " << m_tg << std::endl;
			out << "ctg: " << m_ctg << std::endl;
			return out;
		}
	};




	class DiskFabrice : public GeometryCollection
	{
	protected:
	public:


		DiskFabrice(int nbDiv, Material * material):
			GeometryCollection(material)
		{
			unsigned int center = addVertex(Math::Vector3f()) ;
			::std::vector<unsigned int> vertices ;
			for(int cpt=0 ; cpt<nbDiv ; cpt++)
			{
				double angle = double((2.0f*M_PI/nbDiv)*cpt) ;
				int i = addVertex(Math::makeVector(cos(angle), sin(angle), 0.0)) ;
				vertices.push_back(i) ;
			}
			for(int cpt=0 ; cpt<nbDiv ; cpt++)
			{
				addTriangle(center, vertices[cpt], vertices[(cpt+1)%nbDiv]) ;
			}

			//Math::Vector3 * center = new Math::Vector3();
			//// Generation des points
			//for(int cpt=0 ; cpt<nbDiv ; cpt++)
			//{
			//	double angle = double((2.0f*M_PI/nbDiv)*cpt) ;
			//	addVertex(new Math::Vector3(cos(angle), sin(angle), 0.0)) ;
			//}
			//addVertex(center) ;
			// Generation des triangles
			//for(int cpt=0 ; cpt<nbDiv ; cpt++)
			//{
			//	Math::Vector3 * pt1 = m_vertices[(cpt)%(nbDiv)] ;
			//	Math::Vector3 * pt2 = m_vertices[(cpt+1)%(nbDiv)] ;
			//	addTriangle(new Triangle(center, pt1, pt2, material)) ;
			//}
		}
	} ;
} ;

#endif
