#pragma once

#include <Geometry/Materials/Material.h>
#include <Math/Constant.h>


namespace Geometry
{

	class Glossy :public Material
	{
	protected:

		RGBColor m_specular;
		double m_shininess;


		static bool is_spicky(double s)
		{
			return s > 50;
		}

	public:

		virtual RGBColor ID_COLOR()const
		{
			return SPECULAR_ID_COLOR;
		}


		Glossy(RGBColor const& specular, double shininess=1, RGBColor const& emissive = 0, std::string const& texture = "") :
			Material(emissive, texture),
			m_specular(specular * (shininess+1) / Math::twoPi),
			m_shininess(shininess)
		{
			m_delta = false;
			m_spicky = is_spicky(m_shininess);
		}



		void setSpecular(RGBColor const& color)
		{
			m_specular = color * (m_shininess + 1) / Math::twoPi;
			//m_albedo = color.grey();
		}

		const RGBColor& getSpecular() const
		{
			return m_specular;
		}

		void setShininess(double s)
		{
			m_specular = m_specular * Math::twoPi / (m_shininess + 1);
			m_shininess = s;
			m_specular = m_specular * (m_shininess + 1) / Math::twoPi;
		}

		const double& getShininess() const
		{
			return m_shininess;
		}

	protected:

		__forceinline RGBColor _BSDF(Hit const& hit, Math::Vector3f const& wi, Math::Vector3f const& wo)const
		{
			double coswi = (hit.orientedPrimitiveNormal() * wi);
			if (coswi <= 0)	return 0;
			double cosr = hit.primitive_normal.reflect(wo) * wi;
			if (cosr <= 0)	return 0;
			return m_specular * pow(cosr, m_shininess);
		}

	public:


		virtual void sampleBSDF(Hit const& hit, DirectionSample& out, Math::Sampler& sampler, bool RADIANCE=true, double *pdf_rev=nullptr)const override
		{
			Math::Vector3f reflected = hit.primitive_reflected();
			Math::RandomDirection direction_sampler(&sampler, reflected, m_shininess);
			out.direction = direction_sampler.generate();

			double cosr = reflected * out.direction;
			//the sampler can generate on the wrong hemisphere
			if (cosr <= 0)
			{
				cosr = -cosr;
				out.direction = -out.direction;
			}
			double cosd = hit.orientedPrimitiveNormal() * out.direction;
			double bsdf = pow(cosr, m_shininess);
			
			out.pdf = (m_shininess + 1) * bsdf / Math::twoPi;

			// no transmitance
			if (cosd <= 0)
			{
				bsdf = 0;
			}
			
			
			out.bsdf = m_specular * bsdf;
			if (pdf_rev)
				*pdf_rev = out.pdf;
		}

		virtual RGBColor BSDF(Hit const& hit, Math::Vector3f const& wi, Math::Vector3f const& wo, bool RADIANCE = false)const override
		{
			return _BSDF(hit, wi, wo);
		}


		virtual double pdf(Hit const& hit, Math::Vector3f const& wi, Math::Vector3f const& wo, bool RADIANCE=false)const override
		{
			double c = hit.primitive_normal.reflect(wo) * wi;
			if (c < 0)
				return 0;
			return (m_shininess + 1) / Math::twoPi * pow(c, m_shininess) ;
		}



	};
}