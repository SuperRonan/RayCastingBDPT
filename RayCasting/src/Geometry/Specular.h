#pragma once

#include <Geometry/Material.h>
#include <Math/Constant.h>


namespace Geometry
{

	class Specular :public Material
	{
	protected:

		RGBColor m_specular;
		double m_shininess;


		bool should_use_direct(double s)
		{
			return false;
			return s <= 1;
		}

	public:

		virtual RGBColor ID_COLOR()const
		{
			return SPECULAR_ID_COLOR;
		}


		Specular(RGBColor const& specular, double shininess=1, RGBColor const& emissive = 0, std::string const& texture = "") :
			Material(emissive, should_use_direct(shininess), false, texture),
			m_specular(specular * (shininess+1) / Math::twoPi),
			m_shininess(shininess)
		{
			//m_albedo = specular.grey();
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
			double coswi = (hit.primitive_normal * wi);
			if (coswi <= 0)	return 0;
			double cosr = hit.primitive_normal.reflect(wo) * wi;
			if (cosr <= 0)	return 0;
			return m_specular * pow(cosr, m_shininess);
		}

	public:


		virtual void sampleBSDF(Hit const& hit, unsigned int diffuse_samples, unsigned int specular_samples, DirectionSample& out, Math::Sampler& sampler)const override
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
			double cosd = hit.primitive_normal * out.direction;
			double bsdf = pow(cosr, m_shininess);
			
			out.pdf = (m_shininess + 1) * bsdf / Math::twoPi;

			// no transmitance
			if (cosd <= 0)
			{
				bsdf = 0;
			}
			
			
			out.bsdf = m_specular * bsdf;
		}

		virtual RGBColor BSDF(Hit const& hit, Math::Vector3f const& wi, Math::Vector3f const& wo)const override
		{
			return _BSDF(hit, wi, wo);
		}


		virtual void sampleBSDF(Hit const& hit, unsigned int diffuse_samples, unsigned int specular_samples, DirectionStack& res, Math::Sampler& sampler)const override
		{
			if (!m_specular.isBlack())
			{
				Math::Vector3f reflected = hit.primitive_reflected();
				Math::RandomDirection direction_sampler(&sampler, reflected, m_shininess);
				for (unsigned int i = 0; i < specular_samples; ++i)
				{
					DirectionSample out;
					out.direction = direction_sampler.generate();
					//TODO transparacy
					double cosr = reflected * out.direction;
					double cosd = hit.normal * out.direction;
					if (cosr <= 0)
					{
						cosr = -cosr;
						out.direction = -out.direction;
					}
					if (cosd <= 0)
					{
						//the bsdf of the sample will be 0 anyway,
						continue;
					}
					double bsdf = pow(cosr, m_shininess);
					out.pdf = (m_shininess + 1) * bsdf / Math::twoPi;
					out.bsdf = m_specular * bsdf;
					//out.weight = 1;
					res.push(out);
				}
			}
		}


		virtual void BSDF(Hit const& hit, ColorStack& res, LightSampleStack const& wis)const override
		{
			for (SurfaceLightSample const& wi : wis)
			{
				res.push(BSDF(hit, wi.vector, hit.to_view));
			}
		}



		virtual double pdf(Hit const& hit, Math::Vector3f const& wi, Math::Vector3f const& wo)const override
		{
			double c = hit.primitive_normal.reflect(wo) * wi;
			if (c < 0)
				return 0;
			return (m_shininess + 1) / Math::twoPi * pow(c, m_shininess) ;
		}



	};
}