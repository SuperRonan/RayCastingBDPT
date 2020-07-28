#pragma once

#include <Geometry/Materials/Material.h>
#include <Math/Constant.h>


namespace Geometry
{
	// A bad implementation of a lambert material, with a pdf on both sides (so half the samples are wasted)
	// For an experiment
	class BadLambertian :public Material
	{
	protected:

		RGBColor m_diffuse;

	public:

		BadLambertian(RGBColor const& diffuse, RGBColor const& emissive = 0, std::string const& texture = "") :
			Material(emissive, texture),
			m_diffuse(diffuse / Math::pi)
		{
			m_delta = false;
			m_spicky = false;
			//m_albedo = diffuse.grey();
		}

		void setDiffuse(RGBColor const& color)
		{
			m_diffuse = color / Math::pi;
			//m_albedo = color.grey();
		}

		const RGBColor& getDiffuse() const
		{
			return m_diffuse;
		}

		virtual RGBColor ID_COLOR()const
		{
			return LAMBERT_ID_COLOR;
		}


		virtual void sampleBSDF(Hit const& hit, DirectionSample& out, Math::Sampler& sampler, bool RADIANCE = false, double* pdf_rev = nullptr)const override
		{
			RGBColor dif = m_diffuse * getTexturePixel(hit.tex_uv);
			Math::Vector3f normal = hit.orientedPrimitiveNormal();
			
			//There is not really a point of sampling a 
			Math::RandomDirection direction_sampler(normal);
			
			double u = sampler.generateContinuous<double>(-1, 1);
			double v = sampler.generateContinuous<double>();
			
			out.direction = direction_sampler.generate(abs(u), v);
			
			double cos_wi = (out.direction * normal);
			out.pdf = cos_wi / Math::twoPi;
			if (u > 0)
			{
				out.bsdf = dif;
			}
			else
			{
				out.direction = -out.direction;
				out.bsdf = 0;
			}
		}

		virtual RGBColor BSDF(Hit const& hit, Math::Vector3f const& wi, Math::Vector3f const& wo, bool RADIANCE = false)const override
		{
			RGBColor res = 0;
			if (((wi * hit.primitive_normal) * (wo * hit.primitive_normal)) > 0)
				res = m_diffuse * getTexturePixel(hit.tex_uv);
			return res;
		}



		virtual double pdf(Hit const& hit, Math::Vector3f const& wi, Math::Vector3f const& wo, bool RADIANCE = false)const override
		{
			double cos_wi = std::abs(wi * hit.primitive_normal);
			return cos_wi / Math::twoPi;
		}
	};
}