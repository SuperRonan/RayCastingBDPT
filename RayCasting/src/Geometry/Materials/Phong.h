#pragma once

#include "Material.h"

namespace Geometry
{
	// A BRDF material:
	// diffuse + glossy
	class Phong : public Material
	{
	protected:

		RGBColor m_diffuse, m_glossy;
		double m_shininess;

		double pmf_choose_diffuse()const
		{
			return m_diffuse.energy() / (m_glossy.energy() + m_diffuse.energy());
		}

	public:

		Phong(RGBColor const& dif = 1, RGBColor const& glossy = 0, double s = 1, RGBColor const& em = 0, std::string const& tex = "") :
			Material(em, tex),
			m_diffuse(dif),
			m_glossy(glossy),
			m_shininess(s)
		{

		}


		virtual RGBColor ID_COLOR()const override
		{
			return PHONG_ID_COLOR;
		}

		virtual void sampleBSDF(Hit const& hit, DirectionSample& out, Math::Sampler& sampler, bool RADIANCE=false, double* pdf_rev=nullptr)const override
		{
			const RGBColor tex = getTexturePixel(hit.tex_uv);
			Math::Vector3f normal = hit.orientedPrimitiveNormal();
			double cos_wo = hit.to_view * normal;
			const double pmf = pmf_choose_diffuse();
			const double xi = sampler.generateContinuous<double>();
			if (xi < pmf) // Sample the diffuse
			{
				Math::RandomDirection diffuse_sampler(&sampler, normal);
				out.direction = diffuse_sampler.generate();
			}
			else // Sample the glossy
			{
				Math::Vector3f reflected = normal.reflect(hit.to_view);
				Math::RandomDirection glossy_sampler(&sampler, reflected, m_shininess);
				out.direction = glossy_sampler.generate();
			}
			out.pdf = Phong::pdf(hit, out.direction, hit.to_view, RADIANCE);
			out.bsdf = Phong::BSDF(hit, out.direction, hit.to_view, RADIANCE);
		}


		virtual RGBColor BSDF(Hit const& hit, Math::Vector3f const& wi, Math::Vector3f const& wo, bool RADIANCE = false)const override
		{
			RGBColor res = 0;
			Math::Vector3f normal = hit.orientedPrimitiveNormal(); // should take wi ? 
			double cos_wi = wi * normal;
			double cos_wo = wo * normal;
			if (cos_wo * cos_wi > 0) // same hemisphere
			{
				assert(cos_wi > 0);
				const double diffuse_rho = 1.0 / Math::twoPi;

				Math::Vector3f reflected = hit.primitive_normal.reflect(wo);
				const double cos_r = reflected * wi;
				const double glossy_rho = cos_r > 0 ? (std::pow(cos_r, m_shininess) * (m_shininess + 1) / Math::twoPi / (RADIANCE ? cos_wo : cos_wi) ) : 0;

				RGBColor texture = getTexturePixel(hit.tex_uv);

				res = m_diffuse * texture * diffuse_rho + m_glossy * glossy_rho;
			}
			return res;
		}

		virtual double pdf(const Hit & hit, Math::Vector3f const& wi, Math::Vector3f const& wo, bool RADIANCE=false)const override
		{
			double res = 0;
			Math::Vector3f normal = hit.orientedPrimitiveNormal(); // should take wi ? 
			double cos_wi = wi * normal;
			double cos_wo = wo * normal;
			if (cos_wo * cos_wi > 0) // same hemisphere
			{
				assert(cos_wi > 0);
				const double pdf_diffuse = cos_wi / Math::pi;

				Math::Vector3f reflected = hit.primitive_normal.reflect(wo);
				const double cos_r = reflected * wi;
				const double pdf_glossy = cos_r > 0 ? (std::pow(cos_r, m_shininess) * (m_shininess + 1) / Math::twoPi) : 0;

				const double pmf = pmf_choose_diffuse();

				res = pmf * pdf_diffuse + (1 - pmf) * pdf_glossy;
			}
			else
			{ // maybe the wi is in the glossy lobe
				Math::Vector3f reflected = hit.primitive_normal.reflect(wo);
				const double cos_r = reflected * wi;
				const double pdf_glossy = cos_r > 0 ? (std::pow(cos_r, m_shininess) * (m_shininess + 1) / Math::twoPi) : 0;

				const double pmf = pmf_choose_diffuse();

				res = (1 - pmf) * pdf_glossy;
			}
			return res;
		}

	};
}