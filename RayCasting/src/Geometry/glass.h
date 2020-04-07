#pragma once

#include <Geometry/Material.h>
#include <Math/RandomDirection.h>

namespace Geometry
{

	class Glass : public Material
	{
	protected:

		RGBColor m_albedo;

		double m_eta;

	public:

		Glass(RGBColor const& albedo, double eta, std::string const& texture = "") :
			Material(0, texture),
			m_albedo(albedo),
			m_eta(eta)
		{
			m_delta = true;
			m_spicky = true;
		}

		virtual RGBColor ID_COLOR()const
		{
			return GLASS_ID_COLOR;
		}

		virtual void sampleBSDF(Hit const& hit, DirectionSample& out, Math::Sampler& sampler, bool RADIANCE=true)const override
		{
			bool entering = hit.facing;
			double n1 = 1, n2 = m_eta;
			if (!entering)	std::swap(n1, n2);

			Math::Vector3f normal = hit.orientedPrimitiveNormal();

			double cos_theta = std::abs(hit.to_view * normal);
			double sin_theta = std::sqrt(1 - cos_theta * cos_theta);
			out.bsdf = m_albedo;
			if (sin_theta == 0)
			{
				out.direction = -normal;
				out.pdf = 1;
			}
			else
			{
				Math::Vector3f tg = (hit.to_view - normal * cos_theta) / (sin_theta);
				double eta_ratio = n1 / n2;
				double next_sin_theta = eta_ratio * sin_theta;

				bool reflect = next_sin_theta >= 1;
				double fresnel_reflectance = 1;
				double pdf_reflect = 1;

				double next_cos_theta;
				
				if (!reflect)
				{
					next_cos_theta = std::sqrt(1 - next_sin_theta * next_sin_theta);
					double rs = (n2 * cos_theta - n1 * next_cos_theta) / (n2 * cos_theta + n1 * next_cos_theta);
					double rp = (n1 * cos_theta - n2 * next_cos_theta) / (n1 * cos_theta + n2 * next_cos_theta);

					fresnel_reflectance = 0.5 * (rs * rs + rp * rp);
					pdf_reflect = fresnel_reflectance;

					double xi = sampler.generateContinuous<double>();

					reflect = xi <= pdf_reflect;
				}


				if (reflect)
				{
					out.direction = hit.primitive_reflected();
					out.bsdf *= fresnel_reflectance;
					out.pdf = pdf_reflect;
				}
				else // transmit
				{
					out.direction = normal * (-next_cos_theta) + tg * (-next_sin_theta);
					out.bsdf *= (1 - fresnel_reflectance);
					/*
					if (n2 > n1)
					{
						out.bsdf *= ((n2 * n2) / (n1 * n1));
					}
					*/
					out.pdf = 1 - pdf_reflect;
				}
			}
			
			if (RADIANCE)
				out.bsdf *= 1.0 / std::abs(normal * hit.to_view);
			else
				out.bsdf *= 1.0 / std::abs(out.direction * normal);
			
		}

		virtual RGBColor BSDF(Hit const& hit, Math::Vector3f const& wi, Math::Vector3f const& wo, bool RADIANCE = false)const override
		{
			return 0;
		}

		virtual double pdf(Hit const& hit, Math::Vector3f const& wi, Math::Vector3f const& wo)const override
		{
			return 0;
		}

	};

	class GlossyGlass : public Material
	{
	protected:

		RGBColor m_albedo;

		double m_eta;

		double m_shininess;

	public:

		GlossyGlass(RGBColor const& albedo, double eta, double shininess, std::string const& texture = "") :
			Material(0, texture),
			m_albedo(albedo),
			m_eta(eta),
			m_shininess(shininess)
		{
			m_delta = false;
			m_spicky = shininess > 3;
		}

		virtual RGBColor ID_COLOR()const
		{
			return GLOSSY_GLASS_ID_COLOR;
		}

		virtual void sampleBSDF(Hit const& hit, DirectionSample& out, Math::Sampler& sampler, bool RADIANCE = true)const override
		{
			bool entering = hit.facing;
			double n1 = 1, n2 = m_eta;
			if (!entering)	std::swap(n1, n2);

			Math::Vector3f normal = hit.orientedPrimitiveNormal();

			double cos_theta = std::abs(hit.to_view * normal);
			double sin_theta = std::sqrt(1 - cos_theta * cos_theta);
			out.bsdf = m_albedo;
			if (sin_theta == 0)
			{
				out.direction = -normal;
				out.pdf = 1;
			}
			else
			{
				Math::Vector3f tg = (hit.to_view - normal * cos_theta) / (sin_theta);
				double eta_ratio = n1 / n2;
				double next_sin_theta = eta_ratio * sin_theta;

				bool reflect = next_sin_theta >= 1;
				double fresnel_reflectance = 1;
				double pdf_reflect = 1;

				double next_cos_theta;

				if (!reflect)
				{
					next_cos_theta = std::sqrt(1 - next_sin_theta * next_sin_theta);
					double rs = (n2 * cos_theta - n1 * next_cos_theta) / (n2 * cos_theta + n1 * next_cos_theta);
					double rp = (n1 * cos_theta - n2 * next_cos_theta) / (n1 * cos_theta + n2 * next_cos_theta);

					fresnel_reflectance = 0.5 * (rs * rs + rp * rp);
					pdf_reflect = fresnel_reflectance;

					double xi = sampler.generateContinuous<double>();

					reflect = xi <= pdf_reflect;
				}


				if (reflect)
				{
					out.direction = hit.primitive_reflected();
					out.bsdf *= fresnel_reflectance;
					out.pdf = pdf_reflect;
				}
				else // transmit
				{
					out.direction = normal * (-next_cos_theta) + tg * (-next_sin_theta);
					out.bsdf *= (1 - fresnel_reflectance);
					/*
					if (n2 > n1)
					{
						out.bsdf *= ((n2 * n2) / (n1 * n1));
					}
					*/
					out.pdf = 1 - pdf_reflect;
				}
			}

			if (RADIANCE)
				out.bsdf *= 1.0 / std::abs(normal * hit.to_view);
			else
				out.bsdf *= 1.0 / std::abs(out.direction * normal);

		}

		virtual RGBColor BSDF(Hit const& hit, Math::Vector3f const& wi, Math::Vector3f const& wo, bool RADIANCE = false)const override
		{
			return 0;
		}

		virtual double pdf(Hit const& hit, Math::Vector3f const& wi, Math::Vector3f const& wo)const override
		{
			return 0;
		}

	};
}