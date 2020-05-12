#pragma once

#include <Geometry/Materials/Material.h>
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

		virtual void sampleBSDF(Hit const& hit, DirectionSample& out, Math::Sampler& sampler, bool RADIANCE=true, double *pdf_rev=nullptr)const override
		{
			bool entering = hit.facing;
			double nr = 1, nt = m_eta;
			if (!entering)	std::swap(nr, nt);

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
				double eta_ratio = nr / nt;
				double next_sin_theta = eta_ratio * sin_theta;

				bool reflect = next_sin_theta >= 1;
				double fresnel_reflectance = 1;
				double pdf_reflect = 1;

				double next_cos_theta;
				
				if (!reflect)
				{
					next_cos_theta = std::sqrt(1 - next_sin_theta * next_sin_theta);
					double rs = (nt * cos_theta - nr * next_cos_theta) / (nt * cos_theta + nr * next_cos_theta);
					double rp = (nr * cos_theta - nt * next_cos_theta) / (nr * cos_theta + nt * next_cos_theta);

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
					if (RADIANCE)
						out.bsdf *= (eta_ratio * eta_ratio);
					out.pdf = 1 - pdf_reflect;
				}
			}
			
			out.bsdf *= 1.0 / std::abs(out.direction * normal);
			if (pdf_rev)
				*pdf_rev = out.pdf;
		}

		virtual RGBColor BSDF(Hit const& hit, Math::Vector3f const& wi, Math::Vector3f const& wo, bool RADIANCE = false)const override
		{
			return 0;
		}

		virtual double pdf(Hit const& hit, Math::Vector3f const& wi, Math::Vector3f const& wo, bool RADIANCE = false)const override
		{
			return 0;
		}

	};

	class DispertionGlass : public Material
	{
	protected:

		RGBColor m_albedo;

		double m_eta_min, m_eta_max;
		double wlmin = 380, wlmax = 380 + 290;

	public:

		DispertionGlass(RGBColor const& albedo, double eta1, double eta2, std::string const& texture = "") :
			Material(0, texture),
			m_albedo(albedo),
			m_eta_min(eta1),
			m_eta_max(eta2)
		{
			m_delta = true;
			m_spicky = true;
		}

		virtual RGBColor ID_COLOR()const
		{
			return { 0.2, 0.8, 0.2 };
		}

		virtual void sampleBSDF(Hit const& hit, DirectionSample& out, Math::Sampler& sampler, bool RADIANCE = true, double* pdf_rev = nullptr)const override
		{
			bool entering = hit.facing;
			double nr = 1;
			double xi = sampler.generateContinuous(0.0, 1.0);
			double nt = m_eta_min + (m_eta_max - m_eta_min) * xi;
			double wl = wlmin + (wlmax - wlmin) * xi;
			
			if (!entering)	std::swap(nr, nt);

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
				double eta_ratio = nr / nt;
				double next_sin_theta = eta_ratio * sin_theta;

				bool reflect = next_sin_theta >= 1;
				double fresnel_reflectance = 1;
				double pdf_reflect = 1;

				double next_cos_theta;

				if (!reflect)
				{
					next_cos_theta = std::sqrt(1 - next_sin_theta * next_sin_theta);
					double rs = (nt * cos_theta - nr * next_cos_theta) / (nt * cos_theta + nr * next_cos_theta);
					double rp = (nr * cos_theta - nt * next_cos_theta) / (nr * cos_theta + nt * next_cos_theta);

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
					if (RADIANCE)
						out.bsdf *= (eta_ratio * eta_ratio);
					out.pdf = 1 - pdf_reflect;
				}
			}

			out.bsdf *= 1.0 / std::abs(out.direction * normal);
			RGBColor wl_color = RGBColor::fromWaveLength(wl);
			out.bsdf *= wl_color * 2.0;
			if (pdf_rev)
				*pdf_rev = out.pdf;
		}

		virtual RGBColor BSDF(Hit const& hit, Math::Vector3f const& wi, Math::Vector3f const& wo, bool RADIANCE = false)const override
		{
			return 0;
		}

		virtual double pdf(Hit const& hit, Math::Vector3f const& wi, Math::Vector3f const& wo, bool RADIANCE = false)const override
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

		virtual void sampleBSDF(Hit const& hit, DirectionSample& out, Math::Sampler& sampler, bool RADIANCE = true, double *pdf_rev=nullptr)const override
		{
			double scaling = (m_shininess + 1) / Math::twoPi;
			Math::Vector3f normal = hit.orientedPrimitiveNormal();
			bool entering = hit.facing;
			double nr = 1, nt = m_eta;
			if (!entering)	std::swap(nr, nt);
			
			double cos_theta = std::abs(hit.to_view * normal);
			double sin_theta = std::sqrt(1 - cos_theta * cos_theta);

			bool reflect = sin_theta != 0;
			
			if (!reflect)
			{
				Math::RandomDirection direction_sampler(&sampler, -normal, m_shininess);
				out.direction = direction_sampler.generate();
				double cos_t = out.direction * -normal;
				double f = scaling * std::pow(cos_t, m_shininess);
				out.pdf = f;
				out.bsdf = m_albedo * f;
				if (!RADIANCE)
					out.bsdf *= 1.0 / std::abs(normal * hit.to_view);
				return;
			}

			
			double eta_ratio = nr / nt;
			double next_sin_theta = eta_ratio * sin_theta;

			bool transmit = next_sin_theta < 1;

			const Math::Vector3f reflected = hit.primitive_reflected();
			

			if (transmit) // and reflects
			{
				const Math::Vector3f tg = (hit.to_view - normal * cos_theta) / (sin_theta);
				double next_cos_theta = std::sqrt(1 - next_sin_theta * next_sin_theta);
				const Math::Vector3f transmited = normal * (-next_cos_theta) + tg * (-next_sin_theta);

				double rs = (nt * cos_theta - nr * next_cos_theta) / (nt * cos_theta + nr * next_cos_theta);
				double rp = (nr * cos_theta - nt * next_cos_theta) / (nr * cos_theta + nt * next_cos_theta);
				double fresnel_reflectance = 0.5 * (rs * rs + rp * rp);
				double pdf_reflect = fresnel_reflectance;
				double xi = sampler.generateContinuous<double>();

				const Math::Vector3f & main_direction = xi < fresnel_reflectance ? reflected : transmited;

				Math::RandomDirection direction_sampler(&sampler, main_direction, m_shininess);
				out.direction = direction_sampler.generate();

				bool has_reflected = out.direction * normal > 0;
				double cfr = std::max(reflected * out.direction, 0.0);
				double cft = std::max(transmited * out.direction, 0.0);

				double pfr = std::pow(cfr, m_shininess);
				double pft = std::pow(cft, m_shininess);

				out.pdf = (pfr * pdf_reflect + pft * (1 - pdf_reflect)) * scaling;
				out.bsdf = m_albedo * (has_reflected ? pfr * fresnel_reflectance : pft * (1 - fresnel_reflectance)) * scaling;
				if (RADIANCE && !has_reflected)
					out.bsdf *= (eta_ratio * eta_ratio);
			}
			else // only reflects
			{
				Math::RandomDirection direction_sampler(&sampler, reflected, m_shininess);
				out.direction = direction_sampler.generate();
				double cos_reflected = reflected * out.direction;
				double f = scaling * std::pow(cos_reflected, m_shininess);
				out.pdf = f;
				out.bsdf = m_albedo * f;
			}
			
			out.bsdf *= 1.0 / std::abs(out.direction * normal);
			if (pdf_rev)
				*pdf_rev = out.pdf;
		}

		virtual RGBColor BSDF(Hit const& hit, Math::Vector3f const& wi, Math::Vector3f const& wo, bool RADIANCE = false)const override
		{
			double scaling = (m_shininess + 1) / Math::twoPi;
			Math::Vector3f normal = hit.normal; if (wo * normal < 0) normal = -normal;
			bool entering = hit.facing;
			double nr = 1, nt = m_eta;
			if (!entering)	std::swap(nr, nt);

			double cos_theta = std::abs(wo * normal);
			double sin_theta = std::sqrt(1 - cos_theta * cos_theta);

			bool reflect = sin_theta != 0;

			if (!reflect)
			{
				double cos_t = wi * -normal;
				double f = scaling * std::pow(cos_t, m_shininess);
				RGBColor res = m_albedo * f;
				if (RADIANCE)
					res *= 1.0 / std::abs(normal * wo);
				return res;
			}

			double eta_ratio = nr / nt;
			double next_sin_theta = eta_ratio * sin_theta;

			bool transmit = next_sin_theta < 1;

			const Math::Vector3f reflected = hit.primitive_normal.reflect(wo);

			RGBColor res;
			if (transmit) // and reflects
			{
				const Math::Vector3f tg = (wo - normal * cos_theta) / (sin_theta);
				double next_cos_theta = std::sqrt(1 - next_sin_theta * next_sin_theta);
				const Math::Vector3f transmited = normal * (-next_cos_theta) + tg * (-next_sin_theta);

				double rs = (nt * cos_theta - nr * next_cos_theta) / (nt * cos_theta + nr * next_cos_theta);
				double rp = (nr * cos_theta - nt * next_cos_theta) / (nr * cos_theta + nt * next_cos_theta);
				double fresnel_reflectance = 0.5 * (rs * rs + rp * rp);

				bool has_reflected = wi * normal > 0;
				double cfr = std::max(reflected * wi, 0.0);
				double cft = std::max(transmited * wi, 0.0);

				double pfr = std::pow(cfr, m_shininess);
				double pft = std::pow(cft, m_shininess);

				res = m_albedo * (has_reflected ? pfr * fresnel_reflectance : pft * (1 - fresnel_reflectance)) * scaling;
				if (RADIANCE && !has_reflected)
					res *= (eta_ratio * eta_ratio);
			}
			else // only reflects
			{
				double cos_reflected = reflected * wi;
				double f = scaling * std::pow(cos_reflected, m_shininess);
				res = m_albedo * f;
			}

			res *= 1.0 / std::abs(wi * normal);
			return res;
		}

		virtual double pdf(Hit const& hit, Math::Vector3f const& wi, Math::Vector3f const& wo, bool RADIANCE=false)const override
		{
			double scaling = (m_shininess + 1) / Math::twoPi;
			Math::Vector3f normal = hit.normal; if (wo * normal < 0) normal = -normal;
			bool entering = hit.facing;
			double nr = 1, nt = m_eta;
			if (!entering)	std::swap(nr, nt);

			double cos_theta = std::abs(wo * normal);
			double sin_theta = std::sqrt(1 - cos_theta * cos_theta);

			bool reflect = sin_theta != 0;

			if (!reflect)
			{
				double cos_t = wi * -normal;
				double f = scaling * std::pow(cos_t, m_shininess);
				return f;
			}

			double eta_ratio = nr / nt;
			double next_sin_theta = eta_ratio * sin_theta;

			bool transmit = next_sin_theta < 1;

			const Math::Vector3f reflected = hit.primitive_normal.reflect(wo);

			RGBColor res;
			if (transmit) // and reflects
			{
				const Math::Vector3f tg = (wo - normal * cos_theta) / (sin_theta);
				double next_cos_theta = std::sqrt(1 - next_sin_theta * next_sin_theta);
				const Math::Vector3f transmited = normal * (-next_cos_theta) + tg * (-next_sin_theta);

				double rs = (nt * cos_theta - nr * next_cos_theta) / (nt * cos_theta + nr * next_cos_theta);
				double rp = (nr * cos_theta - nt * next_cos_theta) / (nr * cos_theta + nt * next_cos_theta);
				double fresnel_reflectance = 0.5 * (rs * rs + rp * rp);
				double pdf_reflect = fresnel_reflectance;
				bool has_reflected = wi * normal > 0;
				double cfr = std::max(reflected * wi, 0.0);
				double cft = std::max(transmited * wi, 0.0);

				double pfr = std::pow(cfr, m_shininess);
				double pft = std::pow(cft, m_shininess);

				return (pfr * pdf_reflect + pft * (1 - pdf_reflect)) * scaling;
			}
			else // only reflects
			{
				double cos_reflected = reflected * wi;
				double f = scaling * std::pow(cos_reflected, m_shininess);
				return f;
			}
		}

	};
}