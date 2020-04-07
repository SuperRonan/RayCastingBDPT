#pragma once

#include <Geometry/Material.h>
#include <Math/Constant.h>


namespace Geometry
{
	enum LAMBERT_MODE {REFLECT, TRANSMIT};
	
	template <LAMBERT_MODE MODE>
	class Lambertian :public Material
	{
	protected:

		RGBColor m_diffuse;

	public:

		Lambertian(RGBColor const& diffuse, RGBColor const& emissive=0, std::string const& texture="") :
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


		virtual void sampleBSDF(Hit const& hit, DirectionSample& out, Math::Sampler& sampler, bool RADIANCE=false)const override
		{
			RGBColor dif = m_diffuse * getTexturePixel(hit.tex_uv);
			Math::Vector3f normal = hit.orientedPrimitiveNormal();
			if constexpr (MODE == TRANSMIT)
				normal = -normal;
			//There is not really a point of sampling a 
			Math::RandomDirection direction_sampler(&sampler, normal);
			out.direction = direction_sampler.generate();
			out.pdf = (out.direction * normal) / Math::pi;
			if (out.pdf == 0)
			{
				//this sould never happen, but it eventually does.
				//I don't really know what to do...
				out.bsdf = 0;
				return;
			}
			if (out.pdf < 0)
			{
				out.direction = -out.direction;
				out.pdf = -out.pdf;
			}
			
			out.bsdf = dif;
		}

		virtual RGBColor BSDF(Hit const& hit, Math::Vector3f const& wi, Math::Vector3f const& wo, bool RADIANCE=false)const override
		{
			if constexpr (MODE == REFLECT)
			{
				if (((wi * hit.primitive_normal) * (wo * hit.primitive_normal)) > 0)
					return m_diffuse * getTexturePixel(hit.tex_uv);
			}
			else
			{
				if (((wi * hit.primitive_normal) * (wo * hit.primitive_normal)) < 0)
					return m_diffuse * getTexturePixel(hit.tex_uv);
			}
			return 0;
		}



		virtual double pdf(Hit const& hit, Math::Vector3f const& wi, Math::Vector3f const& wo)const override
		{
			if constexpr (MODE == REFLECT)
			{
				if ((wi * hit.primitive_normal) * (wo * hit.primitive_normal) > 0)
				{
					double res = std::abs(hit.primitive_normal * wi) / Math::pi;
					return res;
				}
			}
			else
			{
				if ((wi * hit.primitive_normal) * (wo * hit.primitive_normal) < 0)
				{
					double res = std::abs(hit.primitive_normal * wi) / Math::pi;
					return res;
				}
			}
			return 0;
		}
	};
}