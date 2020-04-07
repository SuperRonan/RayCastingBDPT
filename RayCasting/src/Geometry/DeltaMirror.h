#pragma once


#include <Geometry/Material.h>
#include <Math/Constant.h>


namespace Geometry
{

	class DeltaMirror :public Material
	{
	protected:

		RGBColor m_specular;


		

	public:

		virtual RGBColor ID_COLOR()const
		{
			return DELTA_MIRROR_ID_COLOR;
		}



		DeltaMirror(RGBColor const& specular, RGBColor const& emissive = 0, std::string const& texture = "") :
			Material(emissive, texture),
			m_specular(specular /* Math::pi*/)
		{
			m_delta = true;
			m_spicky = true;
		}



		void setSpecular(RGBColor const& color)
		{
			m_specular = color /* Math::pi*/;
		}

		const RGBColor& getSpecular() const
		{
			return m_specular;
		}

		

	protected:


		__forceinline RGBColor _BSDF(Hit const& hit, Math::Vector3f const& wi, Math::Vector3f const& wo)const
		{
			return 0;
			return m_specular * (hit.primitive_normal.reflect(wo) == wi);
		}

	public:


		virtual void sampleBSDF(Hit const& hit, DirectionSample& out, Math::Sampler& sampler, bool RADIANCE=true)const override
		{
			out.direction = hit.primitive_reflected();
			double cosd = hit.normal * out.direction;
			out.pdf = 1;
			out.bsdf = m_specular;

			if (RADIANCE)
				out.bsdf *= 1.0 / std::abs(hit.to_view * hit.normal);
			else
				out.bsdf *= 1.0 / std::abs(out.direction * hit.normal);
		}

		virtual RGBColor BSDF(Hit const& hit, Math::Vector3f const& wi, Math::Vector3f const& wo, bool RADIANCE = false)const override
		{
			return _BSDF(hit, wi, wo);
		}

		virtual double pdf(Hit const& hit, Math::Vector3f const& wi, Math::Vector3f const& wo)const override
		{
			//or maybe always 0?
			return hit.primitive_normal.reflect(wo) == wi ? 1.0 : 0.0;
		}
	};
}