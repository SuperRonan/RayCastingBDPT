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
			Material(emissive, false, true, texture),
			m_specular(specular /* Math::pi*/)
		{}



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


		virtual void sampleBSDF(Hit const& hit, unsigned int diffuse_samples, unsigned int specular_samples, DirectionSample& out, Math::Sampler& sampler)const override
		{
			out.direction = hit.primitive_reflected();
			double cosd = hit.normal * out.direction;
			out.pdf = 1;
			out.bsdf = m_specular;
		}

		virtual RGBColor BSDF(Hit const& hit, Math::Vector3f const& wi, Math::Vector3f const& wo)const override
		{
			return _BSDF(hit, wi, wo);
		}


		virtual void sampleBSDF(Hit const& hit, unsigned int diffuse_samples, unsigned int specular_samples, DirectionStack& res, Math::Sampler& sampler)const override
		{
			DeltaMirror::sampleBSDF(hit, diffuse_samples, specular_samples, *res.end(), sampler);
			res.grow();
		}


		virtual void BSDF(Hit const& hit, ColorStack& res, LightSampleStack const& wis)const override
		{
			for (SurfaceLightSample const& wi : wis)
			{
				res.push(_BSDF(hit, wi.vector, hit.to_view));
			}
		}

		virtual double pdf(Hit const& hit, Math::Vector3f const& wi, Math::Vector3f const& wo)const override
		{
			//or maybe always 0?
			return hit.primitive_normal.reflect(wo) == wi ? 1.0 : 0.0;
		}
	};
}