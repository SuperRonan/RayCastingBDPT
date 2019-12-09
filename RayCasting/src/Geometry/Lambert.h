#pragma once

#include <Geometry/Material.h>
#include <Math/Constant.h>


namespace Geometry
{

	class Lambertian :public Material
	{
	protected:

		RGBColor m_diffuse;



	public:

		Lambertian(RGBColor const& diffuse, RGBColor const& emissive=0, std::string const& texture="") :
			Material(emissive, true, false, texture),
			m_diffuse(diffuse / Math::pi)
		{
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


		virtual void sampleBSDF(Hit const& hit, unsigned int diffuse_samples, unsigned int specular_samples, DirectionSample& out, Math::Sampler& sampler, bool _=true)const override
		{
			RGBColor dif = m_diffuse * getTexturePixel(hit.tex_uv);
			Math::Vector3f normal = hit.orientedPrimitiveNormal();
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

		virtual RGBColor BSDF(Hit const& hit, Math::Vector3f const& wi, Math::Vector3f const& wo)const override
		{
			if (((wi * hit.primitive_normal) * (wo * hit.primitive_normal)) > 0)
				return m_diffuse * getTexturePixel(hit.tex_uv);
			return 0;
		}
		

		virtual void sampleBSDF(Hit const& hit, unsigned int diffuse_samples, unsigned int specular_samples, DirectionStack& res, Math::Sampler& sampler)const override
		{
			RGBColor dif = m_diffuse * getTexturePixel(hit.tex_uv);
			if (dif.isBlack())	return;
			Math::RandomDirection direction_sampler(&sampler, hit.orientedPrimitiveNormal());
			for (unsigned int i = 0; i < diffuse_samples; ++i)
			{
				DirectionSample out;
				out.direction = direction_sampler.generate();
				out.pdf = (out.direction * hit.orientedPrimitiveNormal()) / Math::pi;
				if (out.pdf == 0)	continue; //let's say the sample is not generated
				if (out.pdf < 0)
				{
					out.direction = -out.direction;
					out.pdf = -out.pdf;
				}
				out.bsdf = dif; 
				res.push(out);
			}
		}


		virtual void BSDF(Hit const& hit, ColorStack& res, LightSampleStack const& wis)const override
		{
			for (SurfaceLightSample const& wi : wis)
			{
				//should be inlined (I hope)
				res.push(BSDF(hit, wi.vector, hit.to_view));
				//res.push(m_diffuse * getTexturePixel(hit.tex_uv) * (wi.vector * hit.primitive_normal > 0));
			}
		}


		virtual double pdf(Hit const& hit, Math::Vector3f const& wi, Math::Vector3f const& wo)const override
		{
			if ((wi * hit.primitive_normal) * (wo * hit.primitive_normal) > 0)
			{
				double res = std::abs(hit.primitive_normal * wi) / Math::pi;
				return res;
			}
			return 0;
		}
	};
}