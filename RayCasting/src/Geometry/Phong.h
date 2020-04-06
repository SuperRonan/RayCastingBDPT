#pragma once

#include "Material.h"
#include <Geometry/Lambert.h>
#include <Geometry/Specular.h>

namespace Geometry
{

	class Phong : public Material
	{
	protected:

	public:

		Phong(RGBColor const& dif = 1, RGBColor const& spec = 0, double s = 1, RGBColor const& em = 0, std::string const& tex = "") 
		{

		}


		virtual RGBColor ID_COLOR()const override
		{
			return PHONG_ID_COLOR;
		}

		virtual void sampleBSDF(Hit const& hit, unsigned int diffuse_samples, unsigned int specular_samples, DirectionSample& out, Math::Sampler& sampler, bool RADIANCE=true)const override
		{
			Material::sampleBSDF(hit, diffuse_samples, specular_samples, out, sampler);
			return;
		}


		virtual RGBColor BSDF(Hit const& hit, Math::Vector3f const& wi, Math::Vector3f const& wo)const override
		{
			return 0;
		}

		virtual double pdf(const Hit & hit, Math::Vector3f const& wi, Math::Vector3f const& wo)const  override
		{
			return Material::pdf(hit, wi, wo);
		}

	};

	/*
	class PhongMaterial : public Material
	{
	protected:
		/// <summary> The ambient color </summary>
		RGBColor m_ambientColor;
		/// <summary> the diffuse color</summary>
		RGBColor m_diffuseColor;
		/// <summary> The specular color</summary>
		RGBColor m_specularColor;
		/// <summary> The shininess</summary>
		double    m_shininess;

	public:
		const medium * m_medium = nullptr;
		bool m_transparent = false;
		virtual RGBColor ID_COLOR()const
		{
			return PHONG_ID_COLOR;
		}
		PhongMaterial(RGBColor const & ambientColor = RGBColor(), RGBColor const & diffuseColor = RGBColor(),
			RGBColor specularColor = RGBColor(), double shininess = 1.0, RGBColor const & emissiveColor = RGBColor(),
			std::string const & textureFile = "") :
			Material(emissiveColor, textureFile),
			m_ambientColor(ambientColor),
			m_diffuseColor(diffuseColor / Math::pi),
			m_specularColor(specularColor),
			m_shininess(shininess)
		{}
		/// <summary>
		/// Sets the ambient color.
		/// </summary>
		/// <param name="color">The color.</param>
		void setAmbient(RGBColor const & color)
		{
			m_ambientColor = color;
		}
		/// <summary>
		/// Gets the ambient color.
		/// </summary>
		/// <returns>The ambiant color</returns>
		const RGBColor & getAmbient() const
		{
			return m_ambientColor;
		}
		/// <summary>
		/// Sets the diffuse color.
		/// </summary>
		/// <param name="color">The color.</param>
		void setDiffuse(RGBColor const & color)
		{
			m_diffuseColor = color / Math::pi;
		}
		/// <summary>
		/// Gets the diffuse color.
		/// </summary>
		/// <returns></returns>
		const RGBColor & getDiffuse() const
		{
			return m_diffuseColor;
		}
		/// <summary>
		/// Sets the specular color.
		/// </summary>
		/// <param name="color">The color.</param>
		void setSpecular(RGBColor const & color)
		{
			m_specularColor = color;
		}
		/// <summary>
		/// Gets the specular color.
		/// </summary>
		/// <returns></returns>
		const RGBColor & getSpecular() const
		{
			return m_specularColor;
		}
		/// <summary>
		/// Sets the shininess.
		/// </summary>
		/// <param name="s">The shininess.</param>
		void setShininess(double s)
		{
			m_shininess = s;
		}
		/// <summary>
		/// Gets the shininess.
		/// </summary>
		/// <returns></returns>
		const double & getShininess() const
		{
			return m_shininess;
		}

		//virtual void shader(Ray const& ray, Hit<double> const& hit, LightStack const& lights, ShaderOut & out)const override
		//{
		//	RGBColor const& emissive = getEmissive();
		//	RGBColor const& ambient = 0;// getAmbient();
		//	RGBColor diffuse = 0;
		//	RGBColor specular = 0;
		//	const RGBColor tex_color = hasTexture() ? getTexture().safe_pixel(hit.tex_uv) : 1;
		//	const Math::Vector3f & reflected = hit.reflected;
		//	for (size_t i = 0; i < lights.size(); ++i)
		//	{
		//		const DirectionalLight & l = lights[i];
		//		const Math::Vector3f & to_light = l.direction();
		//		diffuse = diffuse + getDiffuse() * l.color() * (hit.normal * to_light);
		//		if (reflected * to_light > 0)
		//		{
		//			specular = specular + getSpecular() * l.color() * pow(to_light * reflected, getShininess());
		//		}
		//	}
		//	out.rest = (ambient + diffuse + specular) * tex_color;
		//	out.emisive = emissive * tex_color;
		//}
		bool is_transparent()const
		{
			return m_transparent;
		}
		/*virtual void apply_micro_facets(Hit<double> & hit)const override
		{

		}
		//virtual void compute_out_dir(DirectionStack& out, Hit<double> const& hit, unsigned int diffuse_samples, unsigned int specular_samples)const override
		//{
		//	//Multi importance sampling ??
		//	//
		//	//TODO chack to proba, and the number of samples, I think it is wierd
		//	//
		//	if (!m_diffuseColor.isBlack())
		//	{
		//		Math::RandomDirection diffuse_sampler(hit.primitive_normal);
		//		for (int i = 0; i < diffuse_samples; ++i)
		//		{
		//			Math::Vector3f direction = diffuse_sampler.generate();
		//			if (direction * hit.primitive_normal < 0)
		//			{
		//				direction = -direction;
		//			}
		//			//TODO unbias
		//			double probability = direction * hit.primitive_normal / Math::pi;
		//			//double bias = 1;
		//			out.push({ probability, direction });
		//		}
		//	}
		//	if (!m_specularColor.isBlack())
		//	{
		//		Math::RandomDirection specular_sampler(hit.reflected, m_shininess);
		//		for (int i = 0; i < specular_samples; ++i)
		//		{
		//			Math::Vector3f direction = specular_sampler.generate();
		//			if (direction * hit.normal < 0)
		//			{
		//				continue;
		//			}
		//			//wrong
		//			double probability = 1.0 / Math::twoPi;
		//			out.push({ probability, direction });
		//		}
		//	}
		//}
		virtual bool sampleBSDF(Hit<double> const& hit, unsigned int diffuse_samples, unsigned int specular_samples, DirectionSample& out, Math::Sampler& sampler)const override
		{
			//for now, it only woks with only diffuse materials, no blending of the two kinds of bsdf yet
			RGBColor tex = hasTexture() ? getTexture().pixel(hit.tex_uv) : 1;
			RGBColor dif = tex * m_diffuseColor;
			RGBColor const& spec = m_specularColor;

			if (!dif.isBlack())
			{
				Math::RandomDirection diffuse_sampler(&sampler, hit.primitive_normal, 1);
				double pdf = 0;
				Math::Vector3f dir;
				while (pdf == 0)
				{
					dir = diffuse_sampler.generate();
					pdf = (dir * hit.primitive_normal) / Math::pi;
				}
				if (pdf < 0)
				{
					dir = -dir;
					pdf = -pdf;
				}
				assert(pdf > 0 && pdf <= 1);
				//TODO MIS
				double weight = 1.0 ;
				RGBColor bsdf = dif + spec * pow(hit.reflected * dir, m_shininess);
				out = { weight, pdf, bsdf, dir };
				return true;
			}
		}
		virtual RGBColor BSDF(Hit<double> const& hit, Math::Vector3f const& wi)const override
		{
			RGBColor tex = hasTexture() ? getTexture().pixel(hit.tex_uv) : 1;
			RGBColor dif = tex * m_diffuseColor;
			RGBColor const& spec = m_specularColor;
			//TODO manage transparancy
			return wi * hit.primitive_normal > 0 ? dif + spec * pow(hit.reflected * wi, m_shininess) : 0;
		}
		virtual void BSDF(Hit<double> const& hit, ColorStack& res, LightSampleStack const& wis)const override
		{
			RGBColor tex = hasTexture() ? getTexture().pixel(hit.tex_uv) : 1;
			RGBColor dif = tex * m_diffuseColor;
			RGBColor const& spec = m_specularColor;
			for (SurfaceLightSample wi : wis)
			{
				Math::Vector3f const& direction = wi.vector;
				//TODO transparancy
				RGBColor bsdf = direction * hit.primitive_normal > 0 ? dif + spec * pow(hit.reflected * direction, m_shininess) : 0;
				res.push(bsdf);
			}
		}


		virtual void sampleBSDF(Hit<double> const& hit, unsigned int diffuse_samples, unsigned int specular_samples, DirectionStack& out, Math::Sampler& sampler)const override
		{
			RGBColor tex = hasTexture() ? getTexture().pixel(hit.tex_uv) : 1;
			RGBColor dif = tex * m_diffuseColor;
			RGBColor const& spec = m_specularColor;
			if (dif.isBlack())
			{
				diffuse_samples = 0;
			}
			if (spec.isBlack())
			{
				specular_samples = 0;
			}
			unsigned int total_samples = diffuse_samples + specular_samples;
			//sampling directions according to a cos pdf
			if (diffuse_samples>0)
			{
				Math::RandomDirection diffuse_sampler(&sampler, hit.primitive_normal, 1);
				for (unsigned int i = 0; i < diffuse_samples; ++i)
				{
					double pdf=0;
					Math::Vector3f dir;
					while (pdf == 0)
					{
						dir = diffuse_sampler.generate();
						pdf = (dir * hit.primitive_normal) / Math::pi;
					}
					if (pdf < 0)
					{
						dir = -dir;
						pdf = -pdf;
					}
					assert(pdf > 0 && pdf <= 1);
					//TODO MIS
					double weight = 1.0 / total_samples;
					RGBColor bsdf = dif + spec * pow(hit.reflected * dir, m_shininess);
					out.push({ weight, pdf, bsdf, dir });
				}
			}
			//sampling directions according to a specular pdf
			if (specular_samples>0)
			{
				Math::RandomDirection spec_sampler(&sampler, hit.reflected, m_shininess);
				for (unsigned int i = 0; i < specular_samples; ++i)
				{
					Math::Vector3f dir = spec_sampler.generate();
					if (dir * hit.normal < 0)
					{
						continue;
					}
					//TODO check if it really works all the time
					double bsdfd = pow(hit.reflected * dir, m_shininess);
					double pdf = (m_shininess + 1) * bsdfd / (Math::twoPi);
					assert(pdf > 0 && pdf <= 1);
					//TODO MIS
					double weight = 1.0 / total_samples;
					RGBColor bsdf = dif + spec * bsdfd;
					out.push({ weight, pdf, bsdf, dir });
				}
			}
			//TODO add MIS
		}
		/*
		virtual Material * new_copy()const
		{
			return new PhongMaterial(*this);
		}

	};
	*/

}