#pragma once

#include <Geometry/RGBColor.h>
#include <Geometry/Texture.h>
#include <Geometry/Sample.h>
#include <string>
#include <Geometry/Hit.h>
#include <Math/RandomDirection.h>
#include <Math/Sampler.h>

#define PURE_EMISSIVE_ID_COLOR {1, 1, 1}
#define LAMBERT_ID_COLOR {1, 0, 0}
#define SPECULAR_ID_COLOR {0, 0, 0.5}
#define DELTA_MIRROR_ID_COLOR {0, 0, 1}
#define PHONG_ID_COLOR {1, 0, 0.5};
#define GLASS_ID_COLOR {0, 1, 0};
#define GLOSSY_GLASS_ID_COLOR {0, 0.5, 0};

namespace Geometry
{


	class Material
	{

		//////////////////////////////////////////////
		//TODO
		// manage copy of material
		/////////////////////////////////////////////



	protected:


		/// <summary> The emissive color</summary>
		RGBColor m_emissiveColor;
		double m_light_spickyness = 0;

		std::string m_textureFile;

		Texture * m_texture;

		bool m_delta=false;
		bool m_spicky = false;

		//double m_albedo=0;

		

		bool clear_texture()
		{
			if (m_texture != nullptr)
			{
				delete m_texture;
				m_texture = nullptr;
				return true;
			}
			return false;
		}

		bool load_texture()
		{
			::std::cout << "Loading texture: " << m_textureFile << "..." << ::std::flush;
			m_texture = new Texture(m_textureFile);
			if (!m_texture->isValid())
			{
				clear_texture();
				::std::cerr << "discarded" << ::std::endl;
			}
			else
			{
				::std::cout << "OK" << ::std::endl;
			}
			return true;
		}

	
		

	public:



		virtual RGBColor ID_COLOR()const
		{
			return PURE_EMISSIVE_ID_COLOR;
		}


		Material(RGBColor const& em = 0, std::string const& path = "", double s=0.0) :
			m_emissiveColor(em),
			m_textureFile(path),
			m_texture(nullptr),
			m_spicky(false),
			m_delta(false),
			m_light_spickyness(s)
		{
			if (!path.empty())
				load_texture();
		}



		virtual ~Material()
		{
			clear_texture();
		}
		

		/////////////////////////
		//Samples a direction acording to the pdf / bsdf, and returns (in out) the sample direction, bsdf, pdf, wight
		// Always generates a sample, but it can be a dummy sample, like this default implementation
		// what matters is that the bsdf is 0, and so there is no need to go further with the sample
		//Do I need the whole hit?
		/////////////////////////
		virtual void sampleBSDF(Hit const& hit, DirectionSample& out, Math::Sampler & sampler, bool RADIANCE=false, double * pdf_rev=nullptr)const
		{
			out.direction = hit.primitive_normal;
			out.pdf = 1;
			out.bsdf = 0;
		}

		
		inline RGBColor BSDF(Hit const& hit, Math::Vector3f const& wi, bool RADIANCE=false)const
		{
			return BSDF(hit, wi, hit.to_view, RADIANCE);
		}

		virtual RGBColor BSDF(Hit const& hit, Math::Vector3f const& wi, Math::Vector3f const& wo, bool RADIANCE=false)const
		{
			return 0;
		}


		//returns the probability of sampling wi knowing that we are comming from wo
		virtual double pdf(Hit const& hit, Math::Vector3f const& wi, Math::Vector3f const& wo, bool RADIANCE=false)const
		{
			return 0; // or 0 or 1, actually this function should not be called, ever?
		}

		__forceinline double pdf(Hit const& hit, Math::Vector3f const& wi, bool RADIANCE=false)const //implicitly, wo = hit.to_view
		{
			return pdf(hit, wi, hit.to_view, RADIANCE);
		}

		/// <summary>
		/// Gets the texture file.
		/// </summary>
		/// <returns></returns>
		const ::std::string & getTextureFile() const
		{
			return m_textureFile;
		}

		/// <summary>
		/// Returns the texture
		/// </summary>
		/// <returns></returns>
		const Texture & getTexture() const
		{
			return *m_texture;
		}

		/// <summary>
		/// Tests if a texture is associated with this material.
		/// </summary>
		/// <returns> true if this materail has a texture, false otherwise</returns>
		bool hasTexture() const
		{
			return m_texture != nullptr && m_texture->isValid();
		}

		__forceinline RGBColor getTexturePixel(Math::Vector2f const& uv)const
		{
			return hasTexture() ? m_texture->pixel(uv) : 1;
		}


		void setTextureFile(const ::std::string & textureFile)
		{
			m_textureFile = textureFile;
			load_texture();
		}




		/// <summary>
		/// Sets the emissive color.
		/// </summary>
		/// <param name="color">The color.</param>
		void setEmissive(RGBColor const & color)
		{
			m_emissiveColor = color;
		}

		/// <summary>
		/// Gets the emissive color.
		/// </summary>
		/// <returns></returns>
		const RGBColor & getEmissive() const
		{
			return m_emissiveColor;
		}


		RGBColor Le(Math::Vector3f const& normal, Math::Vector2f const uv, Math::Vector3f const& dir)const
		{
			RGBColor res = getEmissive();
			if (hasTexture())
			{
				res = res * getTexture().pixel(uv);
			}

			double cosine = normal * dir;
			if (cosine < 0)	return 0;
			double scaling = (m_light_spickyness + 1) / Math::twoPi;
			res *= std::pow(cosine, m_light_spickyness) * scaling * 4;
			return res;
		}


		bool is_emissive()const
		{
			return !m_emissiveColor.isBlack();
		}


		bool spicky()const
		{
			return m_spicky;
		}

		bool delta()const
		{
			return m_delta;
		}


		virtual DirectionSample sampleLightDirection(SurfaceSample const& sls, Math::Sampler& sampler)const
		{
			DirectionSample res;
			Math::OldRandomDirection diffuse_sampler(&sampler, sls.normal, m_light_spickyness+1);
			res.direction = diffuse_sampler.generate();
			res.bsdf = Le(sls.normal, sls.uv, res.direction);
			res.pdf = std::pow(res.direction * sls.normal, m_light_spickyness+1.0) * ((m_light_spickyness+2.0) / Math::twoPi);
			assert(res.direction * sls.normal > 0);
			return res;
		}

		//assumes dir is normalized
		virtual double pdfLight(Hit const& hit, Math::Vector3f const& dir)const
		{
			return std::pow(std::max(dir * hit.primitive_normal, 0.0), m_light_spickyness + 1.0) * ((m_light_spickyness + 2.0) / Math::twoPi);;
		}



		//double albedo()const
		//{
		//	return m_albedo;
		//}
		

	};

}
