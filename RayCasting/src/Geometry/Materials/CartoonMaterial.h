#pragma once


namespace Geometry
{
	//class CartoonMaterial : public Phong
	//{

	//public:


	//	CartoonMaterial(RGBColor const & ambientColor = RGBColor(), RGBColor const & diffuseColor = RGBColor(),
	//		RGBColor specularColor = RGBColor(), double shininess = 1.0, RGBColor const & emissiveColor = RGBColor(),
	//		std::string const & textureFile = ""):
	//		Phong(ambientColor, diffuseColor, specularColor, shininess, emissiveColor, textureFile)
	//	{}

	//	virtual void shader(Ray const& ray, Hit<double> const& hit, LightStack const& lights, ShaderOut & out)const
	//	{
	//		RGBColor const& emissive = getEmissive();
	//		RGBColor const& ambient = getAmbient();
	//		RGBColor diffuse = 0;
	//		RGBColor specular = 0;

	//		const RGBColor tex_color = hasTexture() ? getTexture().safe_pixel(hit.tex_uv) : 1;

	//		const Math::Vector3f & reflected = hit.reflected;
	//		for (size_t i = 0; i < lights.size(); ++i)
	//		{
	//			const DirectionalLight & l = lights[i];
	//			const Math::Vector3f & to_light = l.direction();


	//			double diffuse_weight = hit.normal * to_light;
	//			if (diffuse_weight < 0.2)
	//			{
	//				diffuse_weight = 0;
	//			}
	//			else if (diffuse_weight < 0.7)
	//			{
	//				diffuse_weight = 0.5;
	//			}
	//			else
	//			{
	//				diffuse_weight = 0.9;
	//			}

	//			diffuse = diffuse + getDiffuse() * l.color() * diffuse_weight;

	//			if (reflected * to_light > 0)
	//			{
	//				double specular_weight = pow(to_light * reflected, getShininess());

	//				if (specular_weight < 0.2)
	//				{
	//					specular_weight = 0;
	//				}
	//				else if (specular_weight < 0.7)
	//				{
	//					specular_weight = 0.5;
	//				}
	//				else
	//				{
	//					specular_weight = 0.9;
	//				}

	//				specular = specular + getSpecular() * l.color() * specular_weight;
	//			}
	//		}
	//		out.rest = (ambient + diffuse + specular) * tex_color;
	//		out.emisive = emissive * tex_color;
	//	}



	//	/*
	//	virtual Material * new_copy()const
	//	{
	//		return new CartoonMaterial(*this);
	//	}
	//	*/
	//};
}