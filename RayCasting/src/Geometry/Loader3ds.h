#pragma once

#include <lib3ds/file.h>
#include <lib3ds/mesh.h>
#include <lib3ds/material.h>
#include <string>
#include <iostream>
#include <vector>
#include <algorithm>
#include <iterator>
#include <map>
#include <Geometry/Materials/Material.h>
#include <Geometry/Shapes/Geometry.h>
#include <Geometry/Materials/Lambert.h>
#include <Geometry/Materials/Glossy.h>
#include <Geometry/Materials/Phong.h>
#include <settings.h>


namespace Geometry
{
	class Loader3ds
	{
	private:
		/// \brief The 3DS file name
		Lib3dsFile * m_file;
		/// \brief The loaded materials
		::std::map<::std::string, Material *> m_materials;
		/// \brief The loaded meshes
		::std::vector<GeometryCollection*> m_meshes ;




		static Material* bestType(const Lib3dsMaterial* material)
		{
			RGBColor dif;
			dif.set(material->diffuse);
#ifdef FORCE_LAMBERT
			return new Lambertian<Geometry::REFLECT>(dif);
#endif
			RGBColor spec;
			spec.set(material->specular);
			if (spec.isBlack())
			{
				if (dif.isBlack())
				{
					//no emissive?
					return new Material;
				}
				else
				{
					return new Lambertian<Geometry::REFLECT>(dif);
				}
			}
			else
			{
				if (dif.isBlack())
				{
					//no emissive?
					return new Glossy(spec, material->shininess);
				}
				else
				{
					RGBColor ambient;
					ambient.set(material->ambient);
					return new Phong(dif, spec, material->shininess);
				}
			}
		}

		/// <summary>
		/// Loads a material.
		/// </summary>
		/// <param name="material">The material.</param>
		/// <param name="texturePath">The texture path.</param>
		void loadMaterial(Lib3dsMaterial* material, const ::std::string& texturePath)
		{
			Material* res = bestType(material);


			::std::string textureName = material->texture1_map.name;
			if (textureName != "")
			{
				res->setTextureFile(texturePath + "\\" + textureName);
			}
			::std::string materialName = material->name;
			m_materials[materialName] = res;
			::std::cout << "Loader3ds: created material " << materialName << ", texture: " << textureName << ::std::endl;
		}

		/// <summary>
		/// Loads the materials.
		/// </summary>
		/// <param name="texturePath">The texture path.</param>
		void loadMaterials(const ::std::string& texturePath)
		{
			for (Lib3dsMaterial* material = m_file->materials; material != NULL; material = material->next)
			{
				loadMaterial(material, texturePath);
			}
		}

		/// <summary>
		/// Loads a mesh.
		/// </summary>
		/// <param name="mesh">The mesh.</param>
		void loadMesh(Lib3dsMesh* mesh)
		{
			//::std::string meshName = mesh->name;
			::std::cout << "Loading mesh " << mesh->name << ::std::flush;
			GeometryCollection* geometry = new GeometryCollection(m_materials[mesh->faceL[0].material]);
			// 1.1 - We register the vertices
			for (unsigned int cpt = 0; cpt < mesh->points; ++cpt)
			{
				Lib3dsPoint pt = mesh->pointL[cpt];
				Math::Vector3f v = Math::makeVector(pt.pos[0], pt.pos[1], pt.pos[2]);
				geometry->addVertex(v);
				//vertices.push_back(v);
			}
			// 1.2 - We register the texture coordinates (if any)
			if (mesh->texels > 0)
			{
				::std::cout << ", has texture coordinates" << ::std::flush;
			}
			for (unsigned int cpt = 0; cpt < mesh->texels; ++cpt)
			{
				Math::Vector2f t = Math::makeVector(mesh->texelL[cpt][0], mesh->texelL[cpt][1]);
				geometry->addTextureCoordinates(t);
			}
			// 1.3 we compute the normals
			Lib3dsVector* normals = new Lib3dsVector[3 * mesh->faces];
			Math::Vector3f* vertexNormals = new Math::Vector3f[3 * mesh->faces];
			lib3ds_mesh_calculate_normals(mesh, normals);
			for (unsigned int cpt = 0; cpt < 3 * mesh->faces; ++cpt)
			{
				vertexNormals[cpt] = Math::makeVector(normals[cpt][0], normals[cpt][1], normals[cpt][2]).normalized();
			}
			// 2 - We register the triangles
			for (unsigned int cpt = 0; cpt < mesh->faces; ++cpt)
			{
				Lib3dsFace tmp = mesh->faceL[cpt];
				if (m_materials.find(tmp.material) == m_materials.end())
				{
					::std::cout << "Problem, material " << tmp.material << " not found..." << ::std::endl;
				}
				else
				{
					geometry->addTriangle(tmp.points[0], tmp.points[1], tmp.points[2], vertexNormals + (cpt * 3));
				}
			}
			// 3 - We register the mesh
			m_meshes.push_back(geometry);
			// 4 - Cleanup
			delete[] normals;
			::std::cout << std::endl;
		}

		/// <summary>
		/// Loads the meshes.
		/// </summary>
		void loadMeshes()
		{
			for(Lib3dsMesh * mesh = m_file->meshes ; mesh!=NULL ; mesh = mesh->next)
			{
				loadMesh(mesh) ;
			}
		}

	public:
		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	Loader3ds::Loader3ds(const ::std::string & filename, const ::std::string & texturePath);
		///
		/// \brief	Constructor which loads the provided 3ds file.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	08/02/2016
		///
		/// \param	filename   	Filename of the file.
		/// \param	texturePath	Full pathname of the texture file.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		Loader3ds(const ::std::string& filename, const ::std::string& texturePath)
		{
			::std::cout << "Loader3ds: loading file " << filename << ::std::endl;
			m_file = lib3ds_file_load(filename.c_str());
			if (m_file == NULL)
			{
				::std::cerr << "Loader3ds: unable to load file " + filename << ::std::endl;
				return;
			}
			loadMaterials(texturePath);
			loadMeshes();
			::std::cout << "Loader3ds: OK" << ::std::endl;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	const ::std::vector<GeometryCollection*> & Loader3ds::getMeshes()
		///
		/// \brief	Gets the loaded meshes.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	29/11/2015
		///
		/// \return	The meshes.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		const ::std::vector<GeometryCollection*>& getMeshes() const
		{
			return m_meshes;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	::std::vector<Material*> Loader3ds::getMaterials()
		///
		/// \brief	Gets the loaded materials.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	29/11/2015
		////////////////////////////////////////////////////////////////////////////////////////////////////
		::std::vector<Material*> getMaterials() const
		{
			::std::vector<Material*> result;
			::std::transform(m_materials.begin(), m_materials.end(), ::std::back_inserter(result), [](::std::pair<::std::string, Material*> const& m) -> Material * { return m.second; });
			return result;
		}

		::std::vector<Material**> getEditMaterials()
		{
			::std::vector<Material**> result;
			::std::transform(m_materials.begin(), m_materials.end(), ::std::back_inserter(result), [](::std::pair<::std::string, Material*>  m) -> Material ** { return std::addressof(m.second); });
			return result;
		}
	};

	
}
