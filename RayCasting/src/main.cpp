#include <settings.h>

#include <Geometry/Texture.h>
#include <Math/Vectorf.h>
#include <stdlib.h>
#include <iostream>
#include <Geometry/RGBColor.h>
#include <Geometry/Materials/Material.h>
#include <Geometry/PointLight.h>
#include <Geometry/Camera.h>
#include <Visualizer/Visualizer.h>
#include <Geometry/Scene.h>
#include <Geometry/Loader3ds.h>
#include <Geometry/BoundingBox.h>
#include <chrono>
#include <array>
#include <tbb/tick_count.h>
#include "Geometry/Materials/Phong.h"
#include "Geometry/Materials/CartoonMaterial.h"
#include <Geometry/Materials/Lambert.h>
#include <Geometry/Materials/Glossy.h>
#include <Geometry/Materials/DeltaMirror.h>
#include <Geometry/Materials/Dielectric.h>
#include <Math/Sampler.h>
#include <System/Parallel.h>

#include <Integrators/DirectIntegrator.h>
#include <Integrators/RegularIntegrators.h>
#include <Integrators/RayTracingIntegrator.h>
#include <Integrators/Integrator.h>
#include <Integrators/PathTracingIntegrator.h>
#include <Integrators/ZIntegrator.h>
#include <Integrators/LightIntegrator.h>
#include <Integrators/MISPathTracingIntegrator.h>
#include <Integrators/BidirectionalIntegrator.h>
#include <Integrators/OptimalDirect.h>
#include <Integrators/OptiMISBDPT.h>
#include <Integrators/BlackHolePath.h>
#include <Integrators/PhotonMapper.h>
#include <Integrators/ProgressivePhotonMapper.h>
#include <Integrators/SimpleVCM.h>
#include <Integrators/UncorellatedBDPT.h>
#include <Integrators/VCM.h>
#include <Integrators/OptiVCM.h>

#include <Auto/Auto.h>
#include <Auto/TestScenes.h>

#define ever (;;)



/// <summary>
/// The directory of the 3D objetcs
/// </summary>
const std::string m_modelDirectory = "..\\..\\Models";

using Geometry::RGBColor;

////////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn	void createGround(GeometryCollection::Scene & scene)
///
/// \brief	Adds a ground to the scene.
///
/// \author	F. Lamarche, Université de Rennes 1
/// \date	03/12/2013
///
/// \param [in,out]	scene	The scene.
////////////////////////////////////////////////////////////////////////////////////////////////////
void createGround(Geometry::Scene & scene)
{
	Geometry::BoundingBox sb = scene.getBoundingBox();
	// Non emissive 
	//GeometryCollection::Material * material = new GeometryCollection::Phong(RGBColor(), 0.5, 0.0 , 0.0f, 0, ""); // Non existing material...
	//GeometryCollection::Material * material = new GeometryCollection::Material(RGBColor(), RGBColor(0., 0., 0.), RGBColor(1., 1., 1.), 10000.0f);	//perfect mirror
	//GeometryCollection::Material * material = new GeometryCollection::Material(RGBColor(), RGBColor(1.0f, 1.0f, 1.0f), RGBColor(0.f, 0.f, 0.f), 100.0f); //not a mirror
	//GeometryCollection::Material * material = new GeometryCollection::Material(RGBColor(), RGBColor(1.0,.4,.4)*0.5, RGBColor(), 1000.0f, RGBColor(0.5, 0.5, 0.5)*0.5); // Non existing material...
	//Geometry::Material* material = new Geometry::Glossy(0.9, 1000);
	//Geometry::Material* material = new GeometryCollection::DeltaMirror(0.5);
	Geometry::Material* material = new Geometry::Lambertian<Geometry::REFLECT>(1);

	Geometry::Square * square = new Geometry::Square(material);
	Math::Vector3f scaleV = (sb.max() - sb.min());
	double scale = ::std::max(scaleV[0], scaleV[1])*2.0;
	square->scaleX(scale);
	square->scaleY(scale);
	Math::Vector3f center = (sb.min() + sb.max()) / 2.0;
	center[2] = sb.min()[2];
	square->translate(center);
	scene.add(square);
	::std::cout << "Bounding box: " << sb.min() << "/ " << sb.max() << ", scale: " << scale << ::std::endl;
	::std::cout << "center: " << (sb.min() + sb.max()) / 2.0 << ::std::endl;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn	void createGround(GeometryCollection::Scene & scene)
///
/// \brief	Adds a sirface area light to the scene
///
/// \author	F. Lamarche, Université de Rennes 1
/// \date	03/12/2013
///
/// \param [in,out]	scene	The scene.
////////////////////////////////////////////////////////////////////////////////////////////////////
void createSurfaceLight(Geometry::Scene & scene, double value)
{
	Geometry::BoundingBox sb = scene.getBoundingBox();
	Geometry::Material * material = new Geometry::Material(RGBColor(value, value, value) * 20);
	
	Geometry::Square * square = new Geometry::Square(material);
	Math::Vector3f scaleV = (sb.max() - sb.min());
	//double scale = ::std::max(scaleV[0], scaleV[1])*0.1;
	double factor = 0.5;
	square->scaleX(scaleV[0] * factor);
	square->scaleY(scaleV[1] * factor);
	Math::Vector3f center = (sb.min() + sb.max()) / 2.0;
	center[2] = sb.max()[2] * 3;
	square->translate(center);
	scene.add(square);
}






////////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn	void initDiffuse(GeometryCollection::Scene & scene)
///
/// \brief	Adds a Cornell Box with diffuse material on each walls to the scene. This Cornel box
/// 		contains two cubes.
///
/// \author	F. Lamarche, Université de Rennes 1
/// \date	03/12/2013
///
/// \param [in,out]	scene	The scene.
////////////////////////////////////////////////////////////////////////////////////////////////////
void initDiffuse(Geometry::Scene & scene, Visualizer::Visualizer const& visu)
{
	Geometry::Material * material = new Geometry::Phong(RGBColor(0, 0, 0.0), RGBColor(0.95f, 0.95f, 0.95f), 1, RGBColor());
	Geometry::Material * material2 = new Geometry::Phong(RGBColor(1.0, 1.0, 1.0), RGBColor(0, 0, 0), 1000, RGBColor());
	Geometry::Material * cubeMat = new Geometry::Phong(RGBColor(1, 0.0, 0.0), RGBColor(0.0, 0.0, 0.0), 20.0f, RGBColor());
	//GeometryCollection::Material * cubeMat = new GeometryCollection::Material(RGBColor(), RGBColor(1.0f,0.0,0.0), RGBColor(0.0,0.0,0.0), 20.0f, RGBColor(10.0,0,0)) ;
	Geometry::Cornel::init_cornell(scene, material2, material2, material, material, material, material, 10);

	Geometry::Cube * tmp = new Geometry::Cube(cubeMat);
	tmp->translate(Math::makeVector(1.5, -1.5, 0.0));
	scene.add(tmp);

	Geometry::Cube * tmp2 = new Geometry::Cube(cubeMat);
	tmp2->translate(Math::makeVector(2, 1, -4));
	scene.add(tmp2);

	

	// 2.2 Adds point lights in the scene 
	{
		Geometry::PointLight pointLight(Math::makeVector(0.0f, 0.f, 2.0f), RGBColor(0.5f, 0.5f, 0.5f));
		scene.add(pointLight);
	}
	{
		Geometry::PointLight pointLight2(Math::makeVector(4.f, 0.f, 0.f), RGBColor(0.5f, 0.5f, 0.5f));
		scene.add(pointLight2);
	}
	{
		Geometry::Camera camera(Math::makeVector(-4.0f, 0.0f, 0.0f), Math::makeVector(0.0f, 0.0f, 0.0f), 0.3f, ((double)visu.width())/((double)visu.height()), 1.0f);
		scene.setCamera(camera);
	}
}




//template <class _Phong=GeometryCollection::Phong>
//void initDiffuse_surface_light(GeometryCollection::Scene & scene, Visualizer::Visualizer const& visu)
//{
//	GeometryCollection::Material * material_green = new _Phong(RGBColor(0, 1, .0), RGBColor(0), 0, RGBColor());
//	GeometryCollection::Material * material_blue = new _Phong(RGBColor(0, 0, 1), RGBColor(0), 0, RGBColor());
//	GeometryCollection::Material * material_yellow = new _Phong(RGBColor(1, 1, 0.3), RGBColor(0), 0, RGBColor());
//	GeometryCollection::Material * material_purple = new _Phong(RGBColor(0.5, 0.2, 0.8)*1.2, RGBColor(0), 0, RGBColor());
//	GeometryCollection::Material * material2 = new _Phong(0.7, 0, 0);
//	GeometryCollection::Material * material3 = new _Phong(0.0, 0.9, 100);
//	GeometryCollection::Material * material4 = new _Phong(0, 0.9, 10000);
//	GeometryCollection::Material * cubeMat = new _Phong(RGBColor(1, 0.0, 0.0), RGBColor(0.0, 0.0, 0.0), 0.0f, RGBColor());
//	//GeometryCollection::Material * cubeMat = new GeometryCollection::Material(RGBColor(), RGBColor(1.0f,0.0,0.0), RGBColor(0.0,0.0,0.0), 20.0f, RGBColor(10.0,0,0)) ;
//	GeometryCollection::Cornel::init_cornell(scene, material2, material2, material2, material2, material2, material2, 10);
//
//
//	GeometryCollection::Cube * tmp = new GeometryCollection::Cube(cubeMat);
//	tmp->scale(2);
//	tmp->translate(Math::makeVector(2.5, -2.5, 0.0));
//	scene.add(tmp);
//
//	GeometryCollection::Cube * tmp2 = new GeometryCollection::Cube(cubeMat);
//	tmp2->translate(Math::makeVector(2, 1, -4));
//	scene.add(tmp2);
//
//
//	GeometryCollection::Sphere s({ 0, -2, -4 }, 1, new _Phong(0, 0, 0.7, 100000));
//	//scene.add(s);
//
//	GeometryCollection::Sphere s2({ 0, 2, -4 }, 1, new _Phong(0, 0, 0.7, 10000));
//	//scene.add(s2);
//
//
//	// 2.2 Adds point lights in the scene 
//	{
//		GeometryCollection::PointLight pointLight(Math::makeVector(0.0f, 0.f, 2.0f), RGBColor(0.5f, 0.5f, 0.5f));
//		scene.add(pointLight);
//	}
//	{
//		GeometryCollection::PointLight pointLight2(Math::makeVector(4.f, 0.f, 0.f), RGBColor(0.5f, 0.5f, 0.5f));
//		scene.add(pointLight2);
//	}
//	//Add the surface light
//	///*
//	{
//		GeometryCollection::Material * surface = new _Phong(0, 0, 0, 0, 1);
//		GeometryCollection::Square * slight = new GeometryCollection::Square(surface);
//		slight->scaleY(3);
//		slight->scaleX(3);
//		slight->translate(Math::makeVector(0.0, 0.0, 4.9));
//		scene.add(slight);
//	}
//	//*/
//
//	/*
//	{
//		GeometryCollection::Material * surface = new GeometryCollection::Material(0, 0, 0, 0, 0.5);
//		GeometryCollection::OneTriangle slight(surface);
//		slight.scaleY(4.0);
//		slight.scaleX(4.0);
//		slight.translate(Math::makeVector(0.0, 0.0, -4.9));
//		scene.add(slight);
//	}
//	//*/
//
//
//	//Add the Camera
//	{
//		GeometryCollection::Camera camera(Math::makeVector(-4.0f, 0.0f, 0.0f), Math::makeVector(1.0f, 0.0f, 0.0f), 0.3f, ((double)visu.width()) / ((double)visu.height()), 1.0f);
//		scene.setCamera(camera);
//	}
//
//	
//}







void initDiffuse_dif(Geometry::Scene & scene, Visualizer::Visualizer const& visu)
{
	Geometry::Material * material = new Geometry::Phong(0, 0.8, 1000, RGBColor());
	Geometry::Material * material2 = new Geometry::Phong(0.5, 0, 0, RGBColor());
	Geometry::Material * cubeMat = new Geometry::Phong(RGBColor(1, 0.0, 0.0), RGBColor(0.0, 0.0, 0.0), 20.0f, RGBColor());
	//GeometryCollection::Material * cubeMat = new GeometryCollection::Material(RGBColor(), RGBColor(1.0f,0.0,0.0), RGBColor(0.0,0.0,0.0), 20.0f, RGBColor(10.0,0,0)) ;
	Geometry::Cornel::init_cornell(scene, material2, material2, material, material, material, material, 10);


	Geometry::Cube * tmp= new Geometry::Cube(cubeMat);
	tmp->translate(Math::makeVector(1.5, -1.5, 0.0));
	scene.add(tmp);

	Geometry::Cube * tmp2 = new Geometry::Cube(cubeMat);
	tmp2->translate(Math::makeVector(2, 1, -4));
	scene.add(tmp2);

	//GeometryCollection::Phong * transparent_mat = new GeometryCollection::Phong(RGBColor(), RGBColor(), RGBColor(0.95, 0.95, 0.95), 1000);
	//transparent_mat->m_transparent = true;

	//GeometryCollection::CylinderFabrice transparent_obj(20, 0.5, 0.5, transparent_mat);
	//transparent_obj.scale(1.5);
	//transparent_obj.m_medium = GeometryCollection::medium(2.0);
	//transparent_mat->m_medium = transparent_obj.get_medium();
	//scene.add(transparent_obj);


	{
		Geometry::Material* emissive = new Geometry::Material(0.5);
		Geometry::Sphere sphere({ -4, 0, 4 }, 1, emissive);

		scene.add(sphere);
	}
	// 2.2 Adds point lights in the scene 
	{
		Geometry::PointLight pointLight(Math::makeVector(-4.0f, 0.f, 4.0f), 0.5);
		//scene.add(pointLight);
	}
	{
		Geometry::PointLight pointLight2(Math::makeVector(4.f, 0.f, 0.f), 0.1);
		//scene.add(pointLight2);
	}
	{
		Geometry::Camera camera(Math::makeVector(-4.0f, 0.0f, 0.0f), Math::makeVector(0.0f, 0.0f, 0.0f), 0.3f, ((double)visu.width()) / ((double)visu.height()), 1.0f);
		scene.setCamera(camera);
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn	void initSpecular(GeometryCollection::Scene & scene)
///
/// \brief	Adds a Cornel box in the provided scene. Walls are specular and the box contains two 
/// 		cubes.
///
/// \author	F. Lamarche, Université de Rennes 1
/// \date	03/12/2013
////////////////////////////////////////////////////////////////////////////////////////////////////
void initSpecular(Geometry::Scene & scene, Visualizer::Visualizer const& visu)
{
	Geometry::Material * material = new Geometry::Phong(RGBColor(0, 0, 0.0), RGBColor(0.7f, 0.7f, 0.7f), 100, RGBColor());
	Geometry::Material * material2 = new Geometry::Phong(RGBColor(0, 0, 1.0f), RGBColor(0, 0, 0), 1000, RGBColor());
	//GeometryCollection::Material * cubeMat = new GeometryCollection::Material(RGBColor(), RGBColor(1.0f,0.0,0.0), RGBColor(0.0,0.0,0.0), 20.0f, RGBColor(10.0,0,0)) ;
	Geometry::Material * cubeMat = new Geometry::Phong(RGBColor(1.0f, 0.0, 0.0), RGBColor(0.0, 0.0, 0.0), 20.0f, RGBColor());
	
	
	Geometry::Cornel::init_cornell(scene, material, material, material, material, material, material, 10); //new GeometryCollection::Cube(material2) ;////new Cone(4, material) ; //new GeometryCollection::CylinderFabrice(5, 1, 1, material) ;////////new GeometryCollection::Cube(material) ;////; //new GeometryCollection::Cube(material) ; //new GeometryCollection::CylinderFabrice(100, 2, 1, material) ; //


	Geometry::Cube* tmp=new Geometry::Cube(cubeMat);
	tmp->translate(Math::makeVector(1.5, -1.5, 0.0));
	scene.add(tmp);

	Geometry::Cube * tmp2 = new Geometry::Cube(cubeMat);
	tmp2->translate(Math::makeVector(2, 1, -4));
	scene.add(tmp2);

	// 2.2 Adds point lights in the scene 
	{
		Geometry::PointLight pointLight(Math::makeVector(0.0f, 0.f, 2.0f), RGBColor(0.5f, 0.5f, 0.5f) * 5);
		scene.add(pointLight);
	}
	{
		Geometry::PointLight pointLight2(Math::makeVector(4.f, 0.f, 0.f), RGBColor(0.5f, 0.5f, 0.5f) * 5);
		scene.add(pointLight2);
	}
	// Sets the camera
	{
		Geometry::Camera camera(Math::makeVector(-4.0f, 0.0f, 0.0f), Math::makeVector(0.0f, 0.0f, 0.0f), 0.3f, ((double)visu.width()) / ((double)visu.height()), 1.0f);
		scene.setCamera(camera);
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn	void initDiffuseSpecular(GeometryCollection::Scene & scene)
///
/// \brief	Adds a Cornel box in the provided scene. The cornel box as diffuse and specular walls and 
/// 		contains two boxes.
///
/// \author	F. Lamarche, Université de Rennes 1
/// \date	03/12/2013
///
/// \param [in,out]	scene	The scene.
////////////////////////////////////////////////////////////////////////////////////////////////////
void initDiffuseSpecular(Geometry::Scene & scene, Visualizer::Visualizer const& visu)
{
	Geometry::Material * material = new Geometry::Phong(RGBColor(0, 0, 0.0), RGBColor(0.7f, 0.71f, 0.7f), 1000, RGBColor());
	Geometry::Material * material2 = new Geometry::Phong(RGBColor(1, 1, 1.0f), RGBColor(0, 0, 0), 1000, RGBColor());
	Geometry::Material * cubeMat = new Geometry::Phong(RGBColor(0., 0., 0.), RGBColor(0.7, 0.4, 0.), 20.0f, RGBColor());//gold
	Geometry::Material * cubeMat3 = new Geometry::Phong(RGBColor(0., 0., 0.), RGBColor(0.5, 0.5, 0.5), 20.0f, RGBColor());//silver
	Geometry::Material * cubeMat2 = new Geometry::Phong(RGBColor(1.0f, 0.0, 0.0), RGBColor(0.0, 0.0, 0.0), 20.0f, RGBColor(1));//red
																																			 //GeometryCollection::Material * cubeMat = new GeometryCollection::Material(RGBColor(), RGBColor(0.0f,0.0,0.0), RGBColor(0.0,0.0,0.0), 20.0f, RGBColor(10.0,0,0)) ;
																																			 //GeometryCollection::Material * cubeMat2 = new GeometryCollection::Material(RGBColor(), RGBColor(0.0f,0.0,0.0), RGBColor(0.0,0.0,0.0), 20.0f, RGBColor(0.0,10,0)) ;
	Geometry::Cornel::init_cornell(scene, material2, material2, material, material, material, material, 10); //new GeometryCollection::Cube(material2) ;////new Cone(4, material) ; //new GeometryCollection::CylinderFabrice(5, 1, 1, material) ;////////new GeometryCollection::Cube(material) ;////; //new GeometryCollection::Cube(material) ; //new GeometryCollection::CylinderFabrice(100, 2, 1, material) ; //

	Geometry::CylinderFabrice * cyl = new Geometry::CylinderFabrice(100, 1., -1., cubeMat3);
	cyl->translate(Math::makeVector(0., -2., -2.));
	
	scene.add(cyl);


	Geometry::Cube * tmp = new Geometry::Cube(cubeMat2);
	tmp->translate(Math::makeVector(1.5, -1.5, 0.0));
	scene.add(tmp);

	Geometry::Cube * tmp2 = new Geometry::Cube(cubeMat);
	tmp2->scaleY(5);
	tmp2->scaleX(4);
	tmp2->translate(Math::makeVector(2, 1, -4));

	scene.add(tmp2);

	// 2.2 Adds point lights in the scene 
	{
		Geometry::PointLight pointLight(Math::makeVector(0.0f, 0.f, 2.0f), RGBColor(0.5f, 0.5f, 0.5f)*5);
		scene.add(pointLight);
	}
	{
		Geometry::PointLight pointLight2(Math::makeVector(4.f, 0.f, 0.f), RGBColor(0.5f, 0.5f, 0.5f)*5);
		scene.add(pointLight2);
	}
	{
		Geometry::Camera camera(Math::makeVector(-4.f, .0f, .0f), Math::makeVector(0.0f, 0.0f, 0.0f), 0.4f, ((double)visu.width()) / ((double)visu.height()), 1.0f);
		scene.setCamera(camera);
	}

	/*
	GeometryCollection::Material * black = new GeometryCollection::Phong(RGBColor(), RGBColor(0.1, 0.1, 0.1), RGBColor(0., 0., 0.), 10, RGBColor());
	{
		GeometryCollection::Cube pillar(black);

		pillar.scaleX(.1);
		pillar.scaleY(.1);
		pillar.scaleZ(10);

		pillar.translate(Math::makeVector(5., 5., 0.));
		scene.add(pillar);

		pillar.translate(Math::makeVector(0., -10., 0.));
		scene.add(pillar);

		pillar.translate(Math::makeVector(-10., 0., 0.));
		scene.add(pillar);

		pillar.translate(Math::makeVector(0., 10., 0.));
		scene.add(pillar);
	}
	{
		GeometryCollection::Cube pillar(black);

		pillar.scaleX(10);
		pillar.scaleY(.1);
		pillar.scaleZ(.1);

		pillar.translate(Math::makeVector(0., 5., 5.));
		scene.add(pillar);

		pillar.translate(Math::makeVector(0., 0., -10.));
		scene.add(pillar);


		pillar.translate(Math::makeVector(0., -10., 0.));
		scene.add(pillar);

		pillar.translate(Math::makeVector(0., 0., 10.));
		scene.add(pillar);
	}
	{
		GeometryCollection::Cube pillar(black);

		pillar.scaleX(0.1);
		pillar.scaleY(10);
		pillar.scaleZ(.1);

		pillar.translate(Math::makeVector(5., 0., 5.));
		scene.add(pillar);

		pillar.translate(Math::makeVector(0., 0., -10.));
		scene.add(pillar);


		pillar.translate(Math::makeVector(-10., 0., 0.));
		scene.add(pillar);

		pillar.translate(Math::makeVector(0., 0., 10.));
		scene.add(pillar);
	}
	*/
}










/// <summary>
/// Intializes a scene containing a garage.
/// </summary>
/// <param name="scene"></param>
void initGarage(Geometry::Scene & scene, Visualizer::Visualizer const& visu)
{
	Geometry::Loader3ds loader(m_modelDirectory + "\\garage.3ds", m_modelDirectory + "");

	for (size_t cpt = 0; cpt < loader.getMeshes().size(); ++cpt)
	{
		scene.add(new Geometry::GeometryCollection(*loader.getMeshes()[cpt]));
	}

	// 2.2 Adds point lights in the scene 
	Geometry::BoundingBox sb = scene.getBoundingBox();
	Math::Vector3f position = sb.max();
	Geometry::PointLight light1(position, RGBColor(1.0, 1.0, 1.0) * 500);
	scene.add(light1);
	position[1] = -position[1];
	Geometry::PointLight light2(position, RGBColor(1.0, 1.0, 1.0) * 500);
	scene.add(light2);
	{
		Geometry::Camera camera(Math::makeVector(750.0f, -1500.f, 1000.f)*0.85f, Math::makeVector(200.0f, 0.0f, 0.0f), 0.3f, ((double)visu.width()) / ((double)visu.height()), 1.0f);
		scene.setCamera(camera);
	}
	createGround(scene);
}

/// <summary>
/// Initializes a scene containing a guitar.
/// </summary>
/// <param name="scene"></param>
void initGuitar(Geometry::Scene & scene, Visualizer::Visualizer const& visu)
{
	Geometry::Loader3ds loader(m_modelDirectory + "\\guitar2.3ds", m_modelDirectory + "");

	for (size_t cpt = 0; cpt < loader.getMeshes().size(); ++cpt)
	{
		scene.add(new Geometry::GeometryCollection(*loader.getMeshes()[cpt]));
	}
	// 2.2 Adds point lights in the scene 
	Geometry::BoundingBox sb = scene.getBoundingBox();
	Math::Vector3f position = sb.max();
	Geometry::PointLight light1(position + Math::makeVector(0.f, 0.f, 70.f/*100.f*/), RGBColor(1.0, 1.0, 1.0)*200 );
	scene.add(light1);
	position[1] = -position[1];
	Geometry::PointLight light2(position + Math::makeVector(0.f, 0.f, 200.f), RGBColor(1.0, 1.0, 1.0)*200 );
	scene.add(light2);
	{
		Geometry::Camera camera(Math::makeVector(-500., -1000., 1000.)*1.05, Math::makeVector(500.f, 0.0f, 0.0f), 0.6f, ((double)visu.width()) / ((double)visu.height()), 1.0f);
		camera.translateLocal(Math::makeVector(100.0f, -100.f, -200.f));
		scene.setCamera(camera);
	}
	createGround(scene);
}

/// <summary>
/// Initializes a scene containing a dog
/// </summary>
/// <param name="scene"></param>
void initDog(Geometry::Scene & scene, Visualizer::Visualizer const& visu)
{
	Geometry::Loader3ds loader(m_modelDirectory + "\\Dog\\dog.3ds", m_modelDirectory + "\\dog");

	Geometry::Material* glass = new Geometry::Dielectric({ 1, 0.8, 0.9 }, 1.25);
	Geometry::Material* Lglass = new Geometry::Lambertian<Geometry::TRANSMIT>({ 1, 0.8, 0.9 });

	for (size_t cpt = 0; cpt < loader.getMeshes().size(); ++cpt)
	{
		Geometry::GeometryCollection * obj = new Geometry::GeometryCollection(*loader.getMeshes()[cpt]);
		obj->set_material(glass);
		scene.add(obj);
	}

	// 2.2 Adds point lights in the scene 
	Geometry::BoundingBox sb = scene.getBoundingBox();
	Math::Vector3f position = sb.max();
	Geometry::PointLight light1(position, RGBColor(1.0, 1.0, 1.0) *.5);
	scene.add(light1);
	position[1] = -position[1];
	Geometry::PointLight light2(position, RGBColor(1.0, 1.0, 1.0) *.5);

	Geometry::Material* lam = new Geometry::Lambertian<Geometry::REFLECT>(0.9);
	for (int i = -1; i < 2; i += 2)
		for(int j=-1; j<2; j+=2)
	{
			Geometry::Sphere s = Geometry::Sphere({5 * (i-j), 5*(i+j), 2}, 2, lam);
			scene.add(s);
	}
	Geometry::Material * emisive = new Geometry::Material(RGBColor( 2, 2, 1.8) * 3, "");
	Geometry::Square * surface = new Geometry::Square(emisive);
	surface->scale(5.0);
	surface->translate(Math::makeVector(0, 0, 10));
	scene.add(surface);

	scene.add(light2);
	{
		Geometry::Camera camera(Math::makeVector(10.f, 10.f, 6.f)*0.5, Math::makeVector(0.f, 0.0f, 2.5f), .7f, ((double)visu.width()) / ((double)visu.height()), 1.0f);
		camera.translateLocal(Math::makeVector(-0.4f, 0.f, 0.9f));
		scene.setCamera(camera);
	}
	createGround(scene);
}

/// <summary>
/// Initializes a scene containing a temple.
/// </summary>
/// <param name="scene"></param>
void initTemple(Geometry::Scene & scene, Visualizer::Visualizer const& visu)
{
	Geometry::BoundingBox box(Math::makeVector(0.0f, 0.0f, 0.0f), Math::makeVector(0.0f, 0.0f, 0.0f));
	Geometry::Loader3ds loader(m_modelDirectory + "\\Temple\\Temple of St Seraphim of Sarov N270116_2.3ds", m_modelDirectory + "\\Temple");

	for (size_t cpt = 0; cpt < loader.getMeshes().size(); ++cpt)
	{
		//loader.getMeshes()[cpt]->translate(Math::makeVector(20.f, 0.f, 40.0f));
		scene.add(new Geometry::GeometryCollection(*loader.getMeshes()[cpt]));
		box.update(loader.getMeshes()[cpt]->box());
	}


	// 2.2 Adds point lights in the scene 
	Geometry::BoundingBox sb = scene.getBoundingBox();
	Math::Vector3f position = sb.max();
	Geometry::PointLight light1(position, RGBColor(1.0, 1.0, 1.0)*20.0);
	//GeometryCollection::PointLight light1(position, RGBColor(1.0, 1.0, 1.0));
	scene.add(light1);
	position[1] = -position[1];
	Geometry::PointLight light2(position, RGBColor(1.0, 1.0, 1.0)*10);
	//GeometryCollection::PointLight light2(position, RGBColor(1.0, 1.0, 1.0) );
	scene.add(light2);
	{
		Geometry::Camera camera(Math::makeVector(20.0f, -100.0f, 15.0f), Math::makeVector(-20.f, 0.f, -40.f), 0.3f, ((double)visu.width()) / ((double)visu.height()), 1.0f);
		camera.translateLocal(Math::makeVector(40.f, 0.f, 0.f));
		scene.setCamera(camera);
	}
	createGround(scene);
	//createSurfaceLight(scene, 1);
	createSurfaceLight(scene, 1);
}

/// <summary>
/// Initializes a scene containing a robot.
/// </summary>
/// <param name="scene"></param>
void initRobot(Geometry::Scene & scene, Visualizer::Visualizer const& visu)
{
	Geometry::BoundingBox box(Math::makeVector(0.0f, 0.0f, 0.0f), Math::makeVector(0.0f, 0.0f, 0.0f));
	Geometry::Loader3ds loader(m_modelDirectory + "\\Robot.3ds", m_modelDirectory + "");

	for (size_t cpt = 0; cpt < loader.getMeshes().size(); ++cpt)
	{
		//loader.getMeshes()[cpt]->translate(Math::makeVector(20.f, 0.f, 40.0f));
		scene.add(new Geometry::GeometryCollection(*loader.getMeshes()[cpt]));
		box.update(loader.getMeshes()[cpt]->box());
	}

	// 2.2 Adds point lights in the scene 
	Geometry::BoundingBox sb = scene.getBoundingBox();
	Math::Vector3f position = sb.max();
	Geometry::PointLight light1(position, RGBColor(1.0, 1.0, 1.0) * 10);
	scene.add(light1);
	position[1] = -position[1];
	Geometry::PointLight light2(position, RGBColor(1.0, 1.0, 1.0) *20);
	scene.add(light2);
	{
		Geometry::Camera camera(Math::makeVector(100.0f, -50.0f, 0.0f), Math::makeVector(0.f, 0.f, 5.f), 0.3f, ((double)visu.width()) / ((double)visu.height()), 1.0f);
		camera.translateLocal(Math::makeVector(0.f, 40.f, 50.f));
		scene.setCamera(camera);
	}
	createGround(scene);
}

/// <summary>
/// Initializes a scene containing a grave stone
/// </summary>
/// <param name="scene"></param>
void initGraveStone(Geometry::Scene & scene, Visualizer::Visualizer const& visu)
{
	Geometry::BoundingBox box(Math::makeVector(0.0f, 0.0f, 0.0f), Math::makeVector(0.0f, 0.0f, 0.0f));
	Geometry::Loader3ds loader(m_modelDirectory + "\\gravestone\\GraveStone.3ds", m_modelDirectory + "\\gravestone");

	::std::vector<Geometry::Material*> materials = loader.getMaterials();
	for (auto it = materials.begin(), end = materials.end(); it != end; ++it)
	{
		//(*it)->setAmbient(RGBColor());
		//(*it)->setSpecular((*it)->getSpecular()*0.05);
	}

	for (size_t cpt = 0; cpt < loader.getMeshes().size(); ++cpt)
	{
		//loader.getMeshes()[cpt]->translate(Math::makeVector(20.f, 0.f, 40.0f));
		scene.add(new Geometry::GeometryCollection(*loader.getMeshes()[cpt]));
		box.update(loader.getMeshes()[cpt]->box());
	}

	// 2.2 Adds point lights in the scene 
	Geometry::BoundingBox sb = scene.getBoundingBox();
	Math::Vector3f position = sb.max();
	Geometry::PointLight light1(position, RGBColor(1.0, 1.0, 1.0)*100 );
	scene.add(light1);
	position[1] = -position[1];
	Geometry::PointLight light2(position, RGBColor(1.0, 1.0, 1.0)*100 );
	scene.add(light2);
	{
		Geometry::Camera camera(Math::makeVector(0.f, -300.0f, 200.0f), Math::makeVector(0.f, 0.f, 125.f), 0.3f, ((double)visu.width()) / ((double)visu.height()), 1.0f);
		camera.translateLocal(Math::makeVector(0.f, 80.f, 120.f));
		scene.setCamera(camera);
	}
	createGround(scene);
	//createSurfaceLight(scene, 400);
}

/// <summary>
/// Initializes a scene containing a boat
/// </summary>
/// <param name="scene"></param>
void initBoat(Geometry::Scene & scene, Visualizer::Visualizer const& visu)
{
	Geometry::BoundingBox box(Math::makeVector(0.0f, 0.0f, 0.0f), Math::makeVector(0.0f, 0.0f, 0.0f));
	Geometry::Loader3ds loader(m_modelDirectory + "\\Boat\\boat.3ds", m_modelDirectory + "\\boat");
	//// We remove the specular components of the materials...
	::std::vector<Geometry::Material*> materials = loader.getMaterials();
	for (auto it = materials.begin(), end = materials.end(); it != end; ++it)
	{
		//(*it)->setAmbient(RGBColor(0.01, 0.01, 0.01));
	}

	for (size_t cpt = 0; cpt < loader.getMeshes().size(); ++cpt)
	{
		//loader.getMeshes()[cpt]->translate(Math::makeVector(20.f, 0.f, 40.0f));
		scene.add(new Geometry::GeometryCollection(*loader.getMeshes()[cpt]));
		box.update(loader.getMeshes()[cpt]->box());
	}

	// 2.2 Adds point lights in the scene 
	Geometry::BoundingBox sb = scene.getBoundingBox();
	Math::Vector3f position = sb.max();
	Geometry::PointLight light1(position, RGBColor(1.0, 1.0, 1.0)*6000);
	scene.add(light1);
	position[1] = -position[1];
	Geometry::PointLight light2(position, RGBColor(1.0, 1.0, 1.0) * 20000);
	scene.add(light2);
	{
		Geometry::Camera camera(Math::makeVector(5000.f, 5000.f, 200.0f), Math::makeVector(0.f, 0.f, 60.f), 0.3f, ((double)visu.width()) / ((double)visu.height()), 1.0f);
		camera.translateLocal(Math::makeVector(2000.f, 3000.f, 700.f));
		scene.setCamera(camera);
	}
	createGround(scene);
}

/// <summary>
/// Initializes a scene containing a tibet house
/// </summary>
/// <param name="scene"></param>
void initTibetHouse(Geometry::Scene & scene, Visualizer::Visualizer const& visu)
{
	Geometry::BoundingBox box(Math::makeVector(0.0f, 0.0f, 0.0f), Math::makeVector(0.0f, 0.0f, 0.0f));
	Geometry::Loader3ds loader(m_modelDirectory + "\\TibetHouse\\TibetHouse.3ds", m_modelDirectory + "\\TibetHouse");
	// We remove the specular components of the materials... and add surface light sources :)
	::std::vector<Geometry::Material*> materials = loader.getMaterials();
	for (auto it = materials.begin(), end = materials.end(); it != end; ++it)
	{
		//((GeometryCollection::Phong*)(*it))->setSpecular(RGBColor());
		//((GeometryCollection::Phong*)(*it))->setAmbient(RGBColor());
		//std::cout << (*it)->getAmbient() << std::endl;
		if ((*it)->getTextureFile() == m_modelDirectory + "\\TibetHouse" + "\\3D69C2DE.png")
		{
			((*it))->setEmissive(RGBColor(1.0, 1.0, 1.0)*10*5);
			
			//((GeometryCollection::Phong*)(*it))->setDiffuse(RGBColor(0, 0, 0));
		}
	}

	for (size_t cpt = 0; cpt < loader.getMeshes().size(); ++cpt)
	{
		//loader.getMeshes()[cpt]->translate(Math::makeVector(20.f, 0.f, 40.0f));
		scene.add(new Geometry::GeometryCollection(*loader.getMeshes()[cpt]));
		box.update(loader.getMeshes()[cpt]->box());
	}

	// 2.2 Adds point lights in the scene 
	Geometry::BoundingBox sb = scene.getBoundingBox();
	Math::Vector3f position = sb.max();
	Geometry::PointLight light1(position, RGBColor(1.0, 1.0, 1.0) * 20);
	scene.add(light1);
	position[1] = -position[1];
	Geometry::PointLight light2(position, RGBColor(1.0, 1.0, 1.0) *20);
	scene.add(light2);
	Geometry::PointLight light3(Math::makeVector(5.f, 35.f, 5.f), RGBColor(1.0, 1.0, 1.0)*20 ); //*50
	scene.add(light3);
	{
		Geometry::Camera camera(Math::makeVector(20.f, 0.f, 0.0f), Math::makeVector(5.f, 35.f, 0.f), 0.3f, ((double)visu.width()) / ((double)visu.height()), 1.0f);
		camera.translateLocal(Math::makeVector(0.f, 5.f, 0.f));// +Math::makeVector(0.0, 5.0, 0.0));
		scene.setCamera(camera);
	}
	createGround(scene);
}

/// <summary>
/// Initializes a scene containing a tibet house. Camera is placed inside the house.
/// </summary>
/// <param name="scene"></param>
void initTibetHouseInside(Geometry::Scene & scene, Visualizer::Visualizer const& visu)
{
	Geometry::BoundingBox box(Math::makeVector(0.0f, 0.0f, 0.0f), Math::makeVector(0.0f, 0.0f, 0.0f));
	Geometry::Loader3ds loader(m_modelDirectory + "\\TibetHouse\\TibetHouse.3ds", m_modelDirectory + "\\TibetHouse");
	// We remove the specular components of the materials... and add surface light sources :)
	::std::vector<Geometry::Material*> materials = loader.getMaterials();
	for (auto it = materials.begin(), end = materials.end(); it != end; ++it)
	{
		
		//(*it)->setAmbient(RGBColor());
		if ((*it)->getTextureFile() == m_modelDirectory + "\\TibetHouse" + "\\3D69C2DE.png")
		{
			((*it))->setEmissive(RGBColor(1.0, .2, 0.0));
			
		}
	}

	for (size_t cpt = 0; cpt < loader.getMeshes().size(); ++cpt)
	{
		//loader.getMeshes()[cpt]->translate(Math::makeVector(20.f, 0.f, 40.0f));
		scene.add(new Geometry::GeometryCollection(*loader.getMeshes()[cpt]));
		box.update(loader.getMeshes()[cpt]->box());
	}

	// 2.2 Adds point lights in the scene 
	Geometry::BoundingBox sb = scene.getBoundingBox();
	Math::Vector3f position = sb.max();
	Geometry::PointLight light1(position, RGBColor(1.0, 1.0, 1.0)*10);
	scene.add(light1);
	position[1] = -position[1];
	Geometry::PointLight light2(position, RGBColor(1.0, 1.0, 1.0)*10);
	scene.add(light2);
	Geometry::PointLight light3(Math::makeVector(5.f, 35.f, 5.f), RGBColor(1.0, 1.0, 1.0)*10 );
	scene.add(light3);
	{
		Geometry::Camera camera(Math::makeVector(20.f, 0.f, 5.0f), Math::makeVector(5.f, 35.f, 5.f), 0.3f, ((double)visu.width()) / ((double)visu.height()), 1.0f);
		camera.translateLocal(Math::makeVector(10.f - 2, 30.f, 0.f));
		scene.setCamera(camera);
	}
	createGround(scene);
}

/// <summary>
/// Initializes a scene containing a medieval city
/// </summary>
/// <param name="scene"></param>
void initMedievalCity(Geometry::Scene & scene, Visualizer::Visualizer const& visu)
{
	Geometry::BoundingBox box(Math::makeVector(0.0f, 0.0f, 0.0f), Math::makeVector(0.0f, 0.0f, 0.0f));
	Geometry::Loader3ds loader(m_modelDirectory + "\\Medieval\\MedievalCity.3ds", m_modelDirectory + "\\Medieval\\texture");
	// We remove the specular components of the materials...
	::std::vector<Geometry::Material*> materials = loader.getMaterials();
	for (auto it = materials.begin(), end = materials.end(); it != end; ++it)
	{
		//((GeometryCollection::Phong*)(*it))->setSpecular(RGBColor());
		//((GeometryCollection::Phong*)(*it))->setAmbient(RGBColor());

	}

	for (size_t cpt = 0; cpt < loader.getMeshes().size(); ++cpt)
	{
		//loader.getMeshes()[cpt]->translate(Math::makeVector(20.f, 0.f, 40.0f));
		//loader.getMeshes()[cpt]->rotate(Math::Quaternion<double>(Math::makeVector(0.0, 1.0, 0.0), Math::pi));
		scene.add(new Geometry::GeometryCollection(*loader.getMeshes()[cpt]));
		box.update(loader.getMeshes()[cpt]->box());
	}

	// 2.2 Adds point lights in the scene 
	Geometry::BoundingBox sb = scene.getBoundingBox();
	Math::Vector3f position = sb.max();
	Geometry::PointLight light1(position + Math::makeVector(0.0, 0.0, 150.0), RGBColor(1.0, 0.6, 0.3) * 100);
	scene.add(light1);
	position[1] = -position[1];
	Geometry::PointLight light2(position + Math::makeVector(0.0, 0.0, 1000.0), RGBColor(1.0, 1.0, 1.0) * 100);
	//scene.add(light2);
	Geometry::PointLight light3(Math::makeVector(5.f, 35.f, 5.f), RGBColor(1.0, 1.0, 1.0) * 50);
	//scene.add(light3);

	{
		//GeometryCollection::Material * sun_material = new GeometryCollection::Material(RGBColor(1, 0.6, 0.3) * 2000);
		Geometry::Material * sun_material = new Geometry::Material(RGBColor(1, 0.6, 0.3) * 2000);

		Geometry::Sphere sun = Geometry::Sphere(Math::Vector3f(1000, 1000, 1000), 100, sun_material);

		scene.add(sun);
	}

	{
		//GeometryCollection::Camera camera(Math::makeVector(0.f, 300.f, 1000.0f), Math::makeVector(0.0, 0.0, 0.0), 0.3f, ((double)visu.width()) / ((double)visu.height()), 1.0f);
		Geometry::Camera camera({ 47.661, 15.3808, 44.9518 }, { 48.3153, 16.094, 44.7 }, 0.3f, ((double)visu.width()) / ((double)visu.height()), 1.0f);
		
		//camera.translateLocal(Math::makeVector(0.0, 800., -100.0));
		scene.setCamera(camera);
	}
	createGround(scene);
}

/// <summary>
/// Initializes a scene containing a sombrero
/// </summary>
/// <param name="scene"></param>
void initSombrero(Geometry::Scene & scene, Visualizer::Visualizer const& visu)
{
	Geometry::BoundingBox box(Math::makeVector(0.0f, 0.0f, 0.0f), Math::makeVector(0.0f, 0.0f, 0.0f));
	Geometry::Loader3ds loader(m_modelDirectory + "\\sombrero\\sombrero.3ds", m_modelDirectory + "\\sombrero");
	// We remove the specular components of the materials...
	::std::vector<Geometry::Material*> materials = loader.getMaterials();
	for (auto it = materials.begin(), end = materials.end(); it != end; ++it)
	{
		//((GeometryCollection::Phong*)(*it))->setSpecular(RGBColor());
		//((GeometryCollection::Phong*)(*it))->setAmbient(RGBColor());
	}

	for (size_t cpt = 0; cpt < loader.getMeshes().size(); ++cpt)
	{
		scene.add(new Geometry::GeometryCollection(*loader.getMeshes()[cpt]));
		box.update(loader.getMeshes()[cpt]->box());
	}

	// 2.2 Adds point lights in the scene 
	Geometry::BoundingBox sb = scene.getBoundingBox();
	Math::Vector3f position = sb.max();
	Geometry::PointLight light1(position, RGBColor(1.0, 1.0, 1.0)*100);
	scene.add(light1);
	position[1] = -position[1];
	Geometry::PointLight light2(position, RGBColor(1.0, 1.0, 1.0)*100);
	scene.add(light2);

	Geometry::Material * emisive = new Geometry::Material(1);
	Geometry::Square * surface = new Geometry::Square(emisive);
	surface->scale(500);
	surface->translate(Math::makeVector(0, 0, 400));
	scene.add(surface);
	
	{
		Geometry::Camera camera(Math::makeVector(300.f, 0.f, 100.0f), Math::makeVector(0.f, 0.f, 0.f), 0.3f, ((double)visu.width()) / ((double)visu.height()), 1.0f);
		//camera.translateLocal(Math::makeVector(2000.f, 3000.f, 700.f));
		scene.setCamera(camera);
	}
	createGround(scene);
}




//~4 000 000 de triangles
void initEngine(Geometry::Scene & scene, Visualizer::Visualizer const& visu)
{
	Geometry::BoundingBox box(Math::makeVector(0.0f, 0.0f, 0.0f), Math::makeVector(0.0f, 0.0f, 0.0f));
	Geometry::Loader3ds loader(m_modelDirectory + "\\engine\\engine.3ds", m_modelDirectory + "\\engine\\Texture");
	// We remove the specular components of the materials...
	::std::vector<Geometry::Material*> materials = loader.getMaterials();
	for (auto it = materials.begin(), end = materials.end(); it != end; ++it)
	{
		//(*it)->setSpecular(RGBColor());
		//((GeometryCollection::Phong*)(*it))->setAmbient(((GeometryCollection::Phong*)(*it))->getAmbient() * 0.01);
		//(*it)->setEmissive(RGBColor());
	}

	for (size_t cpt = 0; cpt < loader.getMeshes().size(); ++cpt)
	{
		scene.add(new Geometry::GeometryCollection(*loader.getMeshes()[cpt]));
		box.update(loader.getMeshes()[cpt]->box());
	}

	// 2.2 Adds point lights in the scene 
	Geometry::BoundingBox sb = scene.getBoundingBox();
	Math::Vector3f position = sb.max();
	Geometry::PointLight light1(position, RGBColor(1.0, 1.0, 1.0) * 1);
	scene.add(light1);
	position[1] = -position[1];
	Geometry::PointLight light2(position, RGBColor(1.0, 1.0, 1.0) * 1);
	scene.add(light2);
	{
		Geometry::Camera camera(Math::makeVector(30.f, 0.f, 10.0f), Math::makeVector(0.f, 0.f, 0.f), 0.3f, ((double)visu.width()) / ((double)visu.height()), 1.0f);
		//camera.translateLocal(Math::makeVector(2000.f, 3000.f, 700.f));
		scene.setCamera(camera);
	}
	createGround(scene);
}






enum RenderMode : int {
	rayTracing = 0, box = 1, normal = 2, uv = 3, materialID = 4, zBuffer = 5, AmbientOcclusion = 6,
	pathTracing = 7, iterativePathTracing = 8, lightTracing = 9, bdpt = 10, MISPathTracing = 11, naivePathTracing = 12,
	naiveBDPT = 13, OptiMISBDPT = 14, OptimalDirect = 15, PhotonMapper = 16, ProgressivePhotonMapper = 17, SimpleVCM = 18,
	UncorellatedBDPT = 19, VCM = 20, OptiVCM = 21,
};

std::vector<Integrator::Integrator*> init_integrators(unsigned int sample_per_pixel, unsigned int maxLen, double alpha, unsigned int lights_divisions, size_t w, size_t h)
{
	std::vector<Integrator::Integrator*> res(256);

	{
		res[RenderMode::rayTracing] = new Integrator::RayTracingIntegrator(sample_per_pixel, w, h);
		res[RenderMode::rayTracing]->setLen(maxLen);
	}

	{
		res[RenderMode::OptimalDirect] = new Integrator::OptimalDirect(sample_per_pixel, w, h);
		res[RenderMode::OptimalDirect]->setLen(maxLen);
	}


	{
		res[RenderMode::naivePathTracing] = new Integrator::NaivePathTracingIntegrator(sample_per_pixel, w, h);
		res[RenderMode::naivePathTracing]->setLen(maxLen);
		res[RenderMode::naivePathTracing]->m_alpha = alpha;
	}

	{
		res[RenderMode::iterativePathTracing] = new Integrator::IterativePathTracingIntegrator(sample_per_pixel, w, h);
		res[RenderMode::iterativePathTracing]->setLen(maxLen);
		res[RenderMode::iterativePathTracing]->m_alpha = alpha;
	}

	{
		res[RenderMode::MISPathTracing] = new Integrator::MISPathTracingIntegrator(sample_per_pixel, w, h);
		res[RenderMode::MISPathTracing]->setLen(maxLen);
		res[RenderMode::MISPathTracing]->m_alpha = alpha;
	}

	{
		res[RenderMode::lightTracing] = new Integrator::LightIntegrator(sample_per_pixel, w, h);
		res[RenderMode::lightTracing]->setLen(maxLen);
		res[RenderMode::lightTracing]->m_alpha = alpha;
	}

	{
		res[RenderMode::bdpt] = new Integrator::BidirectionalIntegrator(sample_per_pixel, w, h);
		res[RenderMode::bdpt]->setLen(maxLen);
		res[RenderMode::bdpt]->m_alpha = alpha;
	}

	{
		res[RenderMode::OptiMISBDPT] = new Integrator::OptiMISBDPT(sample_per_pixel, w, h);
		res[RenderMode::OptiMISBDPT]->setLen(maxLen);
		res[RenderMode::OptiMISBDPT]->m_alpha = alpha;
	}

	{
		res[RenderMode::UncorellatedBDPT] = new Integrator::UncorellatedBDPT(sample_per_pixel, w, h);
		res[RenderMode::UncorellatedBDPT]->setLen(maxLen);
		res[RenderMode::UncorellatedBDPT]->m_alpha = alpha;
	}

	{
		res[RenderMode::PhotonMapper] = new Integrator::PhotonMapper(sample_per_pixel, w, h);
		res[RenderMode::PhotonMapper]->setLen(maxLen);
	}

	{
		res[RenderMode::ProgressivePhotonMapper] = new Integrator::ProgressivePhotonMapper(sample_per_pixel, w, h);
		res[RenderMode::ProgressivePhotonMapper]->setLen(maxLen);
	}

	//{
	//	res[RenderMode::SimpleVCM] = new Integrator::SimpleVCM(sample_per_pixel, w, h);
	//	res[RenderMode::SimpleVCM]->setLen(maxLen);
	//}

	{
		res[RenderMode::VCM] = new Integrator::VCM(sample_per_pixel, w, h);
		res[RenderMode::VCM]->setLen(maxLen);
	}

	{
		res[RenderMode::OptiVCM] = new Integrator::OptiVCM(sample_per_pixel, w, h);
		res[RenderMode::OptiVCM]->setLen(maxLen);
	}

	//{
	//	res[RenderMode::naiveBDPT] = new Integrator::NaiveBidirectionalIntegrator(sample_per_pixel, w, h);
	//	res[RenderMode::naiveBDPT]->setLen(maxLen);
	//	res[RenderMode::naiveBDPT]->m_alpha = alpha;
	//}

	{
		res[RenderMode::box] = new Integrator::BOXIntegrator(sample_per_pixel, w, h);
	}

	{
		res[RenderMode::normal] = new Integrator::NormalIntegrator(sample_per_pixel, w, h);
	}

	{
		res[RenderMode::uv] = new Integrator::UVIntegrator(sample_per_pixel, w, h);
	}

	{
		res[RenderMode::materialID] = new Integrator::MIDIntegrator(sample_per_pixel, w, h);
	}

	{
		res[RenderMode::zBuffer] = new Integrator::ZIntegrator(sample_per_pixel, w, h);
	}

	{
		res[RenderMode::AmbientOcclusion] = new Integrator::AmbientOcclusionIntegrator(sample_per_pixel, w, h);
	}


	return res;
}





void waitKeyPressed()
{
	
	SDL_Event event;
	bool done = false;
	while (!done) {
		while (SDL_PollEvent(&event)) {
			switch (event.type) {
			case SDL_KEYDOWN:
				/*break;*/
			case SDL_QUIT:
				done = true;
				break;
			default:
				break;
			}
		}/*while*/
	}/*while(!done)*/
}




const std::string help_message =
"---REAL TIME RENDERING---\n"
"--Controls:\n"
"- Z Q S D Space LCtrl: Camera mouvement\n"
"- Arrows: Camera rotation\n"
"- LShift: Going faster\n"
"- CapsLock: Lock Going faster\n"
"- A: Launch multi pass mode (Stop with RETURN)\n"
"- I: Print infos (camera, visualizer and rendering technique)\n"
"- 0 - 9: set resolution scale\n"
"- F: Print FrameRate\n"
"- Rendering modes:\n"
"	- R: Ray tracing\n"
"	- T: Naive Path Tracing\n"
"	- Y: Iterative Path Tracing\n"
"	- G: MIS Path Tracing\n"
"	- L: Light Tracing\n"
"	- K: Bidirectional Path Tracing\n"
"	- B: Boxes of the acceleration structure\n"
"	- P: Depth\n"
"	- N: Normals\n"
"	- U: UV coordinates\n"
"   - V: Ambiant occlusion\n"
"\n"
"- O: save image (in PPM)\n"
"- M: Mute (TODO)\n"
"- H: print Help\n"
;


enum RenderOption { RealTime, Pass, Debug };


void get_input(std::vector<SDL_Event> const& events,bool * keys, RenderMode & render_mode, int & scale, RenderOption & rt)
{
	for(SDL_Event const& e: events)
	{
		if (e.type == SDL_KEYDOWN)
		{
			switch (e.key.keysym.sym)
			{
			case SDLK_UP:
				keys[0]=1;
				break;
			case SDLK_DOWN:
				keys[1] = 1;
				break;
			case SDLK_LEFT:
				keys[2] = 1;
				break;
			case SDLK_RIGHT:
				keys[3] = 1;
				break;
			case SDLK_z:
				keys[4] = 1;
				break;
			case SDLK_s:
				keys[5] = 1;
				break;
			case SDLK_d:
				keys[6] = 1;
				break;
			case SDLK_q:
				keys[7] = 1;
				break;
			case SDLK_SPACE:
				keys[8] = 1;
				break;
			case SDLK_LCTRL:
				keys[9] = 1;
				break;
			case SDLK_LSHIFT:
				keys[10] = 1;
				break;
			case SDLK_CAPSLOCK:
				keys[11] = ! keys[11];
				break;
			case SDLK_i:
				keys[12] = 1;
				break;
			case SDLK_f:
				keys[13] = !keys[13];
				break;
			case SDLK_r:
				if (render_mode <= RenderMode::AmbientOcclusion)
				{
					render_mode = static_cast<RenderMode>((render_mode + 1) % (RenderMode::AmbientOcclusion+1));
				}
				else
				{
					render_mode = RenderMode::rayTracing;
				}
				break;
			case SDLK_t:
				if (render_mode != RenderMode::naivePathTracing)
					std::cout << "switching to naive Path Tracing" << std::endl;
				render_mode = RenderMode::naivePathTracing;
				break;
			case SDLK_y:
				if (render_mode != RenderMode::iterativePathTracing)
					std::cout << "switching to iterative Path Tracing" << std::endl;
				render_mode = RenderMode::iterativePathTracing;
				break;
			case SDLK_g:
				if (render_mode != RenderMode::MISPathTracing)
					std::cout << "switching to MIS Path Tracing" << std::endl;
				render_mode = RenderMode::MISPathTracing;
				break;
			case SDLK_l:
				if (render_mode != RenderMode::lightTracing)
					std::cout << "switching to Light Tracing" << std::endl;
				render_mode = RenderMode::lightTracing;
				break;
			case SDLK_k:
				if (render_mode != RenderMode::bdpt)
					std::cout << "switching to BDPT" << std::endl;
				render_mode = RenderMode::bdpt;
				break;
			case SDLK_x:
				if (render_mode != RenderMode::OptimalDirect)
					std::cout << "switching to Optimis Direct" << std::endl;
				render_mode = RenderMode::OptimalDirect;
				break;
			case SDLK_w:
				if (render_mode != RenderMode::OptiMISBDPT)
					std::cout << "switching to Optimis BDPT" << std::endl;
				render_mode = RenderMode::OptiMISBDPT;
				break;
			case SDLK_m:
				if (render_mode != RenderMode::PhotonMapper)
					std::cout << "switching to Photon Mapping" << std::endl;
				render_mode = RenderMode::PhotonMapper;
				break;
			case SDLK_p:
				if (render_mode != RenderMode::ProgressivePhotonMapper)
					std::cout << "switching to Progressive Photon Mapping" << std::endl;
				render_mode = RenderMode::ProgressivePhotonMapper;
				break;
			case SDLK_v:
				if (render_mode != RenderMode::VCM)
					std::cout << "switching to VCM" << std::endl;
				render_mode = RenderMode::VCM;
				break;
			case SDLK_c:
				if (render_mode != RenderMode::OptiVCM)
					std::cout << "switching to Optimal VCM" << std::endl;
				render_mode = RenderMode::OptiVCM;
				break;
			case SDLK_u:
				if (render_mode != RenderMode::UncorellatedBDPT)
					std::cout << "switching to UBDPT" << std::endl;
				render_mode = RenderMode::UncorellatedBDPT;
				break;
			case SDLK_a:
				rt = RenderOption::Pass;
				break;
			case SDLK_h:
				std::cout << help_message << std::endl;
				break;
			case SDLK_KP_PLUS:
				keys[14] = 1;
				break;
			case SDLK_KP_MINUS:
				keys[15] = 1;
				break;
			default:
				break;
			}

			//NUMBERS
			if (e.key.keysym.sym >= SDLK_0 && e.key.keysym.sym <= SDLK_9)
			{
				scale = e.key.keysym.sym == SDLK_0 ? 10 : e.key.keysym.sym - SDLK_0;
			}
		}
		else if (e.type == SDL_KEYUP)
		{
			switch (e.key.keysym.sym)
			{
			case SDLK_UP:
				keys[0] = 0;
				break;
			case SDLK_DOWN:
				keys[1] = 0;
				break;
			case SDLK_LEFT:
				keys[2] = 0;
				break;
			case SDLK_RIGHT:
				keys[3] = 0;
				break;
			case SDLK_z:
				keys[4] = 0;
				break;
			case SDLK_s:
				keys[5] = 0;
				break;
			case SDLK_d:
				keys[6] = 0;
				break;
			case SDLK_q:
				keys[7] = 0;
				break;
			case SDLK_SPACE:
				keys[8] = 0;
				break;
			case SDLK_LCTRL:
				keys[9] = 0;
				break;
			case SDLK_LSHIFT:
				keys[10] = 0;
				break;
			case SDLK_KP_PLUS:
				keys[14] = 0;
				break;
			case SDLK_KP_MINUS:
				keys[15] = 0;
				break;
			default:
				break;
			}
		}
	}
}


template <typename floot>
void testPrecision()
{
	std::cout.precision(20);
	floot f = 1;
	for (unsigned int i = 0; i < 1024; ++i)
	{
		floot tmp = f + 1;
		std::cout << i << " : ";
		if (f != tmp)
		{
			std::cout << "Ok for " << f << std::endl;
		}
		else
		{
			std::cout<<"WRONG for " << f << std::endl;
			break;
		}
		f *= 2;
	}
}

void testWaveLength()
{
	Image::Image<RGBColor> img(290, 50);
	Parallel::ParallelFor(0, img.width(), [&img](const int i)
		{
			double l = i + 380;
			RGBColor c = RGBColor::fromWaveLength(l);
			for (int j = 0; j < img.height(); ++j)
			{
				img(i, j) = c;
			}
		});
	Image::ImWrite::write(img);

	img.fill(img.mean());
	Image::ImWrite::write(img);
}



int main(int argc, char** argv)
{
	Parallel::init();
	if (argc > 1)
	{
		return Auto::Auto::__main__(argc, argv);
	}

	int nthread = 4 * 2 * 2;
	Parallel::setNumThreads(nthread);


#ifdef _DEBUG
	int scale = 10;
#else
	int scale = 1;
#endif

	// 1 - Initializes a window for rendering
	//Visualizer::Visualizer visu(2000, 2000, scale);// pour les ecrans 4K
	Visualizer::Visualizer visu(1000, 1000, scale);
	//Visualizer::Visualizer visu(2000, 1000, scale);
	//Visualizer::Visualizer visu(1900, 1000, scale);
	//Visualizer::Visualizer visu(1000, 500, scale);
	//Visualizer::Visualizer visu(500, 500, scale);
	//Visualizer::Visualizer visu(300, 300, scale) ;
	//Visualizer::Visualizer visu(250, 250, scale) ;
	//Visualizer::Visualizer visu(200, 200, scale) ;
	//Visualizer::Visualizer visu(150, 150, scale) ;
	//Visualizer::Visualizer visu(100, 100, scale) ;
	//Visualizer::Visualizer visu(1, 1, 1);

	//testSqaureSample(visu);
	//testSDiskSample(visu);
	//testSpecular(visu);

	// 2 - Initializes the scene
	Geometry::Scene scene;

	// 2.1 initializes the geometry (choose only one initialization)
	//Auto::initRealCornell(scene, visu.width(), visu.height(), 0, 1, 0);
	//Auto::initCausticCornell(scene, visu.width(), visu.height(), 0, 1, 0);
	//Auto::initCausticCornell(scene, visu.width(), visu.height(), 1, 1, 0);
	//Auto::initCornellLamp(scene, visu.width(), visu.height());
	//Auto::initSimpleCornell(scene, visu.width(), visu.height(), 0);
	Auto::initVeach(scene, visu.width(), visu.height());
	//Auto::initVeach(scene, visu.width(), visu.height(), 5);
	//Auto::initSDSCornell(scene, visu.width(), visu.height());
	//Auto::initCornellLaserPrism(scene, visu.width(), visu.height());
	//Auto::initTest(scene, visu.width(), visu.height());
	//Auto::initTestNonSymmetry(scene, visu.width(), visu.height(), 0);
	
	//initGuitar(scene, visu);
	//initDog(scene, visu);
	//initGarage(scene, visu);
	//initRobot(scene, visu);
	//initTemple(scene, visu);
	//initGraveStone(scene, visu);
	//initBoat(scene, visu);
	//initSombrero(scene, visu);
	//initTibetHouse(scene, visu);
	//initTibetHouseInside(scene, visu);
	//initMedievalCity(scene, visu);
	//initEngine(scene, visu);
	

	std::cout << "Building the acceleration structure" << std::endl;
	tic();
	scene.preCompute(8, 1, 1);
	std::cout << "Done! ";
	toc();
	
	// Shows stats
	scene.printStats();


	// 3 - Computes the scene
	unsigned int sample_per_pixel = 16*16;
										
	// max lenght is included
	unsigned int maxLen = 3;

	unsigned int lights_divisions = sample_per_pixel;


	double alpha = 1;

	//scene.setBackgroundColor({ 0.065, 0.065, 0.088 });
	//Auto::setBackground(scene);

	//Compute the light samplers with only 1 division to get a better result in real time
	scene.compute_light_samplers(1);

	scene.check_capacity();


	RenderOption render_option = RenderOption::RealTime;
	RenderMode render_mode = RenderMode::rayTracing;

	std::vector<Integrator::Integrator*> integrators = init_integrators(sample_per_pixel, maxLen, alpha, lights_divisions, visu.width(), visu.height());

	Integrator::Integrator* integrator = integrators[render_mode];


	((Integrator::PhotonMapper*)integrators[RenderMode::PhotonMapper])->setParams(scene, 0.01, 1000000);
	((Integrator::ProgressivePhotonMapper*)integrators[RenderMode::ProgressivePhotonMapper])->setParams(scene, 0.01, 1000000);
	//((Integrator::SimpleVCM*)integrators[RenderMode::SimpleVCM])->setParams(scene, 0.01, 1000);
	((Integrator::VCM*)integrators[RenderMode::VCM])->setParams(scene, 0.0075, 0.25);
	((Integrator::OptiVCM*)integrators[RenderMode::OptiVCM])->setParams(scene, 0.0075, 0.25);

	std::cout << help_message << std::endl;


	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
	Geometry::Camera & cam = scene.m_camera;

	Math::Vector3f pos = cam.getPosition();
	Math::Vector3f dir = cam.m_front;

	//UP DOWN LEFT RIGHT Z S D Q SPACE LCTRL LSHIFT CAPSLOCK C F + -
	const size_t num_keys = 16;
	bool keys[num_keys];
	for (int i = 0; i < num_keys; ++i)
	{
		keys[i] = 0;
	}
		


	const double speed = 2;
	const double angle_speed = 2;

	double forward(0), rightward(0), upward(0), inclination, azimuth;

	Math::Vector3f speed_vec, cam_speed_vec;

	Math::Vector3f dir_sphere = spherical_coordinate(dir);

	inclination = dir_sphere[1];
	azimuth = dir_sphere[2];

	cam.resolution = visu.width() * visu.height();

		
	std::chrono::high_resolution_clock::time_point t2;
	unsigned int frame_count=0;
	for ever
	{
		t2 = std::chrono::high_resolution_clock::now();
		get_input(visu.events, keys, render_mode, scale,render_option);

		integrator = integrators[render_mode];

		visu.update_scale(scale);

		if (render_option == RenderOption::RealTime)
		{

			std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);

			double dt = time_span.count();

			t1 = t2;


			if (keys[0])
			{
				inclination -= angle_speed * dt;
			}
			if (keys[1])
			{
				inclination += angle_speed * dt;
			}
			if (keys[2])
			{
				azimuth += angle_speed * dt;
			}
			if (keys[3])
			{
				azimuth -= angle_speed * dt;
			}
			if (keys[4])
			{
				forward += speed * dt;
			}
			if (keys[5])
			{
				forward -= speed * dt;
			}
			if (keys[6])
			{
				rightward += speed * dt;
			}
			if (keys[7])
			{
				rightward -= speed * dt;
			}
			if (keys[8])
			{
				upward += speed * dt;
			}
			if (keys[9])
			{
				upward -= speed * dt;
			}

			if (inclination > Math::pi)
			{

				inclination = Math::pi - 0.00000001;
			}
			else if (inclination < 0)
			{

				inclination = 0.00000001;
			}

			double sin_inc = sin(inclination);
			double cos_inc = cos(inclination);
			if (inclination > Math::pi)
			{
				cos_inc = -cos_inc;
			}
			dir = Math::makeVector(sin_inc * cos(azimuth - Math::pi), sin_inc * sin(azimuth - Math::pi), cos_inc);

			double caps_speed = 1;
			speed_vec = cam.m_right*rightward + cam.m_front*forward - cam.m_down*upward;
			if (keys[10])
			{
				caps_speed *= 15;
			}
			if (keys[11])
			{
				caps_speed *= 15;
			}
			speed_vec *= caps_speed;


			pos = pos + speed_vec;
			dir = dir + speed_vec;

			const double cam_zoom_speed = 0.01 * caps_speed;
			if (keys[14])
			{
				scene.m_camera.scale_screen(1.0 - cam_zoom_speed);
			}
			if (keys[15])
			{
				scene.m_camera.scale_screen(1.0 + cam_zoom_speed);
			}

			cam.update_both(pos + speed_vec, pos + dir);

			cam.resolution = visu.width() * visu.height();

#ifdef COUNT_RAYS
			tbb::tick_count tic = tbb::tick_count::now();
#endif

			integrator->fastRender(scene, visu);
			//scene.realtime_compute(render_mode);




#ifdef COUNT_RAYS
			tbb::tick_count toc = tbb::tick_count::now();
			std::cout << "rays / sec: " << scene.ray_counter / (toc - tic).seconds() << std::endl;
			scene.ray_counter = 0;
#endif

			forward = 0;
			rightward = 0;
			upward = 0;

			//framerate

			if (keys[13] && !(frame_count % 50))
			{
				std::cout << "fps: " << 1.0 / dt << std::endl;
				//std::cout << "mean: " << mean_framerate << std::endl;
			}
			++frame_count;

			if (keys[12])
			{
				std::cout << "Camera position: " << cam.m_position << "\t Camera target: " << cam.m_target << std::endl;
				std::cout << "focal width: " << cam.m_planeDistance << std::endl;
				visu.print_info(std::cout);
				keys[12] = 0;
				
			}
		}
		else if(render_option == RenderOption::Pass)
		{
			std::cout<<
			"---All pass Mode---\n"
			"Press RETURN to stop\n"
			<< std::endl;

			scene.compute_light_samplers(lights_divisions);

			integrator->render(scene, visu);
			//scene.compute(render_mode, subPixelSampling, passPerPixel);
			
			scene.compute_light_samplers(1);

			render_option = RenderOption::RealTime;
		}
		else if (render_option == RenderOption::Debug)
		{
			visu.clean();
			visu.update();
			
			integrator->debug(scene, visu);
			//scene.debug(x, y, visu.width(), visu.height(), render_mode, subPixelSampling, passPerPixel);
			
		}
	}


	// 4 - waits until a key is pressed
	waitKeyPressed();

	//clean the integrators
	{
		for (Integrator::Integrator* integ : integrators)
		{
			if (integ != nullptr)
			{
				delete integ;
				integ = nullptr;
			}
		}
	}
	
	return 0;
}
