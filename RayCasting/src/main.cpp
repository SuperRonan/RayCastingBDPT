#include <settings.h>

#include <Geometry/Texture.h>
#include <Math/Vectorf.h>
#include <stdlib.h>
#include <iostream>
#include <Geometry/RGBColor.h>
#include <Geometry/Material.h>
#include <Geometry/PointLight.h>
#include <Geometry/Camera.h>
#include <Geometry/Cube.h>
#include <Geometry/Disk.h>
#include <Geometry/Cylinder.h>
#include <Geometry/Cone.h>
#include <Visualizer/Visualizer.h>
#include <Geometry/Scene.h>
#include <Geometry/Cornel.h>
#include <Geometry/Loader3ds.h>
#include <Geometry/BoundingBox.h>
#include <omp.h>
#include <chrono>
#include <array>
#include <tbb/tick_count.h>
#include "Geometry/Phong.h"
#include "Geometry/CartoonMaterial.h"
#include <Geometry/Lambert.h>
#include <Geometry/Specular.h>
#include <Geometry/DeltaMirror.h>
#include <Math/Sampler.h>

#include <Integrators/DirectIntegrator.h>
#include <Integrators/RegularIntegrators.h>
#include <Integrators/RayTracingIntegrator.h>
#include <Integrators/Integrator.h>
#include <Integrators/PathTracingIntegrator.h>
#include <Integrators/ZIntegrator.h>
#include <Integrators/LightIntegrator.h>
#include <Integrators/MISPathTracingIntegrator.h>
#include <Integrators/BidirectionalIntegrator.h>

#include <Auto/Auto.h>
#include <Auto/TestScenes.h>

#define ever (;;)

std::ostream& __clk_out = std::cout;
std::chrono::high_resolution_clock __clk;
std::stack<std::chrono::time_point<std::chrono::high_resolution_clock>> __tics;


void tic()
{
	__tics.push(__clk.now());
}

void toc()
{
	std::chrono::time_point<std::chrono::high_resolution_clock> __toc = __clk.now(), __tic = __tics.top();
	__tics.pop();
	std::chrono::duration<double>  __duration = std::chrono::duration_cast<std::chrono::duration<double>>(__toc - __tic);
	__clk_out << __duration.count() << "s" << std::endl;
}

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
/// \author	F. Lamarche, Universit� de Rennes 1
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
	Geometry::Material* material = new Geometry::Specular(0.5, 1000);
	//GeometryCollection::Material* material = new GeometryCollection::DeltaMirror(0.5);
	//GeometryCollection::Material* material = new GeometryCollection::Lambertian(0.7);

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
/// \author	F. Lamarche, Universit� de Rennes 1
/// \date	03/12/2013
///
/// \param [in,out]	scene	The scene.
////////////////////////////////////////////////////////////////////////////////////////////////////
void createSurfaceLight(Geometry::Scene & scene, double value)
{
	Geometry::BoundingBox sb = scene.getBoundingBox();
	Geometry::Material * material = new Geometry::Phong(RGBColor(), RGBColor(), 100.0f, RGBColor(value, value, value));
	
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
/// \author	F. Lamarche, Universit� de Rennes 1
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
/// \author	F. Lamarche, Universit� de Rennes 1
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
/// \author	F. Lamarche, Universit� de Rennes 1
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

	for (size_t cpt = 0; cpt < loader.getMeshes().size(); ++cpt)
	{
		scene.add(new Geometry::GeometryCollection(*loader.getMeshes()[cpt]));
	}

	// 2.2 Adds point lights in the scene 
	Geometry::BoundingBox sb = scene.getBoundingBox();
	Math::Vector3f position = sb.max();
	Geometry::PointLight light1(position, RGBColor(1.0, 1.0, 1.0) *.5);
	scene.add(light1);
	position[1] = -position[1];
	Geometry::PointLight light2(position, RGBColor(1.0, 1.0, 1.0) *.5);


	Geometry::Material * emisive = new Geometry::Phong(0, 0, 0, 0, { 2, 1, 1 });
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
			((*it))->setEmissive(RGBColor(1.0, 1.0, 1.0)*10*3);
			
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






enum RenderMode { rayTracing = 0, pathTracing = 1, iterativePathTracing = 2, lightTracing = 3, box = 4, normal = 5, uv = 6, materialID = 7, zBuffer = 8, 
	bdpt = 9, MISPathTracing=10, naivePathTracing=11, AmbientOcclusion=12, naiveBDPT=13, OptiMISBDPT=14 };

std::vector<Integrator::Integrator*> init_integrators(unsigned int sample_per_pixel, unsigned int maxBounce, double alpha, unsigned int lights_divisions, size_t w, size_t h)
{
	std::vector<Integrator::Integrator*> res(256);

	{
		res[RenderMode::rayTracing] = new Integrator::RayTracingIntegrator(sample_per_pixel, w, h);
		res[RenderMode::rayTracing]->setDepth(maxBounce);
	}

	{
		res[RenderMode::pathTracing] = new Integrator::RecursivePathTracingIntegrator(sample_per_pixel, w, h);
		res[RenderMode::pathTracing]->setDepth(maxBounce);
		res[RenderMode::pathTracing]->m_alpha = alpha;
	}

	{
		res[RenderMode::naivePathTracing] = new Integrator::NaivePathTracingIntegrator(sample_per_pixel, w, h);
		res[RenderMode::naivePathTracing]->setDepth(maxBounce);
		res[RenderMode::naivePathTracing]->m_alpha = alpha;
	}

	{
		res[RenderMode::iterativePathTracing] = new Integrator::IterativePathTracingIntegrator(sample_per_pixel, w, h);
		res[RenderMode::iterativePathTracing]->setDepth(maxBounce);
		res[RenderMode::iterativePathTracing]->m_alpha = alpha;
	}

	{
		res[RenderMode::MISPathTracing] = new Integrator::MISPathTracingIntegrator(sample_per_pixel, w, h);
		res[RenderMode::MISPathTracing]->setDepth(maxBounce);
		res[RenderMode::MISPathTracing]->m_alpha = alpha;
	}

	{
		res[RenderMode::lightTracing] = new Integrator::LightIntegrator(sample_per_pixel, w, h);
		res[RenderMode::lightTracing]->setDepth(maxBounce);
		res[RenderMode::lightTracing]->m_alpha = alpha;
	}

	{
		res[RenderMode::bdpt] = new Integrator::BidirectionalIntegrator(sample_per_pixel, w, h);
		res[RenderMode::bdpt]->setDepth(maxBounce);
		res[RenderMode::bdpt]->m_alpha = alpha;
	}

	//{
	//	res[RenderMode::naiveBDPT] = new Integrator::NaiveBidirectionalIntegrator(sample_per_pixel, w, h);
	//	res[RenderMode::naiveBDPT]->setDepth(maxBounce);
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









////////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn	void waitKeyPressed()
///
/// \brief	Waits until a key is pressed.
/// 		
/// \author	F. Lamarche, Universit� de Rennes 1
/// \date	03/12/2013
////////////////////////////////////////////////////////////////////////////////////////////////////
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
			case SDLK_b:
				if (render_mode != RenderMode::box)
					std::cout << "switching to box rendering" << std::endl;
				render_mode = RenderMode::box;
				break;
			case SDLK_r:
				if (render_mode != RenderMode::rayTracing)
					std::cout << "switching to Ray Tracing" << std::endl;
				render_mode = RenderMode::rayTracing;
				break;
			case SDLK_u:
				if (render_mode != RenderMode::uv)
					std::cout << "switching to uv rendering" << std::endl;
				render_mode = RenderMode::uv;
				break;
			case SDLK_n:
				if (render_mode != RenderMode::normal)
					std::cout << "switching to normal rendering" << std::endl;
				render_mode = RenderMode::normal;
				break;
			case SDLK_m:
				if (render_mode != RenderMode::materialID)
					std::cout << "switching to materialID rendering" << std::endl;
				render_mode = RenderMode::materialID;
				break;
			case SDLK_p:
				if (render_mode != RenderMode::zBuffer)
					std::cout << "switching to Z buffer rendering" << std::endl;
				render_mode = RenderMode::zBuffer;
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
			case SDLK_j:
				if (render_mode != RenderMode::naiveBDPT)
					std::cout << "switching to Naive BDPT" << std::endl;
				render_mode = RenderMode::naiveBDPT;
				break;
			case SDLK_v:
				if (render_mode != RenderMode::AmbientOcclusion)
					std::cout << "switching to Ambient Occlusion" << std::endl;
				render_mode = RenderMode::AmbientOcclusion;
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



//void testSquareSample(Visualizer::Visualizer& visu)
//{
//	GeometryCollection::Material mat(1);
//	GeometryCollection::Square sqare(&mat, { 0.5, 0.5, 0 });
//	sqare.divide(100);
//	sqare.build_lights();
//	Image::Image<RGBColor> img(visu.width(), visu.height());
//	img.fill();
//	int s = GeometryCollection::LightSampleStack::capacity;
//	GeometryCollection::Camera cam({ 0, 0, 1 }, { 0, 0, 0 });
//	cam.m_right = { 1, 0, 0 };
//	cam.m_down = { 0, 1, 0 };
//	Math::Sampler sampler(1ull);
//	for (size_t i = 1; i < 100; ++i)
//	{
//		GeometryCollection::LightSampleStack lss;
//		sqare.sampleLights(lss, sampler, s);
//		for (GeometryCollection::SurfaceLightSample const& sls : lss)
//		{
//			Math::Vector2f camuv = cam.screen_position(sls.vector) - 0.5;
//			Math::Vector2f directuv = { sls.vector[0], sls.vector[1] };
//			/*std::cout << "-----------------" << std::endl;
//			std::cout << camuv << std::endl;
//			std::cout << directuv << std::endl;*/
//			Math::Vector2f uv = camuv;
//			int x = uv[0] * visu.width();
//			int y = uv[1] * visu.height();
//			img[x][y] += 1;
//			visu.plot(x, y, img[x][y]);
//			
//		}
//		visu.update();
//		sqare.increment_offset(s);
//	}
//	while (visu.update() != Visualizer::Visualizer::done)
//	{
//
//	}
//}
//
//
//
//void testSpecular(Visualizer::Visualizer & visu)
//{
//	GeometryCollection::Specular spec(1.0, 10.0);
//
//	visu.clean();
//	Math::Sampler sampler;
//
//	GeometryCollection::Hit hit;
//	hit.primitive_normal = { 0, 0, 1 };
//	hit.reflected = Math::Vector3f(0, 1, 1).normalized();
//
//	for (int sample = 0; sample < 1000; ++sample)
//	{
//		GeometryCollection::DirectionSample dir;
//		spec.sampleBSDF(hit, 1, 1, dir, sampler);
//
//		Math::Vector3f wi = dir.direction;
//
//		RGBColor bsdf = dir.bsdf;
//
//		std::cout << bsdf - spec.BSDF(hit, wi) << std::endl;
//		std::cout << dir.pdf - spec.pdf(hit, wi)<<std::endl;
//
//		
//	}
//
//	exit(0);
//}


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





int main(int argc, char ** argv)
{
	if (argc > 1)
	{
		return Auto::Auto::__main__(argc, argv);
	}



	int nthread = 4*2*2;
	omp_set_num_threads(nthread);

#ifdef _DEBUG
	int scale = 10;
#else
	int scale = 1;
#endif

	// 1 - Initializes a window for rendering
	//Visualizer::Visualizer visu(2000, 2000, scale);// pour les ecrans 4K
	Visualizer::Visualizer visu(1000, 1000, scale);
	//Visualizer::Visualizer visu(1900, 1000, scale);
	//Visualizer::Visualizer visu(500, 500, scale);
	//Visualizer::Visualizer visu(300, 300, scale) ;
	//Visualizer::Visualizer visu(250, 250, scale) ;
	//Visualizer::Visualizer visu(200, 200, scale) ;
	//Visualizer::Visualizer visu(150, 150, scale) ;
	//Visualizer::Visualizer visu(100, 100, scale) ;

	//testSqaureSample(visu);
	//testSDiskSample(visu);
	//testSpecular(visu);

	// 2 - Initializes the scene
	Geometry::Scene scene;

	// 2.1 initializes the geometry (choose only one initialization)
	Auto::initRealCornell(scene, visu.width(), visu.height(), 1, 1, 0);
	//Auto::initCornellLamp(scene, visu.width(), visu.height());
	//Auto::initSimpleCornell(scene, visu.width(), visu.height(), 2);
	//Auto::initVeach(scene, visu.width(), visu.height());
	//Auto::initTest(scene, visu.width(), visu.height());
	
	//initDiffuse(scene, visu);
	//initDiffuse_dif(scene, visu);
	//initDiffuse_surface_light<GeometryCollection::Phong>(scene, visu);
	//initDiffuseSpecular(scene, visu) ;//custom
	//initSpecular(scene, visu) ;
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
	
	


	/*
	std::cout << "Pre computing the shadows" << std::endl;
	scene.pre_compute_shadows(0.1, 0.1, 10, 10);
	std::cout << "Done!" << std::endl;
	std::cout << scene.num_shadow << std::endl;
	//*/
	
	// Shows stats
	scene.printStats();

	// 3 - Computes the scene
	unsigned int sample_per_pixel = 1024/4;
										
	unsigned int maxBounce = 10;			// Maximum number of bounces

	unsigned int lights_divisions = 16;



	double alpha = 0.9;

	//scene.setBackgroundColor({ 0.065, 0.065, 0.088 });
	//Auto::setBackground(scene);

	//Compute the light samplers with only 1 division to get a better result in real time
	scene.compute_light_samplers(1);

	scene.check_capacity();

	RenderOption render_option = RenderOption::RealTime;
	RenderMode render_mode = RenderMode::rayTracing;

	std::vector<Integrator::Integrator*> integrators = init_integrators(sample_per_pixel, maxBounce, alpha, lights_divisions, visu.width(), visu.height());

	Integrator::Integrator* integrator = integrators[render_mode];


	std::cout << help_message << std::endl;


	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
	Geometry::Camera & cam = scene.m_camera;

	Math::Vector3f pos = cam.getPosition();
	Math::Vector3f dir = cam.getRay(0.5, 0.5).direction().normalized();

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
