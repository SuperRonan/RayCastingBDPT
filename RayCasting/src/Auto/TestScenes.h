#pragma once

#include <Geometry/Scene.h>
#include <Geometry/Cube.h>
#include <Geometry/Disk.h>
#include <Geometry/DeltaMirror.h>
#include <Geometry/Cornel.h>
#include <Geometry/Cylinder.h>
#include <Image/ImWrite.h>
#include <Geometry/glass.h>
#include <Image/ImRead.h>

namespace Auto
{
	using RGBColor = Geometry::RGBColor;




	void setBackground(Geometry::Scene& scene, size_t w = 100, size_t h = 100)
	{
		Image::Image<RGBColor> img(w, h);
		Math::Vector2f center = { 0.5, 0.3 };
		OMP_PARALLEL_FOR
			for (long i = 0; i < w; ++i)
			{
				for (size_t j = 0; j < h; ++j)
				{

					double u = double(i) / double(w);
					double v = double(j) / double(h);
					Math::Vector2f uv = { u, v };
					double d = (center - uv).norm2();
					//double ed = 50*exp(-d / 0.001);
					double r = 0;
					double g = 0;
					double b = 0;
					if (10 * i % w == 0 || 10 * j % h == 0)
					{
						r = 10 * (u);
						b = 10 * (v);
					}

					RGBColor& pixel = img(i, j);
					pixel = RGBColor(r, g, b);
				}
			}

		Image::ImWrite::write(img, "skybox.exr");
		scene.setBackgroundPixels(img);
	}

	void filterGalaxy(Image::Image<RGBColor>& image, int iterations)
	{
		const auto hasEnoughNeighboor = [&](int x, int y)
		{
			int sum = 0;
			if (image.inBounds<int>({ x - 1, y }) && !image(x - 1, y).isBlack())
				++sum;
			if (image.inBounds<int>({ x + 1, y }) && !image(x + 1, y).isBlack())
				++sum;
			if (image.inBounds<int>({ x, y - 1 }) && !image(x, y - 1).isBlack())
				++sum;
			if (image.inBounds<int>({ x, y + 1 }) && !image(x, y + 1).isBlack())
				++sum;
			return sum > 1;
		};
		for (int iter = 0; iter < iterations; ++iter)
		{
			OMP_PARALLEL_FOR
				for (int i = 0; i < image.width(); ++i)
				{
					for (int j = 0; j < image.height(); ++j)
					{
						if (!image(i, j).isBlack() && !hasEnoughNeighboor(i, j))
						{
							image(i, j) = 0;
						}
					}
				}
		}
	}

	void setGalaxyBackground(Geometry::Scene& scene)
	{
		Image::Image<RGBColor> img;
		Image::ImRead::readEXR("../../galaxy.exr", img);
		OMP_PARALLEL_FOR
			for (int i = 0; i < img.width(); ++i)
			{
				for (int j = 0; j < img.height(); ++j)
				{
					
				}
			}
		//filterGalaxy(img, 2);
		//Image::ImWrite::write(img, "skybox.exr");
		scene.setBackgroundPixels(img);
		scene.setBackgroundColor(1);
	}




	void initBlackHole(Geometry::Scene& scene, size_t width, size_t height)
	{
		{
			Geometry::Camera camera(Math::Vector3f(2.f, 0.0f, 0.1), Math::makeVector(0.0f, 0.0f, 0.0f), 0.3f, ((double)width) / ((double)height), 1.0f);
			scene.setCamera(camera);
		}
		setGalaxyBackground(scene);
	}



	void initTest(Geometry::Scene& scene, size_t width, size_t height)
	{
		if (true)
		{
			Geometry::Material* white_diffuse = new Geometry::Lambertian<Geometry::REFLECT>(0.8);
			Geometry::Material* purple_diffuse = new Geometry::Lambertian<Geometry::REFLECT>({ 0.6, 0.1, 0.4 });
			Geometry::Material* green_diffuse = new Geometry::Lambertian<Geometry::REFLECT>({ 0.2, 0.7, 0.3 });
			Geometry::Material* white_emisive = new Geometry::Material({ 15, 11, 9 }, "");

			{
				Geometry::Square* square = new Geometry::Square(white_diffuse);
				square->scale(10);
				square->translate({ 0, 0, 0 });

				scene.add(square);
			}
			{
				Geometry::Square* square = new Geometry::Square(white_diffuse);
				square->scale(1);
				square->translate({ 0, 0, 0.1 });

				//scene.add(square);
			}
			if constexpr (false)//square light
			{
				Geometry::Square* square = new Geometry::Square(white_emisive);
				square->scale(1.5);
				square->translate({ 0, 0, 3.5 });

				scene.add(square);
			}
			else//sphere light
			{
				Geometry::Sphere sphere = Geometry::Sphere({ -1, 0, 3.5 }, 0.4, white_emisive);

				//scene.add(sphere);
			}

			{
				Geometry::OneTriangle* tri = new Geometry::OneTriangle(purple_diffuse);
				tri->rotate(Math::Quaternion<double>({ 0, 1, 0 }, -0.5));
				tri->translate({ 0, 0, 1.2 });


				scene.add(tri);
			}

			{
				Geometry::Cube* cube = new Geometry::Cube(green_diffuse);

				cube->scale(0.4);
				cube->rotate(Math::Quaternion<double>({ 1, 1, 1 }, 1));
				cube->translate({ -1, -1, 1.2 });

				scene.add(cube);
			}


			// 2.2 Adds point lights in the scene 
			{
				Geometry::PointLight pointLight(Math::makeVector(5.0f, 0.f, 2.0f), RGBColor(1.f, 1.f, 1.f));
				scene.add(pointLight);
			}
			{
				Geometry::Camera camera(Math::makeVector(2.f, 0.0f, 2.f), Math::makeVector(0.0f, 0.0f, 0.0f), 0.3f, ((double)width) / ((double)height), 1.0f);
				scene.setCamera(camera);
			}

			//scene.setBackgroundColor({ 0.065, 0.065, 0.088 });
			//setBackground(scene);
			setGalaxyBackground(scene);
		}
		else
		{

		}
	}










	void initRealCornell(Geometry::Scene& scene, size_t width, size_t height, int mode, bool colors, bool cylinder)
	{
		Geometry::Material* white = new Geometry::Lambertian<Geometry::REFLECT>(0.7);
		Geometry::Material* black = new Geometry::Lambertian<Geometry::REFLECT>(0);
		Geometry::Material* red = new Geometry::Lambertian<Geometry::REFLECT>({ 0.62, 0.061, 0.061 });
		Geometry::Material* green = new Geometry::Lambertian<Geometry::REFLECT>({ 0.122, 0.406, 0.1 });

		Geometry::Material* blue = new Geometry::Lambertian<Geometry::REFLECT>({ 0.1, 0.2, 0.75 });
		Geometry::Material* orange = new Geometry::Lambertian<Geometry::REFLECT>({ 0.8, 0.4, 0.1 });

		Geometry::Material* spec = new Geometry::Specular(1, 1000);
		Geometry::Material* mirror = new Geometry::DeltaMirror(1.0);

		double scale = 5;
		Geometry::Cornel::init_cornell(scene, white, white, white, nullptr, red, green, scale);

		double light_size = 1;
		Geometry::Material* light = new Geometry::Lambertian<Geometry::REFLECT>(0.78, RGBColor(16, 10, 5) / (light_size * light_size));//0.78, RGBColor(16, 10, 5)
		Geometry::Square* up_light = new Geometry::Square(light);
		up_light->scale(scale * 0.2 * light_size);
		up_light->translate({ 0, 0, scale / 2 - 0.001 });
		scene.add(up_light);


		//tall block
		{
			Geometry::GeometryCollection* cube;
			Geometry::Material* cube_mat;

			if (mode == 0)
			{
				if (colors)
					cube_mat = blue;
				else
					cube_mat = white;
			}
			else if (mode == 1)
				cube_mat = spec;
			else
				cube_mat = mirror;


			if (cylinder)//cylinder
			{
				//GeometryCollection::Cube* cube = 
				cube = new Geometry::CylinderFabrice(10, 1, 1, cube_mat, false, false);
				cube->scaleX(0.15 * scale);
				cube->scaleY(0.15 * scale);

			}
			else //cube
			{
				cube = new Geometry::Cube(cube_mat);
				cube->scaleX(0.3 * scale);
				cube->scaleY(0.3 * scale);
			}
			cube->scaleZ(0.6 * scale);
			cube->rotate(Math::Quaternion<double>({ 0, 0, 1 }, 0.3));
			cube->translate({ scale * 0.2, scale * 0.15, -scale * 0.21 });
			scene.add(cube);
		}

		//little block
		{
			Geometry::Cube* cube = new Geometry::Cube(colors ? orange : white);
			cube->scaleX(0.3 * scale);
			cube->scaleY(0.3 * scale);
			cube->scaleZ(0.3 * scale);
			cube->rotate(Math::Quaternion<double>({ 0, 0, 1 }, -0.3));
			cube->translate({ -scale * 0.1, -scale * 0.2, -scale * 0.351 });


			scene.add(cube);
		}


		// Sets the camera
		{
			Geometry::Camera camera({ -scale * 2.4, 0 ,0 }, { 0, 0, 0 }, 2.f, ((double)width) / ((double)height), 1.0f);
			scene.setCamera(camera);
		}
	}


	void initComplexCausticCornell(Geometry::Scene& scene, size_t width, size_t height)
	{
		Geometry::Material* white = new Geometry::Lambertian<Geometry::REFLECT>(0.7);
		Geometry::Material* black = new Geometry::Lambertian<Geometry::REFLECT>(0);
		Geometry::Material* red = new Geometry::Lambertian<Geometry::REFLECT>({ 0.62, 0.061, 0.061 });
		Geometry::Material* green = new Geometry::Lambertian<Geometry::REFLECT>({ 0.122, 0.406, 0.1 });

		Geometry::Material* blue = new Geometry::Lambertian<Geometry::REFLECT>({ 0.1, 0.2, 0.75 });
		Geometry::Material* orange = new Geometry::Lambertian<Geometry::REFLECT>({ 0.8, 0.4, 0.1 });

		Geometry::Material* spec = new Geometry::Specular(1, 1000);
		Geometry::Material* mirror = new Geometry::DeltaMirror(1.0);
		Geometry::Material* glass = new Geometry::Glass({ 1, 1, 1.1 }, 1.3);
		Geometry::Material* blurry = new Geometry::Lambertian<Geometry::LAMBERT_MODE::TRANSMIT>(0.9);

		double light_size = 1;
		Geometry::Material* light = new Geometry::Lambertian<Geometry::REFLECT>(0.78, RGBColor(16, 10, 5) / (light_size * light_size));//0.78, RGBColor(16, 10, 5)

		double scale = 5;

		double ratio = (double)width / (double)height;
		
		for(double y : {0.0})
		{
			Geometry::Square* up_light = new Geometry::Square(light);
			up_light->scale(scale * 0.2 * light_size);
			up_light->scaleY(ratio);
			up_light->translate({ 0, y*scale, scale / 2 - 0.001 });
			scene.add(up_light);
		}

		{
			Geometry::Square* ground = new Geometry::Square(white, Math::Vector3f(0, 0, -0.5) * scale, Math::Vector3f(scale, 0, 0), Math::Vector3f(0, ratio * scale, 0));
			scene.add(ground);

			Geometry::Square* ceiling = new Geometry::Square(white, Math::Vector3f(0, 0, 0.5) * scale, Math::Vector3f(scale, 0, 0), Math::Vector3f(0, ratio * scale, 0));
			scene.add(ceiling);

			Geometry::Square* left = new Geometry::Square(red, Math::Vector3f(0, ratio / 2.0 * scale, 0), Math::Vector3f(scale, 0, 0), Math::Vector3f(0, 0, scale));
			scene.add(left);

			Geometry::Square* right = new Geometry::Square(green, Math::Vector3f(0, -ratio / 2.0 * scale, 0), Math::Vector3f(scale, 0, 0), Math::Vector3f(0, 0, scale));
			scene.add(right);

			Geometry::Square* front = new Geometry::Square(white, Math::Vector3f(0.5 * scale, 0, 0), Math::Vector3f(0, ratio * scale, 0), Math::Vector3f(0, 0, scale));
			scene.add(front);
		}
		
		Geometry::Material* cube_mats[] = { glass, blurry};
		double cube_size = 0.3;
		Math::Vector3f centers[] = { Math::Vector3f(0.0, 0.2, -0.4 + cube_size/2), Math::Vector3f(0.0, -0.2, -0.4 + cube_size / 2) };
		double radii[] = { 0.4, 0.49 };
		for (int i = 0; i < 2; ++i)
		{
			Geometry::Cube* cube = new Geometry::Cube(cube_mats[i], centers[i] * scale, Math::Vector3f(scale * cube_size, 0, 0), Math::Vector3f(0, scale * cube_size, 0), Math::Vector3f(0, 0, scale * cube_size));
			scene.add(cube);

			Geometry::Sphere sphere = Geometry::Sphere(centers[i] * scale, cube_size * radii[i] * scale, orange);
			scene.add(sphere);
		}

		Geometry::Camera camera({ -scale * 2.4, 0 ,0 }, { 0, 0, 0 }, 2.f, ((double)width) / ((double)height), 1.0f);
		scene.setCamera(camera);
	}




	void initCausticCornell(Geometry::Scene& scene, size_t width, size_t height, int mode, bool colors, bool cylinder)
	{
		double scale = 5;
		double light_size = 1;
		Geometry::Material* white = new Geometry::Lambertian<Geometry::REFLECT>(0.7);
		Geometry::Material* black = new Geometry::Lambertian<Geometry::REFLECT>(0);
		Geometry::Material* red = new Geometry::Lambertian<Geometry::REFLECT>({ 0.62, 0.061, 0.061 });
		Geometry::Material* green = new Geometry::Lambertian<Geometry::REFLECT>({ 0.122, 0.406, 0.1 });
		Geometry::Material* light = new Geometry::Lambertian<Geometry::REFLECT>(0.78, RGBColor(16, 10, 5) / (light_size * light_size));//0.78, RGBColor(16, 10, 5)
		Geometry::Material* blue = new Geometry::Lambertian<Geometry::REFLECT>({ 0.1, 0.2, 0.75 });
		Geometry::Material* orange = new Geometry::Lambertian<Geometry::REFLECT>({ 0.8, 0.4, 0.1 });

		Geometry::Material* spec = new Geometry::Specular(1, 1000);
		Geometry::Material* mirror = new Geometry::DeltaMirror(1.0);
		Geometry::Material* glass = new Geometry::Glass({ 1, 1, 1.1}, mode ? 1.3 : 1.3);
		Geometry::Material* glass2 = new Geometry::Glass({ 1, 1, 1.1}, 1);


		if(mode)
			Geometry::Cornel::init_cornell(scene, white, white, mirror, white, red, green, scale);
		else
			Geometry::Cornel::init_cornell(scene, white, white, white, white, red, green, scale);
		
		
		Geometry::Square* up_light = new Geometry::Square(light);
		up_light->scale(scale * 0.2 * light_size);
		up_light->translate({ 0, 0, scale / 2 - 0.001 });
		scene.add(up_light);


		if (true)
		{
			double rad = mode ? 0.75 : 0.75;
			Geometry::Sphere sphere = Geometry::Sphere(0.0, rad * scale / 5.0, glass);
			if (mode)
			{
				//sphere.setCenter(Math::Vector3f(0, 0, -scale / 2));
			}
			scene.add(sphere);
		}

		if (false)
		{
			Geometry::Sphere sphere = Geometry::Sphere(0.0, 0.25 * scale / 5.0, glass);
			scene.add(sphere);
		}

		//tall block
		if(false)
		{
			Geometry::GeometryCollection* cube;
			Geometry::Material* cube_mat = spec;

			cube = new Geometry::Cube(cube_mat);
			cube->scaleX(0.3 * scale);
			cube->scaleY(0.3 * scale);
			cube->scaleZ(0.6 * scale);
			cube->rotate(Math::Quaternion<double>({ 0, 0, 1 }, 0.3));
			cube->translate({ scale * 0.2, scale * 0.15, -scale * 0.21 });
			scene.add(cube);
		}

		

		// Sets the camera
		{
			//Geometry::Camera camera = mode ? Geometry::Camera({ -scale * 2.4, 0 ,0 }, { 0, 0, 0 }, 2.f, ((double)width) / ((double)height), 1.0f) : Geometry::Camera({ -scale * 0.49, 0 ,0 }, { 0, 0, 0 }, 0.45, ((double)width) / ((double)height), 1.0f);
			Geometry::Camera camera = Geometry::Camera({ -scale * 0.49, 0 ,0 }, { 0, 0, 0 }, 0.45, ((double)width) / ((double)height), 1.0f);
			//Geometry::Camera camera({ -scale * 0.49, 0 ,0 }, { 0, 0, 0 }, 0.35, ((double)width) / ((double)height), 1.0f);
			scene.setCamera(camera);
		}
	}







	///////////////////////////////////////////////////
	//mode (type of the light:
	//	-0: a sphere just under the ceiling 
	//	-1: a cube in the center of the scene
	//	-2: a square pointing the ceiling
	///////////////////////////////////////////////////
	void initSimpleCornell(Geometry::Scene& scene, size_t width, size_t height, int mode = 3)
	{
		Geometry::Material* white = new Geometry::Lambertian<Geometry::REFLECT>(0.7);
		Geometry::Material* black = new Geometry::Lambertian<Geometry::REFLECT>(0);
		Geometry::Material* red = new Geometry::Lambertian<Geometry::REFLECT>({ 0.62, 0.061, 0.061 });
		Geometry::Material* green = new Geometry::Lambertian<Geometry::REFLECT>({ 0.122, 0.406, 0.1 });

		Geometry::Material* spec = new Geometry::Specular(0.9, 100);
		Geometry::Material* mirror = new Geometry::DeltaMirror(0.9);

		Geometry::Material* blue = new Geometry::Lambertian<Geometry::REFLECT>({ 0.1, 0.2, 0.75 });
		Geometry::Material* orange = new Geometry::Lambertian<Geometry::REFLECT>({ 0.8, 0.4, 0.1 });

		double scale = 5;
		Geometry::Cornel::init_cornell(scene, white, white, white, nullptr, red, green, scale);

		double light_size = 1;
		Geometry::Material* light = new Geometry::Material(RGBColor(16, 10, 5) / (light_size * light_size));//0.78, RGBColor(16, 10, 5)



		if (mode == 0)
		{
			Geometry::Sphere light_sphere(Math::Vector3f(0, scale * 0., scale * 0.), light_size / 2, light);
			scene.add(light_sphere);
		}
		else if (mode == 1)
		{
			Geometry::Cube* cube = new Geometry::Cube(light);
			cube->translate({ 0, 0, 0.15 * scale });
			scene.add(cube);
		}
		else if (mode == 2)//best with light tracing (or bdpt)
		{
			Geometry::Square* square = new Geometry::Square(light, { 0, 0, 0.4 * scale }, { -light_size, 0, 0 }, { 0, light_size, 0 });
			scene.add(square);
		}
		else if (mode == 3)
		{
			Geometry::Square* square = new Geometry::Square(light, { 0, 0.3 * scale, 0.4 * scale }, { light_size, 0, 0 }, { 0, light_size, 0 });
			scene.add(square);
		}
		else if (mode == 4)
		{
			Geometry::Square* square = new Geometry::Square(light, { 0, 0.0 * scale, 0.4 * scale }, { light_size, 0, 0 }, { 0, light_size, 0 });
			scene.add(square);

			Geometry::Sphere diffuse = Geometry::Sphere({ 0, 0.25 * scale, -0.3 * scale }, 0.1 * scale, white);
			scene.add(diffuse);

			Geometry::Sphere specular = Geometry::Sphere({ 0, 0.0 * scale, -0.3 * scale }, 0.1 * scale, spec);
			scene.add(specular);

			Geometry::Sphere delta = Geometry::Sphere({ 0, -0.25 * scale, -0.3 * scale }, 0.1 * scale, mirror);
			scene.add(delta);
		}


		if (mode != 4)
		{
			Geometry::Disk disk(orange, 0.0, 1, { 1, 1, 1 });
			//scene.add(disk);
		}

		if (mode != 4)
		{
			Geometry::Material* spec = new Geometry::Specular(0.8, 10);
			Geometry::Sphere ball({ 0, 0, -0.3 * scale }, light_size, spec);
			//scene.add(ball);
		}

		// Sets the camera
		{
			Geometry::Camera camera({ -scale * 2.4, 0 ,0 }, { 0, 0, 0 }, 2.f, ((double)width) / ((double)height), 1.0f);

			scene.setCamera(camera);


		}
	}










	void initCornellLamp(Geometry::Scene& scene, size_t width, size_t height)
	{
		Geometry::Material* white = new Geometry::Lambertian<Geometry::REFLECT>(0.7);
		Geometry::Material* gray = new Geometry::Lambertian<Geometry::REFLECT>(0.45);
		Geometry::Material* black = new Geometry::Lambertian<Geometry::REFLECT>(0);
		Geometry::Material* red = new Geometry::Lambertian<Geometry::REFLECT>({ 0.62, 0.061, 0.061 });
		Geometry::Material* green = new Geometry::Lambertian<Geometry::REFLECT>({ 0.122, 0.406, 0.1 });
		Geometry::Material* blue = new Geometry::Lambertian<Geometry::REFLECT>({ 0.1, 0.2, 0.75 });
		Geometry::Material* orange = new Geometry::Lambertian<Geometry::REFLECT>({ 0.8, 0.4, 0.1 });

		Geometry::Material* spec = new Geometry::Specular(0.8, 1000);
		Geometry::Material* mirror = new Geometry::DeltaMirror(0.8);

		double scale = 5;
		Geometry::Cornel::init_cornell(scene, white, white, white, nullptr, red, green, scale);

		double light_size = 0.5;
		Geometry::Material* light = new Geometry::Lambertian<Geometry::REFLECT>(0, RGBColor(16, 10, 5) / (light_size * light_size));//0.78, RGBColor(16, 10, 5)

		//add the lamps
		{
			std::vector<Math::Vector3f> cylinder_positions = { {0.3, -0.4, 0.4}, {0.3, 0.4, 0.4}, {-0.4, 0, 0.4} };
			double angle = 1;
			std::vector<Math::Quaternion<double>> rotations = { {{ 1, 0.5, 0 }, angle}, {{ -1, 0.5, 0 }, angle}, {{0, 1, 0}, -angle} };
			double radius = 0.05 * scale;

			for (int i = 0; i < 3; ++i)
			{
				Geometry::CylinderFabrice* cylinder = new Geometry::CylinderFabrice(100, radius, radius, gray, true, false);
				Math::Vector3f cylinder_pos = cylinder_positions[i] * scale;

				Math::Quaternion<double> rotation = rotations[i];
				cylinder->rotate(rotation);
				cylinder->translate(cylinder_pos);

				scene.add(cylinder);

				Math::Vector3f sphere_node = rotation.rotate({ 0, 0, 0.27 });
				Geometry::Sphere sphere((cylinder_pos + sphere_node), radius * 0.9, light);

				scene.add(sphere);
			}

		}

		{
			Geometry::Material* spec = new Geometry::Specular(0.8, 10);
			Geometry::Sphere ball({ 0, 0, -0.3 * scale }, 0.75, spec);
			scene.add(ball);
		}

		// Sets the camera
		{
			Geometry::Camera camera({ -scale * 2.4, 0 ,0 }, { 0, 0, 0 }, 2.f, ((double)width) / ((double)height), 1.0f);

			scene.setCamera(camera);
		}
	}
















	void initVeach(Geometry::Scene& scene, size_t width, size_t height, double multiplier = 25 * 10)
	{
		Geometry::Material* gray = new Geometry::Lambertian<Geometry::REFLECT>(1);

		double scale = 1.0;

		std::vector<RGBColor> sphere_colors = { {1.0, 0.1, 0.1}, {0.2, 1, 0.4}, {1, 0.15, 1}, {0.1, 0.3, 1} };
		std::vector<double> colors_multipliers = { 10, 4, 3, 2.8 };

		std::vector<double> sphere_sizes = { 0.1, 0.4, 1, 2 };

		std::vector<double> shininesses = { 1, 10, 100, 1000 };



		//add the wall and the ground
		{
			Geometry::Cornel::init_cornell(scene, nullptr, gray, gray, nullptr, nullptr, nullptr, scale * 50);
		}


		//add the shperes
		{
			Math::Vector3f d = { 0, -scale * 7, 0 };
			Math::Vector3f base = { scale * 18, 0, -5 * scale };
			base -= d * 1.5;
			for (int i = 0; i < sphere_sizes.size(); ++i)
			{
				Geometry::Material* mat = new Geometry::Material(sphere_colors[i] * colors_multipliers[i] * 5 * scale / (sphere_sizes[i]));
				if (true)//spheres
				{
					Geometry::Sphere sphere = Geometry::Sphere(base + d * i, sphere_sizes[i], mat);
					scene.add(sphere);
				}
				else
				{
					Geometry::Square* square = new Geometry::Square(mat, base + d * i, Math::Vector3f(1, 0, -0.5) * sphere_sizes[i], Math::Vector3f(0, 1, 0) * sphere_sizes[i]);
					scene.add(square);
				}
			}
		}

		std::vector<double> angles = { 0, 20, 35, 50 };
		std::vector<double> corrections = { 30, 0, -20, 0 };

		for (int i = 0; i < angles.size(); ++i)
		{
			Math::Vector3f base = { -1 * scale, 0, 5 * scale };
			Math::Vector3f T = { 0, 0, -26 * scale };
			Geometry::Square* cube = new Geometry::Square(new Geometry::Specular(1, shininesses[i] * multiplier));
			//Geometry::Square* cube = new Geometry::Square(new Geometry::Lambertian<Geometry::REFLECT>(1));
			cube->rotate(Math::Quaternion<double>({ 0, -1, 0 }, rad(corrections[i])));
			cube->translate(T);

			cube->scaleY(25 * scale);
			cube->scaleX(5 * scale);
			cube->rotate(Math::Quaternion<double>({ 0, -1, 0 }, rad(angles[i])));
			cube->translate(base);


			scene.add(cube);
		}


		//add a light
		{
			Geometry::Material* light = new Geometry::Material(10 * scale);
			Geometry::Square* square = new Geometry::Square(light);

			square->rotate(Math::Quaternion<double>({ 0, -1, 0 }, rad(45)));
			square->scale(scale * 10 * 2);
			square->translate({ -100 * scale, 0, 15 * scale });

			//scene.add(square);
		}


		// Sets the camera
		{
			Geometry::Camera camera({ -scale * 25, 0 ,-10 * scale }, { 0, 0, -15 * scale }, 1.f, ((double)width) / ((double)height), 1.0f);
			scene.setCamera(camera);
		}
	}


}