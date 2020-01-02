#pragma once

#include <string>
#include <vector>
#include <Image/Image.h>
#include <iostream>
#include <map>
#include <Auto/TestScenes.h>
#include <Auto/RenderResult.h>
#include <omp.h>
#include <filesystem>
#include <Image/ImWrite.h>

#include <Integrators/DirectIntegrator.h>
#include <Integrators/RegularIntegrators.h>
#include <Integrators/RayTracingIntegrator.h>
#include <Integrators/Integrator.h>
#include <Integrators/PathTracingIntegrator.h>
#include <Integrators/ZIntegrator.h>
#include <Integrators/LightIntegrator.h>
#include <Integrators/MISPathTracingIntegrator.h>
#include <Integrators/BidirectionalIntegrator.h>
#include <ctime>

namespace Auto
{
	class Auto
	{
	protected:


		class TestIntegrator
		{
		public:
			Integrator::Integrator* integrator;
			bool active;
			std::string name;
			RenderResult res;

			TestIntegrator() :
				integrator(nullptr),
				active(false),
				name(),
				res()
			{}

			TestIntegrator(Integrator::Integrator * ptr, bool b, std::string const& namae):
				integrator(ptr),
				active(b),
				name(namae),
				res()
			{}



			bool compute(Geometry::Scene const& scene, size_t w, size_t h)
			{
				if (active)
				{
					std::cout << "Computing " << name << std::endl;
					integrator->render(scene, w, h, res);
					return true;
				}
				return false;
			}
		};


		class TestResult
		{
		public:
			std::string name;
			Image::Image<Geometry::RGBColor> result;
			Image::Image<Geometry::RGBColor> dif;
			double time;
		};

		

		struct Option
		{

			int nthreads = 8;
			int width = 500, height = 500;

			std::string main_path = RESULT_FOLDER;

			std::string scene_name;

			
			unsigned int sample_per_pixel=1;
			unsigned int max_bounce=10;
			unsigned int light_division=1;
			double alpha = 0.9;

			bool use_npt = false;
			unsigned int naive_times = 1;

			bool use_pt = false;
			bool use_ptMIS = false;
			bool use_lt = false;
			bool use_bdpt = false;
			bool use_nbdpt = false;


			Option(std::vector<std::string> const& args)
			{
				for (int i = 0; i < args.size(); ++i)
				{
					std::string arg = args[i];
					if (arg == "-resolution")
					{
						std::string next = args[i + 1];
						width = std::stoi(next);
						++i;

						next = args[i + 1];
						height = std::stoi(next);
						++i;
					}
					else if (arg == "-t")
					{
						std::string const& next = args[i + 1];
						nthreads = std::stoi(next);
						++i;
					}
					else if (arg == "-path")
					{
						std::string const& next = args[i + 1];
						main_path = next;
						++i;
					}
					else if (arg == "-scene")
					{
						std::string const& next = args[i + 1];
						scene_name = next;
						++i;
					}
					else if (arg == "-spp")
					{
						std::string const& next = args[i + 1];
						sample_per_pixel = std::stoi(next);
						++i;
					}
					else if (arg == "-mb")
					{
						std::string const& next = args[i + 1];
						max_bounce = std::stoi(next);
						++i;
					}
					else if (arg == "-ld")
					{
						std::string const& next = args[i + 1];
						light_division = std::stoi(next);
						++i;
					}
					else if (arg == "-a")
					{
						std::string const& next = args[i + 1];
						alpha = std::stod(next);
						++i;
					}
					else if (arg == "-pt")
					{
						use_pt = true;
					}
					else if (arg == "-npt")
					{
						use_npt = true;
						if (i != args.size() - 1)
						{
							std::string const& next = args[i + 1];
							//extra option for npt
							if (next[0] != '-')
							{
								naive_times = std::stoi(next);
								++i;
							}
							
						}
						
					}
					else if (arg == "-ptMIS")
					{
						use_ptMIS = true;
					}
					else if (arg == "-lt")
					{
						use_lt = true;
					}
					else if (arg == "-bdpt")
					{
						use_bdpt = true;
					}
					else if (arg == "-nbdpt")
					{
						use_nbdpt = true;
					}
					else
					{
						std::cout << "Unknown option: " << arg << std::endl;
					}

				}
				
				if (main_path.back() == '\\')
				{
					main_path.pop_back();
				}
			}



			template <class out_t>
			void print(out_t& out)const
			{
				out << "General settings:" << std::endl;
				std::cout << "Rendering in " << width << "x" << height << " with " << nthreads << " threads" << std::endl;
				std::cout << "Scene: " << scene_name << std::endl;
				std::cout << "Results in " << main_path << std::endl;
				std::cout << "spp " << sample_per_pixel << std::endl;
				std::cout << "mb " << max_bounce << std::endl;
				std::cout << "ld " << light_division << std::endl;
				std::cout << "alpha " << alpha << std::endl;
			}
		};


		static void writeResult(std::vector<TestResult> const& res, RenderResult const& ref, Option const& op)
		{
			std::string str = strResults(res, ref, op);

			if (!std::filesystem::is_directory(op.main_path) || !std::filesystem::exists(op.main_path))
			{
				std::filesystem::create_directory(op.main_path);
			}

			time_t now = time(NULL);
			tm ltm = *localtime(&now);
			std::string nows;
			nows += '_';
			nows += ltm.tm_year;
			nows += '_';
			nows += ltm.tm_mon+1;
			nows += '_';
			nows += ltm.tm_mday;
			nows += '_';
			nows += ltm.tm_hour;
			nows += '_';
			nows += ltm.tm_min;
			nows += '_';
			nows += ltm.tm_sec;

			std::filesystem::path path = op.main_path + "\\" + op.scene_name + nows;

			
			
			if (!std::filesystem::is_directory(path) || !std::filesystem::exists(path))
			{
				std::filesystem::create_directory(path);
			}
			
			
			std::ofstream file(path.string() + "\\" + "Report.txt");

			file << str;

			file.close();

			std::cout << str << std::endl;


			Image::ImWrite::write(ref.image, path.string() + "\\" + "ref.exr");
			
			for (TestResult const& tr : res)
			{
				Image::ImWrite::write(tr.result, path.string() + "\\" + tr.name + ".exr");
			}


		}


		static std::string strResults(std::vector<TestResult> const& res, RenderResult const& ref, Option const& op)
		{
			std::stringstream str;
			str << "Computing with " << op.sample_per_pixel << " samples" << std::endl;
			std::string line = "----------------------------------------------------";
			for (TestResult const& tr : res)
			{
				str << line << std::endl;
				str << tr.name << std::endl;
				str << "mean: "<<tr.result.mean()<<", ref mean: "<<ref.image.mean() << ", dif mean: "<< tr.dif.mean() << std::endl;
				str << "eqm: " << Image::Image<Geometry::RGBColor>::eqm(tr.result, ref.image) << std::endl;
				str << "etr: " << Image::Image<Geometry::RGBColor>::etm(tr.result, ref.image) << std::endl;
				str << "time: " << tr.time << ", ref: " << ref.time << ", relative: " << tr.time / ref.time << std::endl;
			}

			return str.str();
		}




	public:
	

		static int __main__(int argc, char** argv)
		{
			std::vector<std::string> args;
			args.reserve(argc);
			for (int i = 1; i < argc; ++i)
			{
				args.push_back(argv[i]);
			}

			Option options(args);

			options.print(std::cout);


			omp_set_num_threads(options.nthreads);

			using SceneInitializer = std::function<void (Geometry::Scene &, size_t, size_t)>;

			std::map<std::string, SceneInitializer> scenes;
			{
				scenes["cornell"] = [](Geometry::Scene & scene, size_t w, size_t h)
				{
					initRealCornell(scene, w, h, 0, 0, 0);
				};
				scenes["cornellColors"] = [](Geometry::Scene& scene, size_t w, size_t h)
				{
					initRealCornell(scene, w, h, 0, 1, 0);
				};
				scenes["cornellColorsSpec"] = [](Geometry::Scene& scene, size_t w, size_t h)
				{
					initRealCornell(scene, w, h, 1, 1, 0);
				};
				scenes["cornellColorsMirror"] = [](Geometry::Scene& scene, size_t w, size_t h)
				{
					initRealCornell(scene, w, h, 2, 1, 0);
				};
				scenes["cornellLamp"] = [](Geometry::Scene & scene, size_t w, size_t h)
				{
					initCornellLamp(scene, w, h);
				};
				scenes["simpleCornell"] = [](Geometry::Scene & scene, size_t w, size_t h)
				{
					initSimpleCornell(scene, w, h);
				};
				scenes["veach"] = [](Geometry::Scene & scene, size_t w, size_t h)
				{
					initVeach(scene, w, h);
				};
				scenes["test"] = [](Geometry::Scene& scene, size_t w, size_t h)
				{
					initTest(scene, w, h);
				};
			}

			if (scenes.find(options.scene_name) == scenes.end())
			{
				std::cerr << "Could not open scene: " << options.scene_name << std::endl;
				return 1;
			}

			Geometry::Scene scene;

			scenes[options.scene_name](scene, options.width, options.height);

			
			std::cout << "Building the acceleration structure" << std::endl;
			scene.preCompute(8, 1, 1);
			std::cout << "Done!" << std::endl;
			scene.printStats();

			scene.compute_light_samplers(options.light_division);

			scene.m_camera.resolution = options.width * options.height;
			
			Integrator::NaivePathTracingIntegrator npti(options.sample_per_pixel * options.naive_times, options.width, options.height);
			Integrator::IterativePathTracingIntegrator ipti(options.sample_per_pixel, options.width, options.height);
			Integrator::MISPathTracingIntegrator mispti(options.sample_per_pixel, options.width, options.height);
			Integrator::LightIntegrator lti(options.sample_per_pixel, options.width, options.height);
			Integrator::BidirectionalIntegrator bdpti(options.sample_per_pixel, options.width, options.height);
			//Integrator::NaiveBidirectionalIntegrator nbdpti(options.sample_per_pixel, options.width, options.height);

			std::vector<TestIntegrator> integrators = { 
				{&npti, options.use_npt, "Naive Path Tracing"},
				{&ipti, options.use_pt, "Path Tracing"},
				{&mispti, options.use_ptMIS, "MIS Path Tracing"},
				{&lti, options.use_lt, "Light Tracing"},
				{&bdpti, options.use_bdpt, "Bidirectional Path Tracing"},
				//{&nbdpti, options.use_nbdpt, "Naive Bidirectional Path Tracing"},
			};

			for (TestIntegrator& inte : integrators)
			{
				inte.integrator->setLen(options.max_bounce);
				inte.integrator->m_alpha = options.alpha;
			}

			

			///////////////////////////////////////////
			// Which integrator result should be considered as the reference? 
			// For now, use the most reliable result, which is, in order
			// npt
			// ipt
			///////////////////////////////////////////

			if (!(options.use_npt || options.use_pt))
			{
				std::cout << "No reference result can be rendered!" << std::endl;
				return 1;
			}
			int ref_index = options.use_npt ? 0 : 1;
			RenderResult& ref = integrators[ref_index].res;

			//compute the iterators results
			for (TestIntegrator& test : integrators)
			{
				test.compute(scene, options.width, options.height);
			}

			std::vector<TestResult> results;
			results.reserve(integrators.size());
			//compare the results to the reference
			int index = 0;
			for (TestIntegrator& test : integrators)
			{
				if (index != ref_index && test.active)
				{
					TestResult res;
					res.result = test.res.image;
					res.name = test.name;
					res.dif = ref.image - test.res.image;
					res.time = test.res.time;
					results.push_back(res);
				}
				++index;
			}


			writeResult(results, ref, options);
			

			return 0;
		}
	};
}