#pragma once

#include <iostream>
#include <chrono>
#include <stack>

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


//roughly counts the number of rays / sec and prints it in the terminal
//#define COUNT_RAYS

//use direct sampling of light sources in path tracing, to be removed with the implementation of naive path tracer
#define SAMPLE_DIRECT

//use a little trick to get less noise with direct lighting, but fully biased
//#define TRICK_DIRECT

//don't use the russian roulette before a given depth (usualy ~2-3)
#define LATE_RUSSIAN

//when computing the lighting, add the contribution of all the lights, if undef: pick one light at random
#define SAMPLE_ALL_LIGHTS

//multiply the final result is the Visualizer x10 for the tone mapping
#define MULT_10

//when loading materials with the 3ds loader, try to convert to Lambert or specular instead of always using Phong
#define CONVERT_MATERIALS

//when converting the materials, for them to be lambertian
#define FORCE_LAMBERT

#include <string>

//shortcut path to save the results
static std::string RESULT_FOLDER = "../../results/";

#define OMP_DYNAMIC_FOR __pragma(omp parallel for schedule(dynamic))

#define OMP_STATIC_FOR __pragma(omp parallel for schedule(static))

#define OMP_PARALLEL_FOR OMP_DYNAMIC_FOR

//in real time mode, use the time as seed, so the noise will change at each frame 
//#define TIME_SEED

//Use the same seed for all the pixels at each frame, which reduces significantly the noise (to 0 naiiiiiiiiiiisu)
//#define SAMPLER_BIAS

//for MIS, show the weight of MIS
//#define MIS_SHOW_WEIGHTS

//for debug of bdpt, shows the result of only one technique, see the bdpt integrator::computeSample for more settings
//#define BDPT_SINGLE_TECHNIQUE

//deactivate MIS for bdpt <=> weight = 1
//#define BDPT_NO_WEIGHT