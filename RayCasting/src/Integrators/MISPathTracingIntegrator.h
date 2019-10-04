#pragma once

#include <Integrators/RayTracingBaseIntegrator.h>
#include <cassert>


namespace Integrator
{
	class MISPathTracingIntegrator: public RayTracingBaseIntegrator
	{
	protected:

	public:

		MISPathTracingIntegrator(unsigned int sample_per_pixel, unsigned int width, unsigned int height) :
			RayTracingBaseIntegrator(sample_per_pixel, width, height)
		{}


		static __forceinline void balance(__out double& w0, __out double& w1, const double p00, const double p01, const double p10, const double p11)
		{
			w0 = p00 / (p00 + p10);
			w1 = p11 / (p01 + p11);			
		}

		static __forceinline void power(__out double& w0, __out double& w1, const double p00, const double p01, const double p10, const double p11, double beta = 2)
		{
			w0 = pow(p00, beta) / (pow(p00, beta)+ pow(p10, beta));
			w1 = pow(p11, beta) / (pow(p01, beta) + pow(p11, beta));
		}


		RGBColor MISAddDirectIllumination(Scene const& scene, Hit const& hit, Math::Sampler& sampler)const
		{
			return 0;
			// TODO
		}

		

		RGBColor sendRay(Scene const& scene, Ray const& pray, Math::Sampler& sampler)const final override
		{
			return 0;
			// TODO
		}


		
	};
}