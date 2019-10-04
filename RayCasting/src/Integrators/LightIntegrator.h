#pragma once

#include <Integrators/LightIntegratorBase.h>
#include <cassert>

namespace Integrator
{
	class LightIntegrator final : public LightIntegratorBase
	{
	protected:

	public:

		LightIntegrator(unsigned int sample_per_pixel, unsigned int width, unsigned int height) :
			LightIntegratorBase(sample_per_pixel, width, height)
		{

		}

		void traceLight(Scene const& scene, LightVertexStack& lvs, Math::Sampler& sampler)const final override
		{
			// TODO
			
		}

	};
}