#pragma once

#include <Visualizer/Visualizer.h>
#include <Geometry/BoundedStack.h>
#include <Math/Vector.h>
#include <Math/Vectorf.h>
#include <Math\Constant.h>
#include <functional>

namespace Geometry
{
	class BRDFViewer
	{
	protected:
		using uint = unsigned int;
		uint m_width, m_height;

		uint m_samples;

		

	public:


		BRDFViewer(uint width=1000, uint height=1000, uint samples=1000) :
			m_width(width),
			m_height(height),
			m_samples(samples)
		{

		}


		


		void run()
		{
			Visualizer::Visualizer visu(m_width, m_height);


			
		}

	};
}