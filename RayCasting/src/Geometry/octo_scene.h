#pragma once

#include "Scene.h"

namespace Geometry
{
	class octo_scene: public Scene
	{
	protected:

	public:

		octo_scene(Visualizer::Visualizer * visu)
			:Scene(visu)
		{}


		void pre_compute(unsigned int max_depth = 4, unsigned int depth = 0)
		{
			if (depth >= max_depth)
			{
				return;
			}
			Math::Vector3f middle = (m_sceneBoundingBox.max + m_sceneBoundingBox) / 2;
			
		}
	};
}