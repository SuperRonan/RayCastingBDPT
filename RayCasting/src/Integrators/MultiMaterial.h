#pragma once

#include <Geometry/Material.h>

namespace Geometry
{
	template <class ... Mat>
	class MultMaterial: public Mat...
	{
		
	};
}