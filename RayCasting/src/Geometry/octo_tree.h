#pragma once

namespace Geometry
{
	class octo_tree
	{
	public:

		virtual ~octo_tree() = 0;

	protected:

	};


	class octo_node: public octo_tree
	{
	protected:

	};

	class octo_leaf : public octo_tree
	{

	};
}