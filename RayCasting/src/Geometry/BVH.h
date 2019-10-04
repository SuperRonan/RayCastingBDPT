#pragma once

#include <struct/static_dual_tree.h>
#include <functional>
#include <vector>
#include <limits>
#include <string>
#include <algorithm>
#include <Geometry/BoundingBox.h>
#include <Geometry/CastedRay.h>
#include <cassert>
#include <Geometry/Primitive.h>
#include <tbb/parallel_sort.h>
#include <deque>
#include <execution>
#include <Geometry/Intersection.h>
#include <tbb/parallel_invoke.h>

namespace Geometry
{
	

	template <class Primitive>
	class BVH
	{
	protected:
		
		using PrimitiveCollection = std::vector<const Primitive*>;

		using TreeType = StaticDualTree<2, BoundingBox, PrimitiveCollection>;

		TreeType * m_tree = nullptr;

		struct Node
		{
			PrimitiveCollection primitives;
			BoundingBox box;
			Node* left = nullptr, *right = nullptr;

			Node(TreeType const& other) :
				primitives(other.getLeafValue()),
				box(other.getNodeValue())
			{}

			Node()
			{}

			bool isLeaf()const
			{
				return !primitives.empty();
			}
		};

		std::vector<Node> m_tree_vector;

	public:

		BVH()
		{}

		bool empty()const
		{
			return m_tree == nullptr;
		}

		BoundingBox boundingBox()const
		{
			return m_tree->getNodeValue();
		}


		bool compute(unsigned int n_max, double ct, double ci, PrimitiveCollection const& primitives)
		{
			if (primitives.empty())
			{
				return false;
			}

			PrimitiveCollection _primitives = primitives;

			computeBody(n_max, ct, ci, _primitives, boundingBox(_primitives), m_tree, 0);

			//fill the vector 
			m_tree_vector.resize(m_tree->size());

			m_tree_vector[0] = *m_tree;
			unsigned int res = buildVector(&m_tree_vector[0], 1, m_tree);
			assert(res == m_tree_vector.size());
			return true;
		}

		bool intersection(Ray const& ray, Intersection<Primitive>& res)const
		{
			if (empty())
			{
				return false;
			}
			CastedRay<Primitive> cray(ray);
			//iterativeIntersection(cray);
			iterativeIntersectionVector(cray);
			res = cray.intersectionFound();
			return cray.validIntersectionFound();
			
			double t_entry, t_exit;
			if (boundingBox().intersect(cray, 0, std::numeric_limits<double>::max(), t_entry, t_exit))
			{
				intersection(cray, m_tree, t_entry, t_exit);
				if (cray.validIntersectionFound())
				{
					res = cray.intersectionFound();
				}
			}
			return cray.validIntersectionFound();
		}
		
		unsigned int countBB(Ray const& ray)const
		{
			if (empty())
			{
				return 0;
			}
			return countBBIterativeVector(ray);
			return countBBIterative(ray);
			return countBB(ray, m_tree);
		}

	protected:

		//index reprensents the first free index in the vector
		unsigned int buildVector(Node * parent, unsigned int index, const TreeType * tree_node)
		{
			unsigned int res = index;
			if (!parent->isLeaf())
			{
				m_tree_vector[index] = *tree_node->operator[](0);
				m_tree_vector[index+1] = *tree_node->operator[](1);

				parent->left = &m_tree_vector[index];
				parent->right = &m_tree_vector[index + 1];

				res += 2;

				res = buildVector(&m_tree_vector[index], res, tree_node->operator[](0));
				res = buildVector(&m_tree_vector[index+1], res, tree_node->operator[](1));
			}
			return res;
		}

		unsigned int countBB(Ray const& ray, const TreeType const* const current_node)const
		{
			double _, __;
			if (current_node->getNodeValue().intersect(ray, 0, std::numeric_limits<double>::max(), _, __))
			{
				unsigned int res = 1 + (current_node->isLeaf() ? 0 : countBB(ray, (*current_node)[0]) + countBB(ray, (*current_node)[1]));
				return res;
			}
			return 0;
		}

		unsigned int countBBIterative(Ray const& ray)const
		{
			double _, __;
			BoundedStack<30, const TreeType*> stack;
			const TreeType* node = m_tree;
			unsigned int res = 0;
			while (true)
			{
				if (node->getNodeValue().intersect(ray, 0, std::numeric_limits<double>::max(), _, __))
				{
					++res;
					if (node->isLeaf())
					{
						if (stack.empty())
						{
							break;
						}
						node = stack.top();
						stack.pop();
						continue;
					}
					stack.push((*node)[1]);
					node = (*node)[0];
				}
				else
				{
					if (stack.empty())
					{
						break;
					}
					node = stack.top();
					stack.pop();
				}
			}
			return res;
		}

		unsigned int countBBIterativeVector(Ray const& ray)const
		{
			double _, __;
			BoundedStack<30, const Node*> stack;
			const Node* node = &m_tree_vector[0];
			unsigned int res = 0;
			while (true)
			{
				if (node->box.intersect(ray, 0, std::numeric_limits<double>::max(), _, __))
				{
					++res;
					if (node->isLeaf())
					{
						if (stack.empty())
						{
							break;
						}
						node = stack.top();
						stack.pop();
						continue;
					}
					stack.push(node->right);
					node = node->left;
				}
				else
				{
					if (stack.empty())
					{
						break;
					}
					node = stack.top();
					stack.pop();
				}
			}
			return res;
		}


		
		void iterativeIntersection(CastedRay<Primitive>& cray)const
		{
			struct Param
			{
				const TreeType* node;
				double t_entry;
				double t_exit;
			};
			Param current = { m_tree, 0, 0 };
			if (current.node->getNodeValue().intersect(cray, 0, std::numeric_limits<double>::max(), current.t_entry, current.t_exit))
			{
				BoundedStack<30, Param> params;
				while (true)
				{
					if (current.node->isLeaf())
					{
						for (const Primitive* pri : current.node->getLeafValue())
						{
							cray.intersect(*pri);
						}
						//stop we can go up
						while (!params.empty())
						{
							if (params.top().t_entry > cray.intersectionFound().t())//discard boxes farther than the current result
							{
								params.pop();
							}
							else
							{
								break;
							}
						}
						if (params.empty())
						{
							break;
						}
						current = params.top();
						current.t_exit = std::min(current.t_exit, cray.intersectionFound().t());
						params.pop();
						continue;
					}
					//else
					//go down the tree
					double t_entry_l, t_exit_l, t_entry_r, t_exit_r;

					bool intersect_bb_l = (*current.node)[0]->getNodeValue().intersect(cray, current.t_entry, current.t_exit, t_entry_l, t_exit_l),
						intersect_bb_r = (*current.node)[1]->getNodeValue().intersect(cray, current.t_entry, current.t_exit, t_entry_r, t_exit_r);

					if (intersect_bb_l)
					{
						if (intersect_bb_r)
						{
							//go down to the closest box, and keep memory of the other one in the stack
							if (t_entry_l < t_entry_r)
							{
								//start by the left one
								{
									//fill the stack
									params.end()->node = (*current.node)[1];
									params.end()->t_entry = t_entry_r;
									params.end()->t_exit = t_exit_r;
									params.grow();
								}
								current.t_entry = t_entry_l;
								current.t_exit = t_exit_l;
								current.node = (*current.node)[0];
							}
							else
							{
								//start by the right one
								{
									//fill the stack
									params.end()->node = (*current.node)[0];
									params.end()->t_entry = t_entry_l;
									params.end()->t_exit = t_exit_l;
									params.grow();
								}
								current.t_entry = t_entry_r;
								current.t_exit = t_exit_r;
								current.node = (*current.node)[1];
							}
						}
						else
						{
							//go down to the left only
							current.node = (*current.node)[0];
							current.t_entry = t_entry_l;
							current.t_exit = t_exit_l;
						}
					}
					else
					{
						if (intersect_bb_r)
						{
							//go down to the right only
							current.node = (*current.node)[1];
							current.t_entry = t_entry_r;
							current.t_exit = t_exit_r;
						}
						else
						{
							//stop we can go up
							while (!params.empty())
							{
								if (params.top().t_entry > cray.intersectionFound().t())
								{
									params.pop();
								}
								else
								{
									break;
								}
							}
							if (params.empty())
							{
								break;
							}
							current = params.top();
							current.t_exit = std::min(current.t_exit, cray.intersectionFound().t());
							params.pop();
						}
					}//intersect bb left
				}//while
			}//if intersection with the surrounding box
		}

		void iterativeIntersectionVector(CastedRay<Primitive>& cray)const
		{
			struct Param
			{
				const Node* node;
				double t_entry;
				double t_exit;
			};
			Param current = { &m_tree_vector[0], 0, 0 };
			if (current.node->box.intersect(cray, 0, std::numeric_limits<double>::max(), current.t_entry, current.t_exit))
			{
				BoundedStack<30, Param> params;
				while (true)
				{
					if (current.node->isLeaf())
					{
						for (const Primitive* pri : current.node->primitives)
						{
							cray.intersect(*pri);
						}
						//stop we can go up
						while (!params.empty())
						{
							if (params.top().t_entry > cray.intersectionFound().t())//discard boxes farther than the current result
							{
								params.pop();
							}
							else
							{
								break;
							}
						}
						if (params.empty())
						{
							break;
						}
						current = params.top();
						current.t_exit = std::min(current.t_exit, cray.intersectionFound().t());
						params.pop();
						continue;
					}
					//else
					//go down the tree
					double t_entry_l, t_exit_l, t_entry_r, t_exit_r;

					bool intersect_bb_l = current.node->left->box.intersect(cray, current.t_entry, current.t_exit, t_entry_l, t_exit_l),
						intersect_bb_r = current.node->right->box.intersect(cray, current.t_entry, current.t_exit, t_entry_r, t_exit_r);

					if (intersect_bb_l)
					{
						if (intersect_bb_r)
						{
							//go down to the closest box, and keep memory of the other one in the stack
							if (t_entry_l < t_entry_r)
							{
								//start by the left one
								{
									//fill the stack
									params.end()->node = current.node->right;
									params.end()->t_entry = t_entry_r;
									params.end()->t_exit = t_exit_r;
									params.grow();
								}
								current.t_entry = t_entry_l;
								current.t_exit = t_exit_l;
								current.node = current.node->left;
							}
							else
							{
								//start by the right one
								{
									//fill the stack
									params.end()->node = current.node->left;
									params.end()->t_entry = t_entry_l;
									params.end()->t_exit = t_exit_l;
									params.grow();
								}
								current.t_entry = t_entry_r;
								current.t_exit = t_exit_r;
								current.node = current.node->right;
							}
						}
						else
						{
							//go down to the left only
							current.node = current.node->left;
							current.t_entry = t_entry_l;
							current.t_exit = t_exit_l;
						}
					}
					else
					{
						if (intersect_bb_r)
						{
							//go down to the right only
							current.node = current.node->right;
							current.t_entry = t_entry_r;
							current.t_exit = t_exit_r;
						}
						else
						{
							//stop we can go up
							while (!params.empty())
							{
								if (params.top().t_entry > cray.intersectionFound().t())
								{
									params.pop();
								}
								else
								{
									break;
								}
							}
							if (params.empty())
							{
								break;
							}
							current = params.top();
							current.t_exit = std::min(current.t_exit, cray.intersectionFound().t());
							params.pop();
						}
					}//intersect bb left
				}//while
			}//if intersection with the surrounding box
		}


		void intersection(
			CastedRay<Primitive>& cray,
			const TreeType const* const current_node,
			const double t_entry,
			const double t_exit
			)const
		{
			if (current_node->isLeaf())
			{
				PrimitiveCollection const& primitives = current_node->getLeafValue();
				//TODO avx sse
				for (const Primitive * pri : primitives)
				{
					cray.intersect(*pri);
				}
			}
			else
			{
				double t_entry_l, t_exit_l, t_entry_r, t_exit_r;

				//TODO remove t_entry and t_exit
				bool intersect_bb_l = (*current_node)[0]->getNodeValue().intersect(cray, t_entry, t_exit, t_entry_l, t_exit_l);
				bool intersect_bb_r = (*current_node)[1]->getNodeValue().intersect(cray, t_entry, t_exit, t_entry_r, t_exit_r);

				if (intersect_bb_l)
				{
					if (intersect_bb_r)
					{
						//begin by the "closest" box
						if (t_entry_l < t_entry_r) // begin by the left box
						{
							intersection(cray, (*current_node)[0], t_entry_l, t_exit_l);
							if (cray.validIntersectionFound())
							{
								if (cray.intersectionFound().t() > t_entry_r)
								{
									//TODO t_exit_r = cray.intersection.t
									intersection(cray, (*current_node)[1], t_entry_r, cray.intersectionFound().t());
								}
							}
							else //no left box intersection
							{
								intersection(cray, (*current_node)[1], t_entry_r, t_exit_r);
							}
						}
						else //begin by the right box
						{
							intersection(cray, (*current_node)[1], t_entry_r, t_exit_r);
							if (cray.validIntersectionFound())
							{
								if (cray.intersectionFound().t() > t_entry_l)
								{
									//TODO t_exit_l = cray.intersection.t
									intersection(cray, (*current_node)[0], t_entry_l, cray.intersectionFound().t());
								}
							}
							else //no right intersection
							{
								intersection(cray, (*current_node)[0], t_entry_l, t_exit_l);
							}
						}
					}
					else
					{
						intersection(cray, (*current_node)[0], t_entry_l, t_exit_l);
					}
				}
				else
				{
					if (intersect_bb_r)
					{
						intersection(cray, (*current_node)[1], t_entry_r, t_exit_r);
					}
					else // none
					{

					}
				}
			}
		}




		BoundingBox boundingBox(PrimitiveCollection const& primitives)const
		{
			BoundingBox res;
			for (const Primitive* primitive : primitives)
			{
				const auto pbox = primitive->box();
				res.update(pbox);
			}
			return res;
		}

		double SAH(double ct, double ci, unsigned int nl, unsigned int nr, double sb, double sbl, double sbr)const
		{
			return ct + ci * (sbl / sb * nl + sbr / sb * nr);
		}


		double SAHPartition(
			double ci, double ct,
			PrimitiveCollection const& sorted_primitives,
			double top_box_surface,
			PrimitiveCollection & p_left,
			PrimitiveCollection & p_right,
			BoundingBox & bb_left,
			BoundingBox & bb_right,
			bool parallel=false
		)const
		{
			assert(sorted_primitives.size() > 1);
			std::deque<const Primitive*> tmp_left, tmp_right;

			double best = std::numeric_limits<double>::max();

			BoundingBox* right_bb = new BoundingBox[sorted_primitives.size()], left_bb, tmp;
			
			//fill the right bb
			for (int i = sorted_primitives.size() - 1; i >= 0; --i)
			{
				tmp_right.push_front(sorted_primitives[i]);
				tmp.update(sorted_primitives[i]->box());
				right_bb[i] = tmp;
			}

			BoundingBox best_left, best_right;

			unsigned int besk = 0;
			for (unsigned int k = 0; k < sorted_primitives.size(); ++k)
			{
				tmp_left.push_back(tmp_right.front());
				left_bb.update(tmp_right.front()->box());
				tmp_right.pop_front();

				double left_bb_surface = left_bb.surface();
				double right_bb_surface = right_bb[k].surface();

				double sah = SAH(ct, ci, tmp_left.size(), tmp_right.size(), top_box_surface, left_bb_surface, right_bb_surface);

				if (sah < best)
				{
					best = sah;
					best_left = left_bb;
					besk = k;
				}
			}
			tmp_left.clear();
			tmp_right.clear();

			best_right = right_bb[besk];
			p_left.resize(besk+1);
			p_right.resize(sorted_primitives.size() - besk-1);
			auto split = sorted_primitives.cbegin() + besk+1;

			if (parallel)
			{
				std::copy(std::execution::par, sorted_primitives.cbegin(), split, p_left.begin());
				std::copy(std::execution::par, split, sorted_primitives.cend(), p_right.begin());
			}
			else
			{
				std::copy(sorted_primitives.cbegin(), split, p_left.begin());
				std::copy(split, sorted_primitives.cend(), p_right.begin());
			}

			bb_left = best_left;
			bb_right = right_bb[besk];

			delete[] right_bb;

			assert(!p_left.empty());
			assert(!p_right.empty());
			assert(p_left.size() + p_right.size() == sorted_primitives.size());

			return best;
		}


		void bestSAHPartition(
			double ci, 
			double ct,
			PrimitiveCollection const& primitives,
			BoundingBox const& bb, 
			PrimitiveCollection & p_left,
			PrimitiveCollection & p_right,
			BoundingBox & bb_left,
			BoundingBox & bb_right,
			bool parallel=false
		)const
		{
			assert(primitives.size() > 1);
			assert(bb.valid());
			double best_sah = std::numeric_limits<double>::max();

			double bb_surface = bb.surface();

			for (unsigned int axis = 0; axis < 3; ++axis)
			{
				PrimitiveCollection axis_sorted = primitives;
				
				auto comparator = [axis](const Primitive* p, const Primitive* q)
				{
					return p->center()[axis] < q->center()[axis];
				};

				if (parallel)
				{
					tbb::parallel_sort(axis_sorted.begin(), axis_sorted.end(), comparator);
				}
				else
				{
					std::sort(axis_sorted.begin(), axis_sorted.end(), comparator);
				}

				PrimitiveCollection axis_p_left, axis_p_right;
				BoundingBox axis_bb_left, axis_bb_right;
				
				double axis_sah = SAHPartition(ci, ct, axis_sorted, bb_surface, axis_p_left, axis_p_right, axis_bb_left, axis_bb_right, parallel);

				if (axis_sah < best_sah)
				{
					best_sah = axis_sah;
					p_left = axis_p_left;
					p_right = axis_p_right;
					bb_left = axis_bb_left;
					bb_right = axis_bb_right;
				}
			}

			assert(!p_left.empty());
			assert(!p_right.empty());
		}


		void computeBody(unsigned int n_max, double ct, double ci, PrimitiveCollection& primitives, 
			BoundingBox const& current_bb, TreeType *& current_tree_pos, unsigned int depth)
		{
			if (primitives.size() <= n_max)
			{
				//make the bounding box just a little bit bigger
				current_tree_pos = TreeType::makeNewLeaf(current_bb.larger(), primitives);
			}
			else
			{
				current_tree_pos = TreeType::makeNewNode(current_bb);

				PrimitiveCollection p_left, p_right;
				BoundingBox bb_left, bb_right;
				bool parallel = depth == 0;
				bestSAHPartition(ci, ct, primitives, current_bb, p_left, p_right, bb_left, bb_right, parallel);

				primitives = PrimitiveCollection();
				if (depth < 4)
				{
					auto left = [&]() {computeBody(n_max, ct, ci, p_left, bb_left, (*current_tree_pos)[0], depth + 1); };
					auto right = [&]() {computeBody(n_max, ct, ci, p_right, bb_right, (*current_tree_pos)[1], depth + 1); };
					tbb::parallel_invoke(left, right);
				}
				else
				{
					
					computeBody(n_max, ct, ci, p_left, bb_left, (*current_tree_pos)[0], depth + 1);
					computeBody(n_max, ct, ci, p_right, bb_right, (*current_tree_pos)[1], depth + 1);
				}
			}
		}

	};
}