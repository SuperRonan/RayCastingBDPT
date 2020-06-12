#pragma once

#include <exception>
#include <string>
#include <cassert>

///////////////////////////
//An exception class dscribing for not implemented functionality
/////////////////////////////
class not_implemented_exception : public std::exception
{
public:

	not_implemented_exception(std::string const& msg) :
		std::exception(), m_msg(msg)
	{}

	virtual ~not_implemented_exception()
	{}

	not_implemented_exception & operator=(std::exception const& e)
	{
		m_msg = e.what();
	}

	virtual const char * what()const noexcept
	{
		return m_msg.c_str();
	}


protected:
	std::string m_msg;
	
};





/////////////////////////////////////
// Static: always the same number of sons for each node
// Dual: Carries two types, one for the nodes and one for the leaf
/////////////////////////////////////
template<unsigned int N, class node_type, class leaf_type>
class StaticDualTree
{
public:

	StaticDualTree(bool leaf, node_type const& n_v, leaf_type const& l_v) :
		m_is_leaf(leaf), node_value(n_v), leaf_value(l_v)
	{
		if (m_is_leaf)
		{

		}
		else
		{
			//sons = new StaticDualTree*[size];
			for (unsigned int i = 0; i < N; ++i)
			{
				sons[i] = nullptr;
			}
		}
	}

	~StaticDualTree()
	{
		if (!m_is_leaf)
		{
			for (unsigned int i = 0; i < N; ++i)
			{
				delete sons[i];
			}
			//delete[] sons;
		}
	}

	static StaticDualTree* makeNewLeaf(node_type const& n_v, leaf_type const& l_v)
	{
		return new StaticDualTree(true, n_v, l_v);
	}

	static StaticDualTree* makeNewNode(node_type const& n_v)
	{
		return new StaticDualTree(false, n_v, leaf_type());
	}

	StaticDualTree *& operator[](unsigned int i)
	{
		assert(i >= 0);
		assert(i < N);
		assert(!m_is_leaf);
		return sons[i];
	}
	StaticDualTree * const& operator[](unsigned int i)const
	{
		assert(i >= 0);
		assert(i < N);
		assert(!m_is_leaf);
		return sons[i];
	}

	bool isLeaf()const
	{
		return m_is_leaf;
	}

	bool isNode()const
	{
		return !m_is_leaf;
	}

	leaf_type & getLeafValue()
	{
		assert(m_is_leaf);
		return leaf_value;
	}
	leaf_type const& getLeafValue()const
	{
		//assert(m_is_leaf);
		return leaf_value;
	}

	node_type & getNodeValue()
	{
		return node_value;
	}
	node_type const& getNodeValue()const
	{
		return node_value;
	}

	size_t size()const
	{
		size_t res = 1;
		if (isNode())
		{
			for (unsigned int i=0; i<N; ++i)
			{
				res += sons[i]->size();
			}
		}
		return res;
	}
	
protected:

	bool m_is_leaf;

	StaticDualTree * sons[N];

	node_type node_value;

	leaf_type leaf_value;

};





//////////////////////////////////
// Deprecated representations of a tree
/////////////////////////////////


/*
template <unsigned int size, class node_type, class leaf_type>
class static_dual_tree
{
public:

	static_dual_tree(node_type const& v=node_type()):
		node_value(v)
	{}

	virtual ~static_dual_tree()
	{}

	virtual bool isLeaf()const = 0;
	virtual bool isNode()const = 0;


	//virtual const static_dual_tree *& operator[](unsigned int index)const = 0;
	virtual static_dual_tree *& operator[](unsigned int index) = 0;

	
	virtual node_type const& getNodeValue()const
	{
		return node_value;
	}
	
	virtual node_type & getNodeValue()
	{
		return node_value;
	}

	//virtual leaf_type const& getLeafValue()const = 0;
	virtual leaf_type & getLeafValue() = 0;


protected:

	node_type node_value;

private:

};

template <unsigned int size, class node_type, class leaf_type>
class static_dual_leaf : public static_dual_tree<size, node_type, leaf_type>
{
public:

	static_dual_leaf(node_type const& n_v=node_type(), leaf_type const& v=leaf_type()):
		static_dual_tree(n_v), leaf_value(v)
	{}

	virtual ~static_dual_leaf()
	{
		
	}

	virtual bool isLeaf()const
	{
		return true;
	}
	virtual bool isNode()const
	{
		return false;
	}

	
	
	virtual static_dual_tree *& operator[](unsigned int index)
	{
		throw not_implemented_exception("Cannot get the son of a leaf!");
	}

	virtual leaf_type & getLeafValue()
	{
		return leaf_value;
	}
protected:

	leaf_type leaf_value;
};

template <unsigned int size, class node_type, class leaf_type>
class static_dual_node : public static_dual_tree<size, node_type, leaf_type>
{
public:

	static_dual_node(node_type const& v) :
		static_dual_tree(v), sons(new static_dual_tree<size, node_type, leaf_type>*[size])
	{
		
		for (unsigned int i = 0; i < size; ++i)
		{
			sons[i] = nullptr;
		}
	}

	

	virtual ~static_dual_node()
	{
		for (unsigned int i = 0; i < size; ++i)
		{
			if (sons[i] != nullptr)
			{
				delete sons[i];
			}
		}
		delete[] sons;
	}

	virtual bool isLeaf()const
	{
		return false;
	}
	virtual bool isNode()const
	{
		return true;
	}

	

	
	
	virtual static_dual_tree *& operator[](unsigned int index)
	{
		assert(index >= 0);
		assert(index < size);

		return sons[index];
	}


	virtual leaf_type & getLeafValue()
	{
		throw not_implemented_exception("Cannot get the leaf value of a node!");
	}


protected:

	static_dual_tree<size, node_type, leaf_type> ** sons;
};
*/