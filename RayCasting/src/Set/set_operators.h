#ifndef _SET_OPERATORS_H
#define _SET_OPERATORS_H

#include <vector>
#include <list>
#include <set>
#include <iostream>

////////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn	template <class T> std::vector<T> & operator+=(std::vector<T> & v, T const & s)
///
/// \brief	Adds an element at the end of a given ::std::vector<T>
///
/// \author	F. Lamarche, Université de Rennes 1
/// \date	03/12/2013
///
/// \tparam	T	Generic type parameter.
/// \param [in,out]	v	The vector in which an element is added.
/// \param	s		 	The added element.
///
/// \return	A reference to the modified vector a.k.a v.
////////////////////////////////////////////////////////////////////////////////////////////////////
template <class T>
std::vector<T> & operator+=(std::vector<T> & v, T const & s)
{
	v.push_back(s) ;
	return v ;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn	template <class T, class Cmp> inline std::vector<T> & operator+=(std::vector<T> & v,
/// 	std::set<T, Cmp> const & s)
///
/// \brief	Adds the elements contained in set s in vector v.
///
/// \author	F. Lamarche, Université de Rennes 1
/// \date	03/12/2013
///
/// \tparam	T  	Generic type parameter.
/// \tparam	Cmp	Type of the compare.
/// \param [in,out]	v	The vector.
/// \param	s		 	The set.
///
/// \return	A reference to v.
////////////////////////////////////////////////////////////////////////////////////////////////////
template <class T, class Cmp>
inline std::vector<T> & operator+=(std::vector<T> & v, std::set<T, Cmp> const & s)
{
	int vSize = v.size() ;
	v.resize(s.size()+vSize) ;
	int cpt=vSize ; 
	for(typename std::set<T, Cmp>::iterator it=s.begin() ; it!=s.end() ; it++, cpt++)
	{
		v[cpt]=*it ;
	}
	return v ;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn	template <class T> inline std::vector<T> & operator+=(std::vector<T> & v,
/// 	std::vector<T> const & s)
///
/// \brief	Contatenates elements in s at the end of v.
///
/// \author	F. Lamarche, Université de Rennes 1
/// \date	03/12/2013
///
/// \tparam	T	Generic type parameter.
/// \param [in,out]	v	The modified vector.
/// \param	s		 	The elements to add.
///
/// \return	The result of the operation.
////////////////////////////////////////////////////////////////////////////////////////////////////
template <class T>
inline std::vector<T> & operator+=(std::vector<T> & v, std::vector<T> const & s)
{
	int vSize = v.size() ;
	v.resize(s.size()+vSize) ;
	int cpt=vSize ; 
	for(typename std::vector<T>::const_iterator it=s.begin() ; it!=s.end() ; it++, cpt++)
	{
		v[cpt]=*it ;
	}
	return v ;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn	template <class T, class Cmp> inline bool in(T const & e, std::set<T, Cmp> const & s)
///
/// \brief	Tests if an element is contained in a given set.
///
/// \author	F. Lamarche, Université de Rennes 1
/// \date	03/12/2013
///
/// \tparam	T  	Generic type parameter.
/// \tparam	Cmp	Type of the compare.
/// \param	e	The set.
/// \param	s	The tested element.
///
/// \return	true if it succeeds, false if it fails.
////////////////////////////////////////////////////////////////////////////////////////////////////
template <class T, class Cmp>
inline bool in(T const & e, std::set<T, Cmp> const & s)
{
	return (s.find(e)!=s.end()) ;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn	template <class T, class Cmp> inline std::set<T, Cmp> & operator+=(std::set<T, Cmp> & s,
/// 	T const & e)
///
/// \brief	Adds an element in a set.
///
/// \author	F. Lamarche, Université de Rennes 1
/// \date	03/12/2013
///
/// \tparam	T  	Generic type parameter.
/// \tparam	Cmp	Type of the compare.
/// \param [in,out]	s	The set.
/// \param	e		 	The element to add.
///
/// \return	The result of the operation.
////////////////////////////////////////////////////////////////////////////////////////////////////
template <class T, class Cmp>
inline std::set<T, Cmp> & operator+=(std::set<T, Cmp> & s, T const & e)
{
	typename std::set<T, Cmp>::iterator pos=s.find(e) ;
	if(pos==s.end())
	{ s.insert(e) ; }
	return s ;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn	template <class T, class Cmp> inline std::set<T, Cmp> & operator-=(std::set<T, Cmp> & s,
/// 	T const & e)
///
/// \brief	Removes an element from a set.
///
/// \author	F. Lamarche, Université de Rennes 1
/// \date	03/12/2013
///
/// \tparam	T  	Generic type parameter.
/// \tparam	Cmp	Type of the compare.
/// \param [in,out]	s	The set.
/// \param	e		 	The element that should be removed.
///
/// \return	The result of the operation.
////////////////////////////////////////////////////////////////////////////////////////////////////
template <class T, class Cmp>
inline std::set<T, Cmp> & operator-=(std::set<T, Cmp> & s, T const & e)
{
	typename std::set<T, Cmp>::iterator pos=s.find(e) ;
	if(pos!=s.end())
	{ s.erase(pos) ; }
	return s ;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn	template <class T, class Cmp> std::set<T, Cmp> & operator+=(std::set<T, Cmp> & s,
/// 	std::list<T> const & l)
///
/// \brief	Adds a collection of elements into a set.
/// 		
/// \author	F. Lamarche, Université de Rennes 1
/// \date	03/12/2013
///
/// \tparam	T  	Generic type parameter.
/// \tparam	Cmp	Type of the compare.
/// \param [in,out]	s	The set.
/// \param	l		 	The elements to add.
///
/// \return	A reference to s.
////////////////////////////////////////////////////////////////////////////////////////////////////
template <class T, class Cmp>
std::set<T, Cmp> & operator+=(std::set<T, Cmp> & s, std::list<T> const & l)
{
	for(typename std::list<T>::const_iterator it=l.begin() ; it!=l.end() ; it++)
	{
		s+=(*it) ;
	}
	return s ;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn	template <class T1, class T2, class Cmp> std::set<T1, Cmp> & operator-=(std::set<T1, Cmp> & s,
/// 	std::list<T2> const & l)
///
/// \brief	Removes a collection of elements from a set.
///
/// \author	F. Lamarche, Université de Rennes 1
/// \date	03/12/2013
///
/// \tparam	T1 	Generic type parameter.
/// \tparam	T2 	Generic type parameter.
/// \tparam	Cmp	Type of the comparison operator.
/// \param [in,out]	s	The set.
/// \param	l		 	The collection of removed elements.
///
/// \return	The result of the operation.
////////////////////////////////////////////////////////////////////////////////////////////////////
template <class T1, class T2, class Cmp>
std::set<T1, Cmp> & operator-=(std::set<T1, Cmp> & s, std::list<T2> const & l)
{
	for(typename std::list<T2>::const_iterator it=l.begin() ; it!=l.end() ; it++)
	{
		s-=(*it) ;
	}
	return s ;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn	template <class T, class Cmp> std::set<T, Cmp> & operator+=(std::set<T, Cmp> & s,
/// 	std::vector<T> const & l)
///
/// \brief	Adds a collection of elements in a set.
///
/// \author	F. Lamarche, Université de Rennes 1
/// \date	03/12/2013
///
/// \tparam	T  	Generic type parameter.
/// \tparam	Cmp	Type of the compare.
/// \param [in,out]	s	The set.
/// \param	l		 	The elements that should be added.
///
/// \return A refernce to s.
////////////////////////////////////////////////////////////////////////////////////////////////////
template <class T, class Cmp>
std::set<T, Cmp> & operator+=(std::set<T, Cmp> & s, std::vector<T> const & l)
{
	for(typename std::vector<T>::const_iterator it=l.begin() ; it!=l.end() ; it++)
	{
		s+=(*it) ;
	}
	return s ;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn	template <class T1, class T2, class Cmp> std::set<T1, Cmp> & operator-=(std::set<T1, Cmp> & s,
/// 	std::vector<T2> const & l)
///
/// \brief	Removes a collection of elements from a given set.
///
/// \author	F. Lamarche, Université de Rennes 1
/// \date	03/12/2013
///
/// \tparam	T1 	Generic type parameter.
/// \tparam	T2 	Generic type parameter.
/// \tparam	Cmp	Type of the compare.
/// \param [in,out]	s	The set.
/// \param	l		 	The elements to remove.
///
/// \return	A reference to s.
////////////////////////////////////////////////////////////////////////////////////////////////////
template <class T1, class T2, class Cmp>
std::set<T1, Cmp> & operator-=(std::set<T1, Cmp> & s, std::vector<T2> const & l)
{
	for(typename std::vector<T2>::const_iterator it=l.begin() ; it!=l.end() ; it++)
	{
		s-=(*it) ;
	}
	return s ;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn	template <class T, class Cmp> std::set<T, Cmp> & operator+=(std::set<T, Cmp> & s1,
/// 	std::set<T, Cmp> const & s2)
///
/// \brief	Adds a collection of elements in a given set.
///
/// \author	F. Lamarche, Université de Rennes 1
/// \date	03/12/2013
///
/// \tparam	T  	Generic type parameter.
/// \tparam	Cmp	Type of the compare.
/// \param [in,out]	s1	The set
/// \param	s2		  	The elements that should be added
///
/// \return	The result of the operation.
////////////////////////////////////////////////////////////////////////////////////////////////////
template <class T, class Cmp>
std::set<T, Cmp> & operator+=(std::set<T, Cmp> & s1, std::set<T, Cmp> const & s2)
{
	for(typename std::set<T, Cmp>::const_iterator it=s2.begin() ; it!=s2.end() ; it++)
	{
		s1+=*it ;
	}
	return s1 ;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn	template <class T, class Cmp> std::set<T, Cmp> & operator-=(std::set<T, Cmp> & s1,
/// 	std::set<T, Cmp> const & s2)
///
/// \brief	Removes a collection of elements from a given set.
///
/// \author	F. Lamarche, Université de Rennes 1
/// \date	03/12/2013
///
/// \tparam	T  	Generic type parameter.
/// \tparam	Cmp	Type of the compare.
/// \param [in,out]	s1	The set
/// \param	s2		  	The elements to remove
///
/// \return	The result of the operation.
////////////////////////////////////////////////////////////////////////////////////////////////////
template <class T, class Cmp>
std::set<T, Cmp> & operator-=(std::set<T, Cmp> & s1, std::set<T, Cmp> const & s2)
{
	for(typename std::set<T, Cmp>::const_iterator it=s2.begin() ; it!=s2.end() ; it++)
	{
		s1-=*it ;
	}
	return s1 ;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn	template <class T, class Cmp> std::set<T, Cmp> & operator*=(std::set<T, Cmp> & s1,
/// 	std::set<T, Cmp> const & s2)
///
/// \brief	Computes the intersection between s1 and s2 and stores the result in s1.
///
/// \author	F. Lamarche, Université de Rennes 1
/// \date	03/12/2013
///
/// \tparam	T  	Generic type parameter.
/// \tparam	Cmp	Type of the compare.
/// \param [in,out]	s1	The first set
/// \param	s2		  	The second set
///
/// \return	A reference to s1
////////////////////////////////////////////////////////////////////////////////////////////////////
template <class T, class Cmp>
std::set<T, Cmp> & operator*=(std::set<T, Cmp> & s1, std::set<T, Cmp> const & s2)
{
	for(typename std::set<T, Cmp>::iterator it=s1.begin() ; it!=s1.end() ; it++)
	{
		if(!in(*it, s2))
		{
			s1.erase(it) ;
		}
	}
	return s1 ; 
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn	template <class T, class Cmp> inline std::set<T, Cmp> operator+(std::set<T, Cmp> const & s1,
/// 	std::set<T, Cmp> const & s2)
///
/// \brief	Computes the union of to sets.
///
/// \author	F. Lamarche, Université de Rennes 1
/// \date	03/12/2013
///
/// \tparam	T  	Generic type parameter.
/// \tparam	Cmp	Type of the compare.
/// \param	s1	The first set
/// \param	s2	The second set
///
/// \return	The union of s1 and s2
////////////////////////////////////////////////////////////////////////////////////////////////////
template <class T, class Cmp>
inline std::set<T, Cmp> operator+(std::set<T, Cmp> const & s1, std::set<T, Cmp> const & s2)
{
	if(s1.size()<s2.size())
		return s2+s1 ;

	std::set<T, Cmp> result(s1) ;
	for(typename std::set<T, Cmp>::const_iterator it=s2.begin() ; it!=s2.end() ; it++)
	{
		if(s1.find(*it)==s1.end())
			result.insert(*it) ;
	}
	return result ;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn	template <class T, class Cmp> inline std::set<T, Cmp> operator-(std::set<T, Cmp> const & s1,
/// 	std::set<T, Cmp> const & s2)
///
/// \brief	Computes the difference between two sets.
///
/// \author	F. Lamarche, Université de Rennes 1
/// \date	03/12/2013
///
/// \tparam	T  	Generic type parameter.
/// \tparam	Cmp	Type of the compare.
/// \param	s1	The first set
/// \param	s2	The second set
///
/// \return	s1-s2.
////////////////////////////////////////////////////////////////////////////////////////////////////
template <class T, class Cmp>
inline std::set<T, Cmp> operator-(std::set<T, Cmp> const & s1, std::set<T, Cmp> const & s2)
{
	std::set<T, Cmp> result(s1) ;

	for(typename std::set<T, Cmp>::const_iterator it=s2.begin() ; it!=s2.end() ; it++)
	{
		if(s1.find(*it)!=s1.end())
			result.erase(*it) ;
	}
	return result ;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn	template <class T, class Cmp> inline std::set<T, Cmp> operator*(std::set<T, Cmp> const & s1,
/// 	std::set<T, Cmp> const & s2)
///
/// \brief	Intersection between two sets.
///
/// \author	F. Lamarche, Université de Rennes 1
/// \date	03/12/2013
///
/// \tparam	T  	Generic type parameter.
/// \tparam	Cmp	Type of the compare.
/// \param	s1	The first set
/// \param	s2	The second set
///
/// \return	The intersection of s1 and s2
////////////////////////////////////////////////////////////////////////////////////////////////////
template <class T, class Cmp>
inline std::set<T, Cmp> operator*(std::set<T, Cmp> const & s1, std::set<T, Cmp> const & s2)
{
	std::set<T, Cmp> result ;

	if(s2.size()<s1.size())
		return s2*s1 ;

	for(typename std::set<T, Cmp>::const_iterator it=s1.begin() ; it!=s1.end() ; it++)
	{
		if(s2.find(*it)!=s2.end())
			result.insert(*it) ;
	}
	return result ;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn	template <class T, class Cmp> inline std::set<T, Cmp> operator^(std::set<T, Cmp> const & s1,
/// 	std::set<T, Cmp> const & s2)
///
/// \brief	Symetric difference between two sets
///
/// \author	F. Lamarche, Université de Rennes 1
/// \date	03/12/2013
///
/// \tparam	T  	Generic type parameter.
/// \tparam	Cmp	Type of the compare.
/// \param	s1	The first set
/// \param	s2	The second set
///
/// \return	The symmetric difference between s1 and s2
////////////////////////////////////////////////////////////////////////////////////////////////////
template <class T, class Cmp>
inline std::set<T, Cmp> operator^(std::set<T, Cmp> const & s1, std::set<T, Cmp> const & s2)
{
	if(s2.size()<s1.size())
		return s2^s1 ;

	std::set<T, Cmp> result(s2) ;
	for(typename std::set<T, Cmp>::const_iterator it=s1.begin() ; it!=s1.end() ; it++)
	{
		if(s2.find(*it)==s2.end())
		{ result.insert(*it) ; }
		else
		{ result.erase(*it) ; }
	}

	return result ;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn	template <class Type> inline std::ostream & operator<<(std::ostream & out, const std::set<Type> & s)
///
/// \brief	Outputs a set in an ostream
///
/// \author	F. Lamarche, Université de Rennes 1
/// \date	03/12/2013
///
/// \param [in,out]	out	The out.
/// \param	s		   	The set
///
/// \return	A reference to out
////////////////////////////////////////////////////////////////////////////////////////////////////
template <class Type>
inline std::ostream & operator<<(std::ostream & out, const std::set<Type> & s)
{
	out<<"{ " ;
	for(std::set<int>::const_iterator it=s.begin() ; it!=s.end() ; it++)
	{
		out<<(*it)<<" " ;
	}
	out<<"}" ;
}


#endif
