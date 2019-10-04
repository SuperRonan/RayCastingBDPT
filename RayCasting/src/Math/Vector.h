#ifndef _Rennes1_Math_Vector_H
#define _Rennes1_Math_Vector_H

#include <iostream>
#include <cmath>
#include <algorithm>
#include <cassert>

namespace Math
{

	////////////////////////////////////////////////////////////////////////////////////////////////////
	/// \class	Vector
	///
	/// \brief	Class representing a N-dimensional vector. 
	///
	/// \param Float The scalar type.
	/// \param dimensions The dimension of the vector
	/// 
	/// \author	Fabrice Lamarche, University Of Rennes 1
	/// \date	23/01/2010
	////////////////////////////////////////////////////////////////////////////////////////////////////
	template <class Float, int dimensions>
	class Vector
	{
	protected:
		//! The scalars. 
		Float m_vector[dimensions] ;

	public:
		//! Iterator on vector coordinates. 
		typedef Float * iterator ;
		//! Const iterator on vector coordinates. 
		typedef const Float * const_iterator ;

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	Vector(Float * floatArray)
		///
		/// \brief	Constructor from an array of scalars. Scalars are copied into the vector.
		///
		/// \author	Fabrice Lamarche, University Of Rennes 1
		/// \date	23/01/2010
		///
		/// \param [in]	floatArray Array of scalars for initialization. 
		////////////////////////////////////////////////////////////////////////////////////////////////////
		Vector(const Float * floatArray)
		{
			for(int cpt=0 ; cpt<dimensions ; ++cpt)
			{
				m_vector[cpt] = floatArray[cpt] ;
			}
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	Vector(Float const & value)
		///
		/// \brief	Constructor that initializes each coordinate with the given value. 
		///
		/// \author	Fabrice Lamarche, University Of Rennes 1
		/// \date	23/01/2010
		///
		/// \param value	the value of each coordinate. 
		////////////////////////////////////////////////////////////////////////////////////////////////////
		Vector(Float const & value = Float())
		{
			for(int cpt=0 ; cpt<dimensions ; ++cpt)
			{
				m_vector[cpt] = value ;
			}
		}

		/*
		template <class ... T>
		Vector(const T & ... args)
		{
			static_assert(sizeof...(args) <= dimensions);
			for (unsigned int i = 0; i < sizeof...(args); ++i)
			{
				m_vector[i] = args[i];
			}
		}
		*/

		template <class T, class Q>
		Vector(T t, Q q)
		{
			m_vector[0] = t;
			m_vector[1] = q;
		}


		template <class T, class Q, class R>
		Vector(T t, Q q, R r)
		{
			m_vector[0] = t;
			m_vector[1] = q;
			m_vector[2] = r;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	Vector(Vector<Float, dimensions> const & v)
		///
		/// \brief	Copy constructor. 
		///
		/// \author	Fabrice Lamarche, University Of Rennes 1
		/// \date	23/01/2010
		///
		/// \param [in,out]	v the v. 
		////////////////////////////////////////////////////////////////////////////////////////////////////
		Vector(Vector<Float, dimensions> const & v)
		{
			for(int cpt=0 ; cpt<dimensions ; cpt++)
			{ m_vector[cpt] = v.m_vector[cpt] ; }
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	Vector(Vector<Float2, dimensions> const & v)
		///
		/// \brief	Generic copy constructor. 
		///
		/// \author	Fabrice Lamarche, University Of Rennes 1
		/// \date	23/01/2010
		///
		/// \param [in,out]	v	the v. 
		////////////////////////////////////////////////////////////////////////////////////////////////////
		template <class Float2>
		Vector(Vector<Float2, dimensions> const & v)
		{
			for(int cpt=0 ; cpt<dimensions ; cpt++)
			{ m_vector[cpt] = Float(v[cpt]) ; }
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	Float * Vector::getBuffer()
		///
		/// \brief	Gets the coordinates as an array of Float.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	28/11/2015
		///
		/// \return	null if it fails, else the buffer.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		Float * getBuffer() 
		{
			return m_vector;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	const Float * Vector::getBuffer() const
		///
		/// \brief	Gets the coordinates as an array of Float.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	28/11/2015
		///
		/// \return	null if it fails, else the buffer.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		const Float * getBuffer() const
		{
			return m_vector ;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	Float norm() const
		///
		/// \brief	Return the norm of the vector. 
		///
		/// \author	Fabrice Lamarche, University Of Rennes 1
		/// \date	23/01/2010
		///
		/// \return	the norm of the vector.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		Float norm() const
		{
			return (Float)sqrt(this->norm2()) ;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	Float norm2() const
		///
		/// \brief	Returns the squared norm of the vector. 
		///
		/// \author	Fabrice Lamarche, University Of Rennes 1
		/// \date	23/01/2010
		///
		/// \return	the squared norm of the vector. 
		////////////////////////////////////////////////////////////////////////////////////////////////////
		Float norm2() const
		{
			return (*this)*(*this) ;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	const Float & operator[] (int index) const
		///
		/// \brief	 Access to the value of each coordinate. 
		///
		/// \author	Fabrice Lamarche, University Of Rennes 1
		/// \date	23/01/2010
		///
		/// \param	index	Zero-based index of the coordinate. 
		///
		/// \return	the value of the 'index' nth coordinate. 
		////////////////////////////////////////////////////////////////////////////////////////////////////
		const Float & operator[] (int index) const
		{ return m_vector[index] ; }


		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	Float & operator[] (int index)
		///
		/// \brief	 Access to the value of each coordinate. 
		///
		/// \author	Fabrice Lamarche, University Of Rennes 1
		/// \date	23/01/2010
		///
		/// \param	index	Zero-based index of the coordinate. 
		///
		/// \return	the value of the 'index' nth coordinate. 
		////////////////////////////////////////////////////////////////////////////////////////////////////
		Float & operator[] (int index) 
		{ return m_vector[index] ; }

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	Vector<Float, dimensions+1> push(Float const & coord) const
		///
		/// \brief	Creates a new vector of dimension 'dimensions'+1 by pushing the coordinate 'coord' on back of the vector. 
		///
		/// \author	Fabrice Lamarche, University Of Rennes 1
		/// \date	25/01/2010
		///
		/// \param coord	the coordinate to push. 
		///
		/// \return	the newly created vector. 
		////////////////////////////////////////////////////////////////////////////////////////////////////
		Vector<Float, dimensions+1> push(Float const & coord) const 
		{
			Math::Vector<Float, dimensions+1> result ;
			for(int cpt=0 ; cpt<dimensions ; ++cpt)
			{
				result[cpt] = m_vector[cpt] ; 
			}
			result[dimensions] = coord ;
			return result ;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	Vector<Float, dimensions> operator+ (Vector<Float, dimensions> const & v) const
		///
		/// \brief	Addition operator. 
		///
		/// \author	Fabrice Lamarche, University Of Rennes 1
		/// \date	23/01/2010
		///
		/// \param [in,out]	v	the v. 
		///
		/// \return	The result of the operation. 
		////////////////////////////////////////////////////////////////////////////////////////////////////
		Vector<Float, dimensions> operator+ (Vector<Float, dimensions> const & v) const
		{
			Vector<Float, dimensions> result ;
			for(int cpt=0 ; cpt<dimensions ; cpt++)
			{ result[cpt] = m_vector[cpt]+v.m_vector[cpt] ;	}
			return result ;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	Vector<Float, dimensions> operator- (Vector<Float, dimensions> const & v) const
		///
		/// \brief	Substraction operator. 
		///
		/// \author	Fabrice Lamarche, University Of Rennes 1
		/// \date	23/01/2010
		///
		/// \param [in,out]	v	the v. 
		///
		/// \return	The result of the operation. 
		////////////////////////////////////////////////////////////////////////////////////////////////////
		Vector<Float, dimensions> operator- (Vector<Float, dimensions> const & v) const
		{
			Vector<Float,dimensions> result ;
			for(int cpt=0 ; cpt<dimensions ; cpt++)
			{ result[cpt] = m_vector[cpt]-v.m_vector[cpt] ;	}
			return result ;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	Vector<Float, dimensions> operator- () const
		///
		/// \brief	Negation operator. 
		///
		/// \author	Fabrice Lamarche, University Of Rennes 1
		/// \date	23/01/2010
		///
		/// \return	The result of the operation. 
		////////////////////////////////////////////////////////////////////////////////////////////////////
		Vector<Float, dimensions> operator- () const
		{
			Vector<Float, dimensions> result ;
			for(int cpt=0 ; cpt<dimensions ; cpt++)
			{ result[cpt] = -m_vector[cpt] ;	}
			return result ;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	Float operator * (Vector<Float, dimensions> const & v) const
		///
		/// \brief	Scalar product. 
		///
		/// \author	Fabrice Lamarche, University Of Rennes 1
		/// \date	23/01/2010
		///
		/// \param [in,out]	v	the v. 
		///
		/// \return	The result of the operation. 
		////////////////////////////////////////////////////////////////////////////////////////////////////
		Float operator * (Vector<Float, dimensions> const & v) const
		{
			Float result ;
			result = m_vector[0]*v.m_vector[0] ;
			for(int cpt=1 ; cpt<dimensions ; cpt++)
			{ result = result + m_vector[cpt]*v.m_vector[cpt] ; }
			return result ;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	Vector<Float, dimensions> operator* (Float const & v) const
		///
		/// \brief	Muliplication operator. 
		///
		/// \author	Fabrice Lamarche, University Of Rennes 1
		/// \date	23/01/2010
		///
		/// \param [in,out]	v	the v. 
		///
		/// \return	The result of the operation. 
		////////////////////////////////////////////////////////////////////////////////////////////////////
		Vector<Float, dimensions> operator* (Float const & v) const
		{
			Vector<Float, dimensions> result ;
			for(int cpt=0 ; cpt<dimensions ; cpt++)
			{ result[cpt] = m_vector[cpt]*v ;	}
			return result ;			
		}


		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	Vector<Float, dimensions> operator/ (Float const & v) const
		///
		/// \brief	Division operator. 
		///
		/// \author	Fabrice Lamarche, University Of Rennes 1
		/// \date	23/01/2010
		///
		/// \param [in,out]	v	the v. 
		///
		/// \return	The result of the operation. 
		////////////////////////////////////////////////////////////////////////////////////////////////////
		Vector<Float, dimensions> operator/ (Float const & v) const
		{
			Vector<Float, dimensions> result ;
			for(int cpt=0 ; cpt<dimensions ; cpt++)
			{ result[cpt] = m_vector[cpt]/v ;	}
			return result ;			
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	Vector<Float, dimensions> & operator+= (Vector<Float, dimensions> const & v)
		///
		/// \brief	Assignment by addition operator. 
		///
		/// \author	Fabrice Lamarche, University Of Rennes 1
		/// \date	23/01/2010
		///
		/// \param [in,out]	v	the v. 
		///
		/// \return	A shallow copy of this object. 
		////////////////////////////////////////////////////////////////////////////////////////////////////
		Vector<Float, dimensions> & operator+= (Vector<Float, dimensions> const & v)
		{
			for(int cpt=0 ; cpt<dimensions ; ++cpt)
			{
				m_vector[cpt] += v.m_vector[cpt] ;
			}
			return (*this) ;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	Vector<Float, dimensions> & operator-= (Vector<Float, dimensions> const & v)
		///
		/// \brief	Assignment by subtraction operator. 
		///
		/// \author	Fabrice Lamarche, University Of Rennes 1
		/// \date	23/01/2010
		///
		/// \param [in,out]	v	the v. 
		///
		/// \return	A shallow copy of this object. 
		////////////////////////////////////////////////////////////////////////////////////////////////////
		Vector<Float, dimensions> & operator-= (Vector<Float, dimensions> const & v)
		{
			for(int cpt=0 ; cpt<dimensions ; ++cpt)
			{
				m_vector[cpt] -= v.m_vector[cpt] ;
			}
			return (*this) ;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	Vector<Float, dimensions> & operator*= (Float const & factor)
		///
		/// \brief	Assignment by muliplication operator. 
		///
		/// \author	Fabrice Lamarche, University Of Rennes 1
		/// \date	23/01/2010
		///
		/// \param [in,out]	factor	the factor. 
		///
		/// \return	A shallow copy of this object. 
		////////////////////////////////////////////////////////////////////////////////////////////////////
		Vector<Float, dimensions> & operator*= (Float const & factor)
		{
			for(int cpt=0 ; cpt<dimensions ; ++cpt)
			{
				m_vector[cpt] *= factor ;
			}
			return (*this) ;				
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	Vector<Float, dimensions> & operator/= (Float const & factor)
		///
		/// \brief	Assignment by division operator. 
		///
		/// \author	Fabrice Lamarche, University Of Rennes 1
		/// \date	23/01/2010
		///
		/// \param [in,out]	factor	the factor. 
		///
		/// \return	A shallow copy of this object. 
		////////////////////////////////////////////////////////////////////////////////////////////////////
		Vector<Float, dimensions> & operator/= (Float const & factor)
		{
			for(int cpt=0 ; cpt<dimensions ; ++cpt)
			{
				m_vector[cpt] /= factor ;
			}
			return (*this) ;				
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	Vector<Float, dimensions> & operator= (Vector<Float, dimensions> const & v)
		///
		/// \brief	Copy operator. 
		///
		/// \author	Fabrice Lamarche, University Of Rennes 1
		/// \date	23/01/2010
		///
		/// \param [in,out]	v	the v. 
		///
		/// \return	A shallow copy of this object. 
		////////////////////////////////////////////////////////////////////////////////////////////////////
		Vector<Float, dimensions> & operator= (Vector<Float, dimensions> const & v) 
		{
			for(int cpt=0 ; cpt<dimensions ; cpt++)
			{
				m_vector[cpt]=v.m_vector[cpt] ;
			}
			return (*this) ;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	Vector<Float, dimensions> & operator= (Vector<Float2, dimensions> const & v)
		///
		/// \brief	Generic copy operator. 
		///
		/// \author	Fabrice Lamarche, University Of Rennes 1
		/// \date	23/01/2010
		///
		/// \param [in,out]	v	the v. 
		///
		/// \return	A shallow copy of this object. 
		////////////////////////////////////////////////////////////////////////////////////////////////////
		template <class Float2>
		Vector<Float, dimensions> & operator= (Vector<Float2, dimensions> const & v) 
		{
			for(int cpt=0 ; cpt<dimensions ; cpt++)
			{
				m_vector[cpt]=(Float)v[cpt] ;
			}
			return (*this) ;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	Vector<Float, dimensions> & operator= (Float const & s)
		///
		/// \brief	Copy operator. 
		///
		/// \author	Fabrice Lamarche, University Of Rennes 1
		/// \date	23/01/2010
		///
		/// \param [in,out]	s	the. 
		///
		/// \return	A shallow copy of this object. 
		////////////////////////////////////////////////////////////////////////////////////////////////////
		Vector<Float, dimensions> & operator= (Float const & s)
		{
			for(int cpt=0 ; cpt<dimensions ; cpt++)
			{ m_vector[cpt] = s ; }
			return *this ;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	Vector<Float, dimensions> inv() const
		///
		/// \brief	Pseudo inverse of the vector. 
		///
		/// \author	Fabrice Lamarche, University Of Rennes 1
		/// \date	23/01/2010
		///
		/// \return	. 
		////////////////////////////////////////////////////////////////////////////////////////////////////
		Vector<Float, dimensions> inv() const
		{ return (*this)/norm2() ; }

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \property	bool operator<(Vector<Float, dimensions> const & v) const
		///
		/// \brief	 Comparison operator (lexicographical order). 
		///
		/// \return	The result of the comparison. 
		////////////////////////////////////////////////////////////////////////////////////////////////////
		bool operator<(Vector<Float, dimensions> const & v) const
		{
			for(int cpt=0 ; cpt<dimensions ; cpt++)
			{
				if(m_vector[cpt]==v[cpt]) continue ;
				return m_vector[cpt]<v[cpt] ;
			}
			return false ;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	bool operator==(Vector<Float, dimensions> const & v) const
		///
		/// \brief	Equality operator. 
		///
		/// \author	Fabrice Lamarche, University Of Rennes 1
		/// \date	23/01/2010
		///
		/// \param [in,out]	v	the v. 
		///
		/// \return	true if the parameters are considered equivalent. 
		////////////////////////////////////////////////////////////////////////////////////////////////////
		bool operator==(Vector<Float, dimensions> const & v) const
		{
			for(int cpt=0 ; cpt<dimensions ; cpt++)
			{
				if(m_vector[cpt]!=v[cpt]) return false ;
			}
			return true ;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	bool operator!=(Vector<Float, dimensions> const & v) const
		///
		/// \brief	Inequality operator. 
		///
		/// \author	Fabrice Lamarche, University Of Rennes 1
		/// \date	23/01/2010
		///
		/// \param [in,out]	v	the v. 
		///
		/// \return	true if the parameters are not considered equivalent. 
		////////////////////////////////////////////////////////////////////////////////////////////////////
		bool operator!=(Vector<Float, dimensions> const & v) const
		{ return !((*this)==v) ; }

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	Vector<Float, dimensions> normalized() const
		///
		/// \brief	Return the associated normalized vector. 
		///
		/// \author	Fabrice Lamarche, University Of Rennes 1
		/// \date	23/01/2010
		///
		/// \return	. 
		////////////////////////////////////////////////////////////////////////////////////////////////////
		Vector<Float, dimensions> normalized() const
		{
			return (*this)/norm() ;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	int size() const
		///
		/// \brief	Returns the number of dimensions of this object. 
		///
		/// \author	Fabrice Lamarche, University Of Rennes 1
		/// \date	23/01/2010
		///
		/// \return	. 
		////////////////////////////////////////////////////////////////////////////////////////////////////
		int size() const
		{ return dimensions ; }

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	iterator begin()
		///
		/// \brief	Begin iterator. 
		///
		/// \author	Fabrice Lamarche, University Of Rennes 1
		/// \date	23/01/2010
		///
		/// \return	. 
		////////////////////////////////////////////////////////////////////////////////////////////////////
		iterator begin() 
		{ return m_vector ; }

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	iterator end()
		///
		/// \brief End iterator. 
		///
		/// \author	Fabrice Lamarche, University Of Rennes 1
		/// \date	23/01/2010
		///
		/// \return	. 
		////////////////////////////////////////////////////////////////////////////////////////////////////
		iterator end()
		{ return m_vector+dimensions ; }

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	const_iterator begin() const
		///
		/// \brief	Begin const iterator. 
		///
		/// \author	Fabrice Lamarche, University Of Rennes 1
		/// \date	23/01/2010
		///
		/// \return	. 
		////////////////////////////////////////////////////////////////////////////////////////////////////
		const_iterator begin() const 
		{ return &(m_vector[0]) ; }

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	const_iterator end() const
		///
		/// \brief	End const iterator. 
		///
		/// \author	Fabrice Lamarche, University Of Rennes 1
		/// \date	23/01/2010
		///
		/// \return	. 
		////////////////////////////////////////////////////////////////////////////////////////////////////
		const_iterator end() const
		{ return m_vector+dimensions ; }

		///////////////////////////////////////////////////////////////////////////////////
		/// \brief Computes absolute value of all vector components
		/// 
		/// \return 
		/// 
		/// \author F. Lamarche, University of Rennes 1.
		///////////////////////////////////////////////////////////////////////////////////
		Vector<Float, dimensions> simdAbs() const
		{
			Vector<Float, dimensions> result ;
			for(int cpt=0 ; cpt<dimensions ; ++cpt)
			{
				result.m_vector[cpt] = fabs(m_vector[cpt]) ;
			}
			return result ;
		}

		///////////////////////////////////////////////////////////////////////////////////
		/// \brief Computes the minimum value component by component
		/// 
		/// \param v
		/// \return 
		/// 
		/// \author F. Lamarche, University of Rennes 1.
		///////////////////////////////////////////////////////////////////////////////////
		Vector<Float, dimensions> simdMin(Vector<Float, dimensions> const & v) const
		{
			Vector<Float, dimensions> result ;
			for(int cpt=0 ; cpt<dimensions ; ++cpt)
			{
				result.m_vector[cpt] = ::std::min(m_vector[cpt], v.m_vector[cpt]) ;
			}
			return result ;
		}

		///////////////////////////////////////////////////////////////////////////////////
		/// \brief Computes the maximum value component by component
		/// 
		/// \param v
		/// \return 
		/// 
		/// \author F. Lamarche, University of Rennes 1.
		///////////////////////////////////////////////////////////////////////////////////
		Vector<Float, dimensions> simdMax(Vector<Float, dimensions> const & v) const
		{
			Vector<Float, dimensions> result ;
			for(int cpt=0 ; cpt<dimensions ; ++cpt)
			{
				result.m_vector[cpt] = ::std::max(m_vector[cpt], v.m_vector[cpt]) ;
			}
			return result ;
		}

		Vector<Float, dimensions> zeroIfNegativeCoordinate() const
		{
			bool zero = false ;
			for(int cpt=0 ; cpt<dimensions ; ++cpt)
			{
				if(m_vector[cpt]<0.0) { return Vector<Float,dimensions>(Float(0.0)) ; }
			}
			return *this ;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	Math::Vector<Float, dimensions+1> :::insert(int index, Float const & value) const
		///
		/// \brief	Creates a new vector by inserting a coordinate at position index.
		///
		/// \author	Fabrice Lamarche, university of Rennes 1
		/// \param	index	Zero-based index of the coordinate to insert.
		/// \param	value	The inserted value.
		///
		/// \return	The new vector.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		Math::Vector<Float, dimensions+1> insert(int index, Float const & value) const
		{
			Math::Vector<Float, dimensions+1> result ;
			for(int cpt=0 ; cpt<index ; ++cpt)
			{
				result[cpt] = (*this)[cpt] ;
			}
			result[index] = value ;
			for(int cpt=index ; cpt<dimensions ; ++cpt)
			{
				result[cpt+1] = (*this)[cpt] ;
			}
			return result ;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	Math::Vector<Float, dimensions-1> :::remove(int index) const
		///
		/// \brief	Creates a new vector by removing the given coordinate.
		///
		/// \author	Fabrice Lamarche, university of Rennes 1
		/// \param	index	The index of the removed coordinate.
		///
		/// \return	The new vector.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		Math::Vector<Float, dimensions-1> remove(int index) const
		{
			Math::Vector<Float, dimensions-1> result ;
			for(int cpt=0 ; cpt<index ; ++cpt)
			{
				result[cpt] = (*this)[cpt] ;
			}
			for(int cpt=index+1 ; cpt<dimensions ; ++cpt)
			{
				result[cpt-1] = (*this)[cpt] ;
			}
			return result ;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	Math::Vector<Float, dimensions+1> :::pushBack(Float const & value) const
		///
		/// \brief	Creates a new vector by pushing the given value at the end of this vector.
		///
		/// \author	Fabrice Lamarche, university of Rennes 1
		/// \param	value	The added value.
		///
		/// \return	The new vector.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		Math::Vector<Float, dimensions+1> pushBack(Float const & value) const
		{
			return insert(dimensions, value) ;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	Math::Vector<Float, dimensions+1> :::pushFront(Float const & value) const
		///
		/// \brief	Creates a new vector by pushing a value at the beginning of this vector.
		///
		/// \author	Fabrice Lamarche, university of Rennes 1
		/// \param	value	The value to push.
		///
		/// \return The new vector.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		Math::Vector<Float, dimensions+1> pushFront(Float const & value) const
		{
			return insert(0, value) ;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	Math::Vector<Float, dimensions-1> :::popBack() const
		///
		/// \brief	Creates a new vector by removing the last coordinate.
		///
		/// \author	Fabrice Lamarche, university of Rennes 1
		/// \return	The new vector.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		Math::Vector<Float, dimensions-1> popBack() const
		{
			return remove(dimensions-1) ;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	Math::Vector<Float, dimensions-1> :::popFront() const
		///
		/// \brief	Creates a new vector by removing the first coordinate.
		///
		/// \author	Fabrice Lamarche, university of Rennes 1
		/// \return	The new vector.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		Math::Vector<Float, dimensions-1> popFront() const
		{
			return remove(0) ;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	Math::Vector<Float, dimensions> simdMul(Math::Vector<Float, dimensions> const & v) const
		///
		/// \brief	Computes the multiplication coordinate by coordinate of two vectors.
		///
		/// \author	Fabrice Lamarche, university of Rennes 1
		/// \return	The new vector.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		Math::Vector<Float, dimensions> simdMul(Math::Vector<Float, dimensions> const & v) const
		{
			Math::Vector<Float, dimensions> result;
			for (int cpt = 0; cpt < dimensions; ++cpt)
			{
				result[cpt] = (*this)[cpt] * v[cpt];
			}
			return result;
		}

		template<class T, class Q>
		Math::Vector<T, dimensions> simdDiv(Math::Vector<Q, dimensions> const & v) const
		{
			Math::Vector<Float, dimensions> result;
			for (int cpt = 0; cpt < dimensions; ++cpt)
			{
				result[cpt] = (*this)[cpt] / v[cpt];
			}
			return result;
		}



		static Vector<Float, 3> make_sphere(Float inclination, Float azimuth)
		{
			return Vector<Float, 3>(sin(inclination) * cos(azimuth), sin(inclination) * sin(azimuth), cos(inclination));
		}


		//Assumes that this is normal
		Vector reflect(Vector const& dir)const
		{
			return (*this) * (dir * (*this)) * 2 - dir;
		}


	} ;

	template <class Float>
	Vector<Float, 3> operator^(const Vector<Float, 3> & v1, const Vector<Float, 3> & v2)
	{
		Vector<Float,3> result ;
		result[0] = v1[1]*v2[2]-v1[2]*v2[1] ;
		result[1] = v1[2]*v2[0]-v1[0]*v2[2] ;
		result[2] = v1[0]*v2[1]-v1[1]*v2[0] ;
		return result ;
	}


	template <class Float>
	Vector<Float, 2> makeVector(Float const & c1, Float const & c2)
	{
		Float table[] = {c1, c2} ;
		return Math::Vector<Float, 2>(table) ;
	}


	template <class Float>
	Vector<Float, 3> makeVector(Float const & c1, Float const & c2, Float const & c3)
	{
		Float table[] = {c1, c2, c3} ;
		return Math::Vector<Float, 3>(table) ;
	}


	template <class Float>
	Vector<Float, 4> makeVector(Float const & c1, Float const & c2, Float const & c3, Float const & c4)
	{
		Float table[] = {c1, c2, c3, c4} ;
		return Math::Vector<Float, 4>(table) ;
	}


	template <class Float>
	Vector<Float, 5> makeVector(Float const & c1, Float const & c2, Float const & c3, Float const & c4, Float const & c5)
	{
		Float table[] = {c1, c2, c3, c4, c5} ;
		return Math::Vector<Float, 5>(table) ;
	}


	template <class Float, int dimensions>
	inline std::ostream & operator<< (std::ostream & out, Vector<Float, dimensions> const & v)
	{
		out << "{";
		for(int cpt=0 ; cpt<dimensions ; cpt++)
		{
			out << v[cpt];
			if (cpt != dimensions-1)
				out << ", ";
		}
		out << "}";
		return out ;
	}

} 


#endif
