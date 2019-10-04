#ifndef _Math_Finite_H
#define _Math_Finite_H

#include <limits>

namespace Math
{
	////////////////////////////////////////////////////////////////////////////////////////////////////
	/// \fn	inline bool is_finite(double value)
	///
	/// \brief	Query if 'value' is finite.
	///
	/// \author	F. Lamarche, Université de Rennes 1
	/// \date	02/03/2016
	///
	/// \param	value	The value.
	///
	/// \return	true if finite, false if not.
	////////////////////////////////////////////////////////////////////////////////////////////////////
	inline bool is_finite(double value)
	{
		return fabsf(value)!=::std::numeric_limits<double>::infinity() ;
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////
	/// \fn	inline bool is_indeterminate(double value)
	///
	/// \brief	Query if 'value' is indeterminate.
	///
	/// \author	F. Lamarche, Université de Rennes 1
	/// \date	02/03/2016
	///
	/// \param	value	The value.
	///
	/// \return	true if indeterminate, false if not.
	////////////////////////////////////////////////////////////////////////////////////////////////////
	inline bool is_indeterminate(double value)
	{
		return !(value==value) ;
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////
	/// \fn	inline bool is_valid(double value)
	///
	/// \brief	Query if 'value' is valid.
	///
	/// \author	F. Lamarche, Université de Rennes 1
	/// \date	02/03/2016
	///
	/// \param	value	The value.
	///
	/// \return	true if valid, false if not.
	////////////////////////////////////////////////////////////////////////////////////////////////////
	inline bool is_valid(double value)
	{
		return is_finite(value) && !is_indeterminate(value) ;
	}
}

#endif