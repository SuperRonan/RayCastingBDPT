#ifndef _Geometry_RGBColor
#define _Geometry_RGBColor

#include <iostream>
#include <algorithm>
#include <cassert>
#include <ImfRgba.h>

namespace Geometry
{
	////////////////////////////////////////////////////////////////////////////////////////////////////
	/// \class	RGBColor
	///
	/// \brief	A Red / Green / Blue color. Each component should be in range [0;1]. If a component 
	/// 		value is greater than 1, rendering should be able to handle HDR.
	///
	/// \author	F. Lamarche, Université de Rennes 1
	/// \date	04/12/2013
	////////////////////////////////////////////////////////////////////////////////////////////////////
	class RGBColor
	{
	protected:
		/// \brief	The red(0), green(1) and blue (2) components.
		double m_color[3] ;

	public:
		//void validateColor()
		//{
		//	//if (m_color[0] < 0.f || m_color[1] < 0.f || m_color[2] < 0.f)
		//	//{
		//	//	for (int cpt = 0; cpt < 3; ++cpt)
		//	//	{
		//	//		//::std::cout<<m_color[cpt]<<"  "<<::std::flush;
		//	//		m_color[cpt] = std::max(m_color[cpt], 0.0);
		//	//	}
		//	//	//char c;
		//	//	//::std::cin >> c;
		//	//	::std::cerr << "Invalid color!" << ::std::endl;
		//	//	//char * ptr = nullptr;
		//	//	//*ptr = 0;
		//	//	//assert(false);
		//	//}
		//}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	RGBColor::RGBColor(double R=0, double G=0, double B=0)
		///
		/// \brief	Constructor.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	04/12/2013
		///
		/// \param	R	Red component.
		/// \param	G	Green component.
		/// \param	B	Blue component.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		RGBColor(double R, double G, double B)
		{
			m_color[0] = R ;
			m_color[1] = G ;
			m_color[2] = B ;
			//validateColor();
		}

		RGBColor(double g=0)
		{
			m_color[0] = g;
			m_color[1] = g;
			m_color[2] = g;
			//validateColor();
		}

		/// <summary>
		/// Determines whether this color is black.
		/// </summary>
		/// <returns>
		///   <c>true</c> if this color is black; otherwise, <c>false</c>.
		/// </returns>
		bool isBlack() const
		{
			return m_color[0] == 0.0f && m_color[1] == 0.0f && m_color[2] == 0.0f;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	void set(double r, double g, double b)
		///
		/// \brief	Sets the color.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	04/12/2013
		///
		/// \param	r	Red component.
		/// \param	g	Green component.
		/// \param	b	Blue component.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		void set(double r, double g, double b)
		{
			m_color[0] = r;
			m_color[1] = g;
			m_color[2] = b;

		}

	
		/// <summary>
		/// Sets the color with the provided array of doubles (at least 3 doubles long : r,g,b components).
		/// </summary>
		/// <param name="rgb">The RGB components.</param>
		void set(const double * rgb)
		{
			m_color[0] = rgb[0];
			m_color[1] = rgb[1];
			m_color[2] = rgb[2];

		}

		/// <summary>
		/// Sets the color with the provided array of floats (at least 3 floats long : r,g,b components).
		/// </summary>
		/// <param name="rgb">The RGB components.</param>
		void set(const float * rgb)
		{
			m_color[0] = rgb[0];
			m_color[1] = rgb[1];
			m_color[2] = rgb[2];

		}

		/// <summary>
		/// Returns the grey value (average of red, green and blue).
		/// </summary>
		/// <returns></returns>
		double grey() const
		{
			return (m_color[0] + m_color[1] + m_color[2]) / 3.0;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	RGBColor RGBColor::operator+ (RGBColor const & c) const
		///
		/// \brief	Addition operator.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	04/12/2013
		///
		/// \param	c	The other color.
		///
		/// \return	The sum of this color with the provided one.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		RGBColor operator+ (RGBColor const & c) const
		{
			return RGBColor(c.m_color[0]+m_color[0], c.m_color[1]+m_color[1], c.m_color[2]+m_color[2]) ;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	RGBColor RGBColor::operator* (RGBColor const & c) const
		///
		/// \brief	Multiplication operator. Multiplication is achieved by multiplying the two red components,
		/// 		the two green components and the two blue components.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	04/12/2013
		///
		/// \param	c	The other rgb color.
		///
		/// \return	The result of the operation.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		RGBColor operator* (RGBColor const & c) const
		{
			return RGBColor(c.m_color[0]*m_color[0], c.m_color[1]*m_color[1], c.m_color[2]*m_color[2]) ;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	RGBColor RGBColor::operator* (double v) const
		///
		/// \brief	Multiplies a color by a scalar value.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	04/12/2013
		///
		/// \param	v	The scalar.
		///
		/// \return	The result of the operation.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		RGBColor operator* (double v) const
		{
			return RGBColor(m_color[0]*v, m_color[1]*v, m_color[2]*v) ;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	RGBColor RGBColor::operator/ (double v) const
		///
		/// \brief	Devides a color by a scalar value.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	04/12/2013
		///
		/// \param	v	The scalar.
		///
		/// \return	The result of the operation.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		RGBColor operator/ (double v) const
		{
			return RGBColor(m_color[0]/v, m_color[1]/v, m_color[2]/v) ;			
		}


		RGBColor& operator+=(RGBColor const& color)
		{
			m_color[0] += color[0];
			m_color[1] += color[1];
			m_color[2] += color[2];
			return *this;
		}

		RGBColor& operator*=(RGBColor const& color)
		{
			m_color[0] *= color[0];
			m_color[1] *= color[1];
			m_color[2] *= color[2];
			return *this;
		}

		RGBColor& operator-=(RGBColor const& other)
		{
			m_color[0] -= other[0];
			m_color[1] -= other[1];
			m_color[2] -= other[2];
			return *this;
		}

		RGBColor operator-(RGBColor const& other)const
		{
			RGBColor res = *this;
			return res -= other;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	double RGBColor::operator[] (int c) const
		///
		/// \brief	Access to each color component. 0: red, 1: green, 2: blue.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	04/12/2013
		///
		/// \param	c	The index of the component.
		///
		/// \return	The value of the component.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		double operator[] (int c) const
		{
			return m_color[c] ;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	double & RGBColor::operator[] (int c)
		///
		/// \brief	Access to each color component. 0: red, 1: green, 2: blue.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	04/12/2013
		///
		/// \param	c	The index of the component.
		///
		/// \return	The value of the component.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		double & operator[] (int c)
		{
			return m_color[c] ;
		}

		/// <summary>
		/// Equality comparator
		/// </summary>
		/// <param name="color"> The other color </param>
		/// <returns></returns>
		bool operator==(RGBColor const & color) const
		{
			return m_color[0]==color[0] && m_color[1]==color[1] && m_color[2]==color[2] ;
		}

		/// <summary>
		/// Inequality comparator
		/// </summary>
		/// <param name="color"> The other color </param>
		/// <returns></returns>
		bool operator!=(RGBColor const & color) const
		{
			return !((*this)==color) ;
		}

		/// <summary>
		/// Lexicographical comparison between colors.
		/// </summary>
		/// <param name="color"> The color </param>
		/// <returns></returns>
		bool operator<(const RGBColor & color) const
		{
			return std::lexicographical_compare(m_color, m_color + 3, color.m_color, color.m_color + 3);
		}
			

		bool any_negative()const
		{
			//use std::algo
			return m_color[0] < 0 || m_color[1] < 0 || m_color[2] < 0;
		}

		double brightness()const
		{
			return m_color[0] * 0.3 + m_color[1] * 0.59 + m_color[2] * 0.11;
		}

		bool anythingWrong()const
		{
			for (int i = 0; i < 3; ++i)
			{
				if (std::isnan(m_color[i]) || std::isinf(m_color[i]))
				{
					return true;
				}
			}
			return false;
		}





		Imf::Rgba convertToImf()const
		{
			Imf::Rgba res;
			res.r = m_color[0];
			res.g = m_color[1];
			res.b = m_color[2];
			res.a = 0;
			return res;
		}



	} ;

	/// <summary>
	/// Outputs a color in a stream.
	/// </summary>
	/// <param name="out"></param>
	/// <param name="color"></param>
	/// <returns></returns>
	inline ::std::ostream & operator << (::std::ostream & out, RGBColor const & color)
	{
		out << "(" << color[0] << "," << color[1] << "," << color[2] << ")";
		return out;
	}
}

#endif
