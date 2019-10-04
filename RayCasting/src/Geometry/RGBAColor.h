#pragma once


#include "RGBColor.h"

namespace Geometry
{
	class RGBAColor
	{
	protected:
		RGBColor m_color;
		double m_alpha;

	public:

		RGBAColor(double gray = 0, double a=1) :
			m_color(gray, gray, gray), m_alpha(a)
		{}

		RGBAColor(double r, double g, double b, double a = 1) :
			m_color(r, g, b), m_alpha(a)
		{}

		RGBAColor(RGBColor const & c, double a = 1):
			m_color(c), m_alpha(a)
		{}


		double alpha()const
		{
			return m_alpha;
		}

		double & alpha()
		{
			return m_alpha;
		}


		bool isBlack() const
		{
			return m_color.isBlack();
		}

		bool  is_transparent() const
		{
			return m_alpha < 1.0;
		}


		

		void set(double r, double g, double b, double a)
		{
			m_color.set(r, g, b);
			m_alpha = a;
		}

		void set(double * rgb, double a)
		{
			m_color.set(rgb);
			m_alpha = a;
		}

		void set(float * rgb, double a)
		{
			m_color.set(rgb);
			m_alpha = a;
		}

		void set(double r, double g, double b)
		{
			m_color.set(r, g, b);
		}

		void set(double * rgb)
		{
			m_color.set(rgb);
		}

		void set(float * rgb)
		{
			m_color.set(rgb);
		}

		double grey() const
		{
			return m_color.grey();
		}


		RGBAColor operator+ (RGBAColor const& c)const
		{
			double new_alpha = m_alpha + c.m_alpha * (1 - m_alpha);
			RGBColor res = (m_color * m_alpha + c.m_color * c.m_alpha * (1 - m_alpha)) / new_alpha;
			return RGBAColor(res, new_alpha);
		}

		//This operator is not comutative anymore
		RGBAColor operator*(RGBAColor const& c)const
		{
			//TODO
			return RGBAColor(m_color * c.m_color, m_alpha);
		}

		RGBAColor operator*(double v) const
		{
			return RGBAColor(m_color * v, m_alpha);
		}

		RGBAColor operator/(double v) const
		{
			return RGBAColor(m_color / v, m_alpha);
		}

		double operator[] (int c) const
		{
			return m_color[c];
		}

		double & operator[] (int c)
		{
			return m_color[c];
		}




		bool operator==(RGBAColor const & c)
		{
			return m_color == c.m_color && m_alpha == c.m_alpha;
		}

		bool operator!=(RGBAColor const & c)
		{
			return m_color != c.m_color || m_alpha != c.m_alpha;
		}

		bool operator<(const RGBAColor & color) const
		{
			return m_color < color.m_color;
		}





	};

	inline ::std::ostream & operator << (::std::ostream & out, RGBAColor const & color)
	{
		out << "(" << color[0] << "," << color[1] << "," << color[2] << "x" << color.alpha() << ")";
		return out;
	}
}