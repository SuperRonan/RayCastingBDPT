#pragma once

namespace Geometry
{

	class medium
	{
	protected:
		double m_index;
		RGBColor m_color;

		//maybe add a color, for attenuation...

	public:



		medium(double index = 1.0, RGBColor const& color=RGBColor(1.0, 1.0, 1.0)) :
			m_index(index), m_color(color)
		{}

		medium(medium const& med):
			m_index(med.m_index), m_color(med.m_color)
		{}


		double index()const
		{
			return m_index;
		}

		double & index()
		{
			return m_index;
		}

		const RGBColor & color()const
		{
			return m_color;
		}

		RGBColor & color()
		{
			return m_color;
		}
	};
}