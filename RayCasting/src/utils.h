#pragma once


#include <Geometry\Triangle.h>
#include <Geometry\BoundingBox.h>
#include <algorithm>
#include <ostream>
#include <cmath>
#include <string>
#include <sstream>
#include <chrono>

///////////////////////////////
// 
//////////////////////////////
inline Geometry::BoundingBox get_bounding_box(Geometry::Triangle const& tri)
{
	Math::Vector3f min(tri.vertex(0)), max(tri.vertex(0));
	for (unsigned int i = 0; i < 3; ++i)
	{
		min[i] = ::std::min(min[i], tri.vertex(1)[i]);
		min[i] = ::std::min(min[i], tri.vertex(2)[i]);

		max[i] = ::std::max(max[i], tri.vertex(1)[i]);
		max[i] = ::std::max(max[i], tri.vertex(2)[i]);
	}
	return Geometry::BoundingBox(min, max);
}


inline double my_atan(double y, double x)
{
	return atan2(y, x) + Math::pi;
}

inline Math::Vector3f spherical_coordinate(Math::Vector3f const& vec)
{
	Math::Vector3f res;

	res[0] = vec.norm();

	res[1] = acos(vec[2] / res[0]);//theta: inclination

	res[2] = my_atan(vec[1], vec[0]);//phi: azimuth

	return res;
}



inline std::string & operator+=(std::string & str, int i)
{
	std::stringstream sstr;
	sstr << i;
	str += sstr.str();
	return str;
}

inline std::string & operator+=(std::string & str, double i)
{
	std::stringstream sstr;
	sstr << i;
	str += sstr.str();
	if (i == (int)i)
	{
		str += ".0";
	}
	return str;
}



inline std::string percent(int current, int total)
{
	std::string res;

	double p = (current * 1000.0) / total;

	p = floor(p);
	p /= 10.0;

	res += p;

	res += '%';

	return res;
}


inline std::string progession_bar(int current, int total, unsigned int samples = 20)
{
	std::string res = "[";
	int n = (current * samples) / total;
	res += std::string(n, '=');
	res += std::string(samples - n, ' ');
	res += "]";
	res += ": " + percent(current, total);
	return res;
}

template <unsigned int N, class T>
inline Math::Vector<T, N> sqrt(Math::Vector<T, N> const& vec)
{
	Math::Vector<T, N> res;
	for (unsigned int i = 0; i < N; ++i)
	{
		res[i] = sqrt(vec[i]);
	}
	return res;
}


template <class Float>
inline bool valid_uv(Math::Vector<Float, 2> const& uv)
{
	return !(uv[0] < 0 || uv[0] > 1 || uv[1] < 0 || uv[1] > 1);
}


unsigned long long nano()
{
	auto ima = std::chrono::high_resolution_clock::now();
	return ima.time_since_epoch().count();
}


double rad(double deg)
{
	return deg * Math::pi / 180.0;
}
