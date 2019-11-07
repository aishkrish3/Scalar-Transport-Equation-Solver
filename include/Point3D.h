#pragma once
#include <gnuplot-iostream.h>
struct Point3D
{
	Point3D();
	Point3D(double _x, double _y, double _z);

	double x, y, z;
};

// Tells gnuplot-iostream how to print objects of class Point3D.
namespace gnuplotio
{
template <>
struct BinfmtSender<Point3D>
{
	static void send(std::ostream &stream)
	{
		BinfmtSender<double>::send(stream);
		BinfmtSender<double>::send(stream);
		BinfmtSender<double>::send(stream);
	}
};

template <>
struct BinarySender<Point3D>
{
	static void send(std::ostream &stream, const Point3D &v)
	{
		BinarySender<double>::send(stream, v.x);
		BinarySender<double>::send(stream, v.y);
		BinarySender<double>::send(stream, v.z);
	}
};

template<>
struct TextSender<Point3D> {
	static void send(std::ostream &stream, const Point3D &v) {
		TextSender<double>::send(stream, v.x);
		stream << " ";
		TextSender<double>::send(stream, v.y);
		stream << " ";
		TextSender<double>::send(stream, v.z);
	}
};
}

