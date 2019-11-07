#include <algorithm>
#include <math.h>
#include "Mesh.h"

constexpr double Mesh::alpha;   //Diffusion coefficient
constexpr double Mesh::tf;		//Final time

class RangeGenerator
{
public:
	RangeGenerator(int size) :
			i(0), size(size)
	{
	}

	double operator()()
	{
		return i++ / size;
	}
private:
	double i;
	const double size;
};

class RangeFunctionGenerator
{
public:
	RangeFunctionGenerator(const std::vector<double> &x) :
			i(0), x(x)
	{
	}

	double operator()()
	{
		return x[i++] + 0.5 * (x[1] - x[0]);
	}
private:
	double i;
	const std::vector<double> &x;
};

MeshData::MeshData() :
		dt(0)
{
}
Mesh::Mesh(int nx, int ny)
{
	/*Generate Mesh*/
	x.resize(nx + 1, 0);
	y.resize(ny + 1, 0);

	xm.resize(nx, 0);
	ym.resize(ny, 0);

	std::generate(x.begin(), x.end(), RangeGenerator(nx));
	std::generate(y.begin(), y.end(), RangeGenerator(ny));

	std::generate(xm.begin(), xm.end(), RangeFunctionGenerator(x));
	std::generate(ym.begin(), ym.end(), RangeFunctionGenerator(y));
}

Mesh::~Mesh()
{
}

double Mesh::get_dx() const
{
	if (xm.size() >= 2)
		return xm[1] - xm[0];
	return 0;
}

double Mesh::get_dy() const
{
	if (ym.size() >= 2)
		return ym[1] - ym[0];
	return 0;
}

double Mesh::get_dxi() const
{
	auto dx = get_dx();

	if (dx == 0)
		return 0;
	return 1. / dx;
}

double Mesh::get_dyi() const
{
	auto dy = get_dy();

	if (dy == 0)
		return 0;
	return 1. / dy;

}

double Mesh::get_dxi_2() const
{
	return pow(get_dxi(), 2.);
}

double Mesh::get_dyi_2() const
{
	return pow(get_dyi(), 2.);
}

const std::vector<double>& Mesh::getXm() const
{
	return xm;
}

const std::vector<double>& Mesh::getYm() const
{
	return ym;
}

MeshData Mesh::get_data(const double CFL) const
{
	MeshData data;
	double maxu = 0.;
	double dt_c, dt_v;

	/*Calculate face velocities*/
	for (auto x_i : x)
	{
		std::vector<double> u1;
		for (auto ym_j : ym)
		{
			u1.push_back(
					1. / 10.
							- pow((sin(M_PI * x_i)), 2)
									* (sin(M_PI * (ym_j - 0.05))
											* cos(M_PI * (ym_j - 0.05))
											- sin(M_PI * (ym_j + 0.05))
													* cos(M_PI * (ym_j + 0.05))));
			maxu = std::max(fabs(u1.back()), maxu);
		}
		data.u.push_back(u1);
	}

	for (auto xm_i : xm)
	{
		std::vector<double> v1;
		for (auto y_j : y)
		{
			v1.push_back(
					sin(M_PI * xm_i) * cos(M_PI * xm_i)
							* (pow((sin(M_PI * (y_j - 0.05))), 2)
									- pow(sin(M_PI * (y_j + 0.05)), 2)));
			maxu = std::max(fabs(v1.back()), maxu);
		}
		data.v.push_back(v1);
	}

	/*Calculate source term at cell center*/
	for (auto xm_i : xm)
	{
		std::vector<double> omega1;
		for (auto ym_j : ym)
		{
			omega1.push_back(
					exp(
							-(pow((xm_i - 0.75), 2) + pow((ym_j - 0.5), 2))
									/ (pow(0.05, 2)))
							- exp(
									-(pow((xm_i - 0.25), 2)
											+ pow((ym_j - 0.5), 2))
											/ (pow(0.05, 2))));
		}
		data.omega.push_back(omega1);
	}

	dt_c = CFL * get_dx() / maxu;
	dt_v = CFL * pow(get_dx(), 2) / (alpha * 4.);
	data.dt = std::min(dt_c, dt_v); //timestep based on CFL
	return data;
}
