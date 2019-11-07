#pragma once

#include <vector>
struct MeshData
{
	MeshData();
	double dt;
	std::vector<std::vector<double> > u;
	std::vector<std::vector<double> > v;
	std::vector<std::vector<double> > omega;
};

class Mesh
{
public:
	Mesh(int nx, int ny);
	virtual ~Mesh();

	double get_dx() const;
	double get_dy() const;

	double get_dxi() const;
	double get_dyi() const;

	double get_dxi_2() const;
	double get_dyi_2() const;

	MeshData get_data(const double CFL) const;

	const std::vector<double>& getXm() const;

	const std::vector<double>& getYm() const;

	static constexpr double alpha = 1e-2;
	static constexpr double tf = 4.;
private:
	std::vector<double> x, y, xm, ym;

};
