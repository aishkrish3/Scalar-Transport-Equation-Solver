#pragma once

#include <Mesh.h>
#include <Point3D.h>
#include <armadillo>

class Analysis
{
protected:
	arma::mat phi;
	arma::mat phi_old;
	arma::vec phi_mid;
	Mesh mesh;

	int get_nx() const
	{
		return phi.n_rows;
	}
	int get_ny() const
	{
		return phi.n_cols;
	}
	virtual MeshData CreateCNConvDiff(double CFL) = 0;
	virtual void CalcPhi(int i, int j, double conv, double diff,
			const MeshData &data) = 0;
	virtual bool SolvePhi() = 0;
public:
	Analysis(int nx, int ny) :
			phi(nx, ny, arma::fill::zeros), phi_old(nx, ny, arma::fill::zeros), phi_mid(
					nx * ny, arma::fill::zeros), mesh(nx, ny)
	{
	}

	virtual ~Analysis()
	{
	}
	std::vector<std::vector<Point3D> > Analyze(double CFL);

};

