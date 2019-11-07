#pragma once

#include <Analysis.h>

class CrankNicolson: public Analysis
{
protected:
	arma::mat A;
	virtual MeshData CreateCNConvDiff(double CFL);
	virtual void CalcPhi(int i, int j, double conv, double diff,
			const MeshData &data);
	virtual bool SolvePhi();
public:
	CrankNicolson(int nx, int ny) :
			Analysis(nx, ny), A(nx * ny, nx * ny, arma::fill::zeros)
	{
	}
	virtual ~CrankNicolson()
	{
	}
};

