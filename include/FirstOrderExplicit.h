#pragma once

#include <Analysis.h>

class FirstOrderExplicit: public Analysis
{
protected:
	virtual MeshData CreateCNConvDiff(double CFL);
	virtual void CalcPhi(int i, int j, double conv, double diff,
			const MeshData &data);
	virtual bool SolvePhi();
public:
	FirstOrderExplicit(int nx, int ny) :
			Analysis(nx, ny)
	{
	}
	virtual ~FirstOrderExplicit()
	{
	}
};
