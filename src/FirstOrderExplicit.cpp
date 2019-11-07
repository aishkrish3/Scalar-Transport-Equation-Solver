#include "FirstOrderExplicit.h"

MeshData FirstOrderExplicit::CreateCNConvDiff(double CFL)
{
	return mesh.get_data(CFL);
}

bool FirstOrderExplicit::SolvePhi()
{
	return false;
}

void FirstOrderExplicit::CalcPhi(int i, int j, double conv, double diff,
		const MeshData &data)
{
	phi(i, j) = phi_old(i, j) + data.dt * (conv + diff + data.omega[i][j]);
}
