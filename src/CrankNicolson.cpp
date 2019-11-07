#include "CrankNicolson.h"
#include <math.h>

MeshData CrankNicolson::CreateCNConvDiff(double CFL)
{
	const int nx = get_nx(), ny = get_ny();
	double dxi_2 = mesh.get_dxi_2();
	double dxi = mesh.get_dxi();
	auto data = mesh.get_data(CFL);

	/*Generate A matrix only for Crank Nicolson using lexicographic ordering*/
	A = arma::mat(nx * ny, nx * ny, arma::fill::zeros);

	auto check_range = [](int &ii, int size)
	{
		if (ii < 0)
		{	ii += size;}
		if (ii >= size)
		{	ii -= size;}
	};

	for (int i = 0; i < nx; i++)
	{
		for (int j = 0; j < ny; j++)
		{
			int n = i + (nx * j);

			A(n, n) = -4. * dxi_2;
			auto modify_A = [&](int ii, int jj)
			{
				check_range(ii, nx);
				check_range(jj, ny);
				int nn = ii + (nx * jj);
				A(n,nn) = dxi_2;
			};

			modify_A(i + 1, j);
			modify_A(i - 1, j);
			modify_A(i, j + 1);
			modify_A(i, j - 1);
		}
	}

	A *= -Mesh::alpha;
	for (int i = 0; i < nx; i++)
	{
		for (int j = 0; j < ny; j++)
		{
			int n = i + (nx * j);
			double ue = data.u[i + 1][j];
			double uw = data.u[i][j];
			double un = data.v[i][j + 1];
			double us = data.v[i][j];

			A(n, n) += 5. / 6. * std::max(0., ue) * dxi
					+ 2. / 6. * std::min(0., ue) * dxi
					- 2. / 6. * std::max(0., uw) * dxi
					- 5. / 6. * std::min(0., uw) * dxi
					+ 5. / 6. * std::max(0., un) * dxi
					+ 2. / 6. * std::min(0., un) * dxi
					- 2. / 6. * std::max(0., us) * dxi
					- 5. / 6. * std::min(0., us) * dxi;

			int ii = i - 1;
			int jj = j;
			check_range(ii, nx);

			int nn = ii + (nx * jj);

			A(n, nn) += -1. / 6. * std::max(0., ue) * dxi
					- 5. / 6. * std::max(0., uw) * dxi
					- 2. / 6. * std::min(0., uw) * dxi;

			ii = i + 1;
			jj = j;
			check_range(ii, nx);

			nn = ii + (nx * jj);
			A(n, nn) += 2. / 6. * std::max(0., ue) * dxi
					- 5. / 6. * std::min(0., ue) * dxi
					+ 1. / 6. * std::min(0., uw) * dxi;

			ii = i;
			jj = j - 1;
			check_range(jj, ny);
			nn = ii + (nx * jj);

			A(n, nn) += -1. / 6. * std::max(0., un) * dxi
					- 5. / 6. * std::max(0., us) * dxi
					- 2. / 6. * std::min(0., us) * dxi;

			ii = i;
			jj = j + 1;
			check_range(jj, ny);
			nn = ii + (nx * jj);

			A(n, nn) += 2. / 6. * std::max(0., un) * dxi
					+ 5. / 6. * std::min(0., un) * dxi
					+ 1. / 6. * std::min(0., us) * dxi;

			ii = i - 2;
			jj = j;
			check_range(ii, nx);
			nn = ii + (nx * jj);

			A(n, nn) += 1. / 6. * std::max(0., uw) * dxi;

			ii = i + 2;
			jj = j;
			check_range(ii, nx);
			nn = ii + (nx * jj);

			A(n, nn) += -1. / 6. * std::min(0., ue) * dxi;

			ii = i;
			jj = j - 2;
			check_range(jj, ny);
			nn = ii + (nx * jj);

			A(n, nn) += 1. / 6. * std::max(0., us) * dxi;

			ii = i;
			jj = j + 2;
			check_range(jj, ny);
			nn = ii + (nx * jj);

			A(n, nn) += -1. / 6. * std::min(0., un) * dxi;
		}
	}

	A *= data.dt / 2.;

	for (int i = 0; i < nx; i++)
	{
		for (int j = 0; j < ny; j++)
		{
			int n = i + (nx * j);

			A(n, n) += 1;
		}
	}

	return data;
}

void CrankNicolson::CalcPhi(int i, int j, double conv, double diff,
		const MeshData &data)
{
	/*Calculate phi explicitly*/
	const int nx = get_nx();
	int n = i + (nx * j);
	phi_mid(n) = phi_old(i, j)
			+ 0.5 * data.dt * (conv + diff + 2. * data.omega[i][j]);

}

bool CrankNicolson::SolvePhi()
{
	/*solve phi using Armadillo*/
	try
	{
		auto phi1D = arma::solve(A, phi_mid);
		phi = arma::reshape(phi1D, get_nx(), get_ny());
		return true;
	} catch (std::exception &e)
	{
		std::cerr << e.what() << std::endl;
	}

	return false;
}

