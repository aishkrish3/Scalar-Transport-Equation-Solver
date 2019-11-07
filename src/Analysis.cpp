#include <Analysis.h>

std::vector<std::vector<Point3D> > Analysis::Analyze(double CFL)
{
	std::vector<std::vector<Point3D> > final_result;

	const int nx = get_nx(), ny = get_ny();
	double dxi = mesh.get_dxi(); /*inverse of dx*/
	double dyi = mesh.get_dyi(); /*inverse of dy*/
	double dxi_2 = mesh.get_dxi_2(); /*inverse of dx^2*/
	double dyi_2 = mesh.get_dyi_2(); /*inverse of dy^2*/

	std::cout << "Calculating Convection and Diffusion terms...\n";
	auto data = CreateCNConvDiff(CFL);

	const double output_freq = Mesh::tf / 4.;
	double otime = output_freq;

	auto check_range = [](int &ii, int size)
	{
		if (ii < 0)
		{	ii += size;}
		if (ii >= size)
		{	ii -= size;}
	};

	std::vector<Point3D> phi_result;

	for (auto &x : mesh.getXm())
	{
		for (auto &y : mesh.getYm())
		{
			phi_result.push_back(Point3D(x, y, 0));
		}
	}

	auto add_result = [&]()
	{
		int k = 0;
		for (int i = 0; i < nx; ++i)
		{
			for (int j = 0; j < ny; ++j)
			{
				phi_result[k++].z = phi(i, j);
			}
		}

		final_result.push_back(phi_result);
	};

	/*Looping through time*/

	for (double t = 0; t <= Mesh::tf; t += data.dt)
	{
		std::cout << "Analyzing time step " << t << "s of " << Mesh::tf
				<< "s...\n";
		/*Looping through space*/
		for (int i = 0; i < nx; i++)
		{
			for (int j = 0; j < ny; j++)
			{
				int ip1 = i + 1;
				int ip2 = i + 2;
				int im1 = i - 1;
				int im2 = i - 2;
				check_range(ip1, nx);
				check_range(ip2, nx);
				check_range(im1, nx);
				check_range(im2, nx);

				int jp1 = j + 1;
				int jp2 = j + 2;
				int jm1 = j - 1;
				int jm2 = j - 2;

				check_range(jp1, ny);
				check_range(jp2, ny);
				check_range(jm1, ny);
				check_range(jm2, ny);
				/*calculate diffusion*/
				double diff = 0;

				diff = Mesh::alpha * dxi_2
						* (phi_old(im1, j) - 2. * phi_old(i, j)
								+ phi_old(ip1, j));
				diff += Mesh::alpha * dyi_2
						* (phi_old(i, jm1) - 2. * phi_old(i, j)
								+ phi_old(i, jp1));
				/*Face velocities*/
				double ue = data.u[i + 1][j];
				double uw = data.u[i][j];
				double un = data.v[i][j + 1];
				double us = data.v[i][j];

				double phi_e, phi_w, phi_n, phi_s;
				if (ue > 0)
					phi_e = (-phi_old(im1, j) + 5. * phi_old(i, j)
							+ 2. * phi_old(ip1, j)) / 6.;
				else
					phi_e = (2. * phi_old(i, j) + 5. * phi_old(ip1, j)
							- phi_old(ip2, j)) / 6.;
				if (uw > 0)
					phi_w = (-phi_old(im2, j) + 5. * phi_old(im1, j)
							+ 2. * phi_old(i, j)) / 6.;
				else
					phi_w = (2. * phi_old(im1, j) + 5. * phi_old(i, j)
							- phi_old(ip1, j)) / 6.;
				if (un > 0)
					phi_n = (-phi_old(i, jm1) + 5. * phi_old(i, j)
							+ 2. * phi_old(i, jp1)) / 6.;
				else
					phi_n = (2. * phi_old(i, j) + 5. * phi_old(i, jp1)
							- phi_old(i, jp2)) / 6.;
				if (us > 0)
					phi_s = (-phi_old(i, jm2) + 5. * phi_old(i, jm1)
							+ 2. * phi_old(i, j)) / 6.;
				else
					phi_s = (2. * phi_old(i, jm1) + 5. * phi_old(i, j)
							- phi_old(i, jp1)) / 6.;

				/*calculate convection*/
				double conv = -dxi * (ue * phi_e - uw * phi_w);
				conv -= dyi * (un * phi_n - us * phi_s);
				CalcPhi(i, j, conv, diff, data);
			}
		}

		/*solve phi matrix using Armadillo*/
		if (!SolvePhi())
		{
			std::cout << "\nFailed to solve for Phi!\n";
			break;
		}
		phi_old = phi;

		if (t > otime)
		{
			add_result();
			otime += output_freq;
		}
	}

	add_result();
	std::cout << "\nDone!\n";
	return final_result;
}
