#include <boost/tuple/tuple.hpp>
#include <boost/array.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/irange.hpp>
#include <boost/bind.hpp>

#include <Mesh.h>
#include <ScalarTransportEquation.h>
#include <CrankNicolson.h>
#include <FirstOrderExplicit.h>

int main()
{
	const int MESH_SIZE = 32;
	const double CFL = 0.8;  //Courant–Friedrichs–Lewy condition
	int n;
	std::cout
			<< "Enter 1 for Crank Nicolson; 2 for First Order Explicit for time integration"
			<< std::endl;
	std::cin >> n;
	auto analysis =
			(n == 1) ?
					dynamic_cast<Analysis*>(new CrankNicolson(MESH_SIZE,
							MESH_SIZE)) :
					dynamic_cast<Analysis*>(new FirstOrderExplicit(MESH_SIZE,
							MESH_SIZE));

	auto result = analysis->Analyze(CFL);
	delete analysis;
	/* Plot scalar with time*/

	Gnuplot gp;

	gp << "set key top" << std::endl;
	gp << "set view map" << std::endl;
	gp << "set xrange [0:1]" << std::endl;
	gp << "set yrange [0:1]" << std::endl;
	gp << "set cbrange [-0.15:0.15]" << std::endl;
	gp << "set contour base" << std::endl;
	gp << "set cntrparam levels incr -0.15,0.02,0.15 " << std::endl;
	gp << "set multiplot layout 2,2 rowsfirst title 'Analysis results'"
			<< std::endl;

	double output_freq = Mesh::tf / 4.;
	double i = 1;
	for (auto plot : result)
	{
		gp << "set title 't = " << std::fixed << std::setprecision(2)
				<< (i++) * output_freq << "s'" << std::endl;
		gp << "plot '-' w image" << std::endl;
		gp.send1d(plot);
	}
}
