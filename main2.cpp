#include <iostream>
#include <vector>
#include <string>
#include <stdexcept>
#include <iomanip>
#include <fstream>
#include <thread>
#include <print>

#include "mdarray.h"
#include "distrib2d.h"

void write_matrix(const std::string& filename, const array2d<double>& mat)
{
    std::ofstream out(filename);

    if (!out)
    {
        throw std::runtime_error("Failed to open " + filename);
    }

    out << std::setprecision(17);

    const int N = mat.dim0();
    const int M = mat.dim1();

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < M; j++)
        {
            out << mat[i, j] << ' ';
        }

        out << '\n';
    }
}

void write_vec(const std::string& filename, const std::vector<double>& vec)
{
    std::ofstream out(filename);

    if (!out)
    {
        throw std::runtime_error("Failed to open " + filename);
    }

    out << std::setprecision(17);

    for (double d : vec)
    {
        out << d << ' ';
    }
}

//arg is coeff dir
int main(int argc, char* argv[])
{
    if (argc != 2)
    {
        std::cerr << "invalid arg count\n";
        return 1; 
    }

    const std::string dir(argv[1]);

    if (dir == "-h" || dir == "--help")
    {
        std::cerr << "The argument is for the folder that contains the coeffs\n";
        return 0;
    }

	const double temp = 1000.0;
	const double dz = 100e3;
	const double mass = 2.66e-26;
    
    const double t1 = 0.0;
	const double t2 = 1000.0;
	const double dt = 2.0;

    const double vmin = -6000.0;
    const double vmax = -vmin;
    const double dv = 5.0;

    Distrib2D dist{vmin, vmax, dv, mass, temp, dir};
    
    std::vector<double> times = dist.calc_times(t1, t2, dt);
    Distrib2D::Moments m = dist.get_moments(t1, t2, dt, dz);
    

    Distrib2D dist2{-3000.0, 3000.0, 10, mass, temp, dir};
    std::vector<array2d<double>> dists(10);
    std::vector<std::string> titles(10);
    
    for (int i = 0; i < 9; i++)
    {
        const double t = 100.0 + 50.0*i;
        titles[i] =  std::to_string((int)t);
        dists[i] = dist2.get_f_vf_dist(t, dz);
    }
   
    std::print("2d distribution dat file info. Rows:{}, Cols:{}\n", dists[0].dim0(), dists[0].dim1());
    std::print("writing outputs to disk\n");

    write_vec("./moment_output/time.dat", times);
    write_vec("./moment_output/density.dat", m.density);
    write_vec("./moment_output/avg_par_vel.dat", m.avg_par_velocity);
    write_vec("./moment_output/ion_temp_par.dat", m.ion_temp_par);
    write_vec("./moment_output/ion_temp_perp.dat", m.ion_temp_perp);

    //should be ok
    #pragma omp parallel for
    for (int i = 0; i < 9; i++)
    {
        write_matrix("./output/" + titles[i] + ".dat", dists[i]);
    }

	return 0;
}