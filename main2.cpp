#include <iostream>
#include <vector>
#include <string>
#include <stdexcept>
#include <iomanip>
#include <fstream>

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
            if (i != 0)
            {
                out << ' ';
            }

            out << mat[i, j];
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

int main()
{
	const double temp = 1000.0;
	const double dz = 100e3;
	const double mass = 2.66e-26;
    
    const double t1 = 0.0;
	const double t2 = 1000.0;
	const double dt = 2.0;

    const double vmin = -6000.0;
    const double vmax = -vmin;
    const double dv = 10;

    Distrib2D dist{vmin, vmax, dv, mass, temp, "./coeffs"};
    
    std::vector<double> times = dist.calc_times(t1, t2, dt);
    Distrib2D::Moments m = dist.get_moments(t1, t2, dt, dz);
    

    std::vector<array2d<double>> dists(10);
    std::vector<std::string> titles(10);
    
    for (int i = 0; i < 9; i++)
    {
        const double t = 100.0 + 50.0*i;
        titles[i] =  std::to_string((int)t);
        dists[i] = dist.get_f_vf_dist(t, dz);
    }
   


    write_vec("./moment_output/time.dat", times);
    write_vec("./moment_output/density.dat", m.density);
    write_vec("./moment_output/avg_par_vel.dat", m.avg_par_velocity);
    write_vec("./moment_output/ion_temp_par.dat", m.ion_temp_par);
    write_vec("./test.dat", m.ion_temp_par);
    write_vec("./moment_output/ion_temp_perp.dat", m.ion_temp_perp);

    for (int i = 0; i < 9; i++)
    {
        write_matrix("./output/" + titles[i] + ".dat", dists[i]);
    }

	return 0;
}