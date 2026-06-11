#include "distrib2d.h"
#include <iostream>
//sloppy and slow at the moment

// mV/m
double Distrib2D::electric_field(double t) const
{
	if (t <= 0)
	{
		return 0.0;
	}

    const double ramp = t*e_field_ramp_rate;
    if (ramp < e_field_peak)// rising edge
    {
    	return ramp;
    }

    return std::max(0.0, e_field_peak - (ramp - e_field_peak));   // falling edge, min of 0
}

//1d or 3d??? temp????
double Distrib2D::maxwell_boltzmann_dist(double v) const
{
	const double T = temp;
    const double a = std::sqrt(mass / (2.0*pi*kb*T));
    return a*std::exp(-mass*v*v / (2.0*kb*T));
}

//TODO: check!!!
//what angle!?
//maybe precompute and interpolate
double Distrib2D::distribution_function(double v, double n, double t, int angle) const
{
	const int a = 45;
	const double efield = electric_field(t);
	if (efield < 20.0)
	{
		const double u = std::clamp(efield/20.0, 0.0, 1.0);
		const double integral = ef.compute_integral(20.0, a);
		const double f_e20 = ef.compute_dist(20.0, a, v)/integral;
		const double f_maxwellian = maxwell_boltzmann_dist(v);
		return n*(f_maxwellian + u * (f_e20 - f_maxwellian));
	}

	const double integral = ef.compute_integral(efield, a);
	return n*ef.compute_dist(efield, angle, v)/integral;
}



double Distrib2D::boundary_density(double t) const
{
	return 1.0;
}

double Distrib2D::calc_dt(double vi, double vf, double t) const
{
	return (vf - vi)/acceleration;
}

//should throw if sqrt arg is less than zero
double Distrib2D::compute_vi_par(double vf_par, double t, double dz) const
{
	return std::sqrt(vf_par*vf_par - 2.0*acceleration*dz);
}


double Distrib2D::compute_f_vf_par(double vf_par, double t, double dz) const
{
	const double vi = compute_vi_par(vf_par, t, dz);
	const double dt = calc_dt(vi, vf_par, t);
	return distribution_function(vi, boundary_density(t - dt), t - dt, 0);
}

double Distrib2D::compute_f_vf_perp(double vf_par, double vf_perp, double t, double dz) const
{
	const double vi = vf_perp;
	const double dt = calc_dt(compute_vi_par(vf_par, t, dz), vf_par, t);
	return distribution_function(vi, boundary_density(t - dt), t - dt, 90);
}

/*
double Distrib2D::compute_f_vf_2d(double vf_par, double vf_perp, double t, double dz) const
{
	const double vi_perp = vf_perp; //no forces acting on it
	const double vi_par = compute_vi_par(vf_par, t, dz);
	const double dt = calc_dt(vi_par, vf_par, t); //same particle

	const double density = boundary_density(t - dt);

	const double f_vf_par = distribution_function(vi_par, density, t - dt, 0);
	const double f_vf_perp = distribution_function(vi_perp, density, t - dt, 90);

	//TODO: is this correct?
	return f_vf_par*f_vf_perp/density;
}*/

std::vector<std::vector<double>> Distrib2D::get_f_vf_dist(double t, double dz) const
{
	std::vector<std::vector<double>> result(nv, std::vector<double>(nv));

	#pragma omp parallel for
	for (int i = 0; i < nv; i++)
	{
		const double vf_par = vf_par_vec[i];
		const double vi_par = compute_vi_par(vf_par, t, dz);
		const double dt = calc_dt(vi_par, vf_par, t); //same particle
		const double bnd_density = boundary_density(t - dt);
		const double f_vf_par = distribution_function(vi_par, bnd_density, t - dt, 0);

		for (int j = 0; j < nv; j++)
		{
			const double vf_perp = vf_perp_vec[j];
			const double f_vf_perp = distribution_function(vf_perp, bnd_density, t - dt, 90);
			result[j][i] = f_vf_par*f_vf_perp/bnd_density;
		}
	}

	return result;
}



double Distrib2D::calc_density(double t, double dz) const
{
	double integral = 0.0;

    for (int i = 0; i < nv; i++)
    {
		const double vf_par = vf_par_vec[i];
		const double vi_par = compute_vi_par(vf_par, t, dz);
		const double dt = calc_dt(vi_par, vf_par, t); //same particle
		const double bnd_density = boundary_density(t - dt);
		const double f_vf_par = distribution_function(vi_par, bnd_density, t - dt, 0);

        const double wi = 1.0 - 0.5 * ((i == 0) + (i == nv - 1));
        for (int j = 0; j < nv; j++)
        {
			const double wj = 1.0 - 0.5 * ((j == 0) + (j == nv - 1));
			const double vf_perp = vf_perp_vec[j];
			const double f_vf_perp = distribution_function(vf_perp, bnd_density, t - dt, 90);
			const double dist2d = f_vf_par*f_vf_perp/bnd_density;
			integral += wi*wj*dist2d;
        }
    }

    return dv*dv*integral;
}

double Distrib2D::calc_avg_velocity_par(double t, double dz, double density) const
{
	double integral = 0.0;
	for (int i = 1; i < nv - 1; i++)
	{
		const double vi = compute_f_vf_par(vf_par_vec[i], t, dz);
		integral += vf_par_vec[i]*vi;
	}

	const double f_0 = vf_par_vec[0]*compute_f_vf_par(vf_par_vec[0], t, dz);
	const double f_n1 = vf_par_vec[nv - 1]*compute_f_vf_par(vf_par_vec[nv - 1], t, dz);

	integral = dv * (integral + (f_0 + f_n1)/2.0);

	return integral / density; //should it be calling the 2d density function??
}

double Distrib2D::calc_avg_velocity_perp(double t, double dz, double density) const
{
	double integral = 0.0;
	for (int i = 1; i < nv - 1; i++)
	{
		const double vi = compute_f_vf_perp(vf_par_vec[i],vf_perp_vec[i], t, dz);
		integral += vf_perp_vec[i]*vi;
	}

	const double f_0 = vf_perp_vec[0]*compute_f_vf_perp(vf_par_vec[0], vf_perp_vec[0], t, dz);
	const double f_n1 = vf_perp_vec[nv - 1]*compute_f_vf_perp(vf_par_vec[nv - 1], vf_perp_vec[nv - 1], t, dz);

	integral = dv * (integral + (f_0 + f_n1)/2.0);

	return integral / density; //should it be calling the 2d density function??
}

//utility function
std::vector<double> Distrib2D::calc_times(double t1, double t2, double dt) const
{
	const int N = std::round((t2 - t1) / dt);
	std::vector<double> times(N);

	for (int i = 0; i < N; i++)
	{
		times[i] = t1 + i*dt;
	}

	return times;
}

double Distrib2D::calc_ion_temp_par(double t, double dz, double density, double avg_velocity) const
{
	const double u = avg_velocity;

	double integral = 0.0;
	for (int i = 1; i < nv - 1; i++)
	{
		integral += vf_par_vec[i]*vf_par_vec[i]*compute_f_vf_par(vf_par_vec[i], t, dz);
	}

	const double f_0 = vf_par_vec[0]*vf_par_vec[0]*compute_f_vf_par(vf_par_vec[0], t, dz);
	const double f_n1 = vf_par_vec[nv - 1]*vf_par_vec[nv - 1]*compute_f_vf_par(vf_par_vec[nv - 1], t, dz);
	integral = dv * (integral + (f_0 + f_n1)/2.0);

	return mass * (integral / density - u*u) / kb;
}

//TODO: complete
double Distrib2D::calc_ion_temp_perp(double t, double dz, double density, double avg_velocity) const
{
	const double u = avg_velocity;

	double integral = 0.0;
	for (int i = 1; i < nv - 1; i++)
	{
		integral += vf_perp_vec[i]*vf_perp_vec[i]*compute_f_vf_perp(vf_par_vec[i], vf_perp_vec[i], t, dz);
	}

	const double f_0 = vf_perp_vec[0]*vf_perp_vec[0]*compute_f_vf_perp(vf_par_vec[0], vf_perp_vec[0], t, dz);
	const double f_n1 = vf_perp_vec[nv - 1]*vf_perp_vec[nv - 1]*compute_f_vf_perp(vf_par_vec[nv - 1], vf_perp_vec[nv - 1], t, dz);
	integral = dv * (integral + (f_0 + f_n1)/2.0);

	return mass * (integral / density - u*u) / kb;
}

//TODO: complete
Distrib2D::Moments Distrib2D::get_moments(double t1, double t2, double dt, double dz) const
{
	const int N = std::round((t2 - t1) / dt);

	Moments m;
	m.density.resize(N);
	m.avg_par_velocity.resize(N);
	m.ion_temp_par.resize(N);
	m.ion_temp_perp.resize(N);

	#pragma omp parallel for
	for (int i = 0; i < N; i++)
	{
		const double t = t1 + dt*i;
		const double density = calc_density(t, dz);
		const double avg_velocity = calc_avg_velocity_par(t, dz, density);
		m.density[i] = density;
		m.avg_par_velocity[i] = avg_velocity;
		m.ion_temp_par[i] = calc_ion_temp_par(t, dz, density, avg_velocity);
		m.ion_temp_perp[i] = calc_ion_temp_perp(t, dz, density, calc_avg_velocity_perp(t, dz, density));
	}

	return m;
}

