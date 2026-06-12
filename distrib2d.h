#pragma once
#include <cmath>
#include <vector>
#include <numbers>


#include "ElectricField.h"

//TODO
//I want to have a trapz function for the integration
//instead of writing the integration every time

//currently a lot of the code is redundant for the parallel/perpendicular functions

class Distrib2D
{
public:
	Distrib2D(double v_min, double v_max, double v_step, double m, double T, std::string coeff_path) : 
		vf_min(v_min), vf_max(v_max), dv(v_step), mass(m), temp(T), ef(coeff_path)
	{
		nv = std::round((vf_max - vf_min) / v_step) + 1;

		//populate vf arrays
		vf_par_vec.resize(nv);
		for (int i = 0; i < nv; i++)
		{
			vf_par_vec[i] = vf_min + dv*i;
		}

		vf_perp_vec = vf_par_vec; //copy assignment. same, I think it is redundant, even in the future
	}
	
	//moments at time steps
	struct Moments
	{
		std::vector<double> density;
		std::vector<double> avg_par_velocity;
		std::vector<double> ion_temp_par;
		std::vector<double> ion_temp_perp;
	};

	std::vector<std::vector<double>> get_f_vf_dist(double t, double dz) const;
	Moments get_moments(double t1, double t2, double dt, double dz) const; // const or not?
	std::vector<double> calc_times(double t1, double t2, double dt) const;
	const std::vector<double>& get_vf_par() const { return vf_par_vec; }
	const std::vector<double>& get_vf_perp() const { return vf_perp_vec; }
	double electric_field(double t) const; // [mV/m]

private:
	double vf_max, vf_min, dv, mass, temp;
	std::vector<double> vf_par_vec, vf_perp_vec;
	int nv;
	ElectricField ef;

	static constexpr double pi = std::numbers::pi_v<double>;
	static constexpr double kb = 1.380649e-23; //https://physics.nist.gov/cgi-bin/cuu/Value?k
	static constexpr double acceleration = -5.0; //m s^-2
	static constexpr double e_field_ramp_rate = 1.0; //1 mV/m per second
	static constexpr double e_field_peak = 100.0; //  [mV/m]

	//double electric_field(double t) const; // [mV/m]

	double maxwell_boltzmann_dist(double v) const;
	double distribution_function(double v, double n, double t, int angle) const;
	double boundary_density(double t) const;

	double calc_dt(double vi, double vf, double t) const;
	double compute_vi_par(double vf_par, double t, double dz) const;
	double compute_f_vf_par(double vf_par, double t, double dz) const;
	double compute_f_vf_perp(double vf_par, double vf_perp, double t, double dz) const;
	//double compute_f_vf_2d(double vf_par, double vf_perp, double t, double dz) const;

	double calc_density(double t, double dz) const;

	double calc_avg_velocity_par(double t, double dz, double density) const;
	double calc_avg_velocity_perp(double t, double dz, double density) const;

	double calc_ion_temp_par(double t, double dz, double density, double avg_vel_par, double avg_vel_perp) const;
	double calc_ion_temp_perp(double t, double dz, double density, double avg_vel_par, double avg_vel_perp) const;
	
};
//2d shit do it ok


/*
	#pragma once
	#include <cmath>
	#include <vector>
	#include <numbers>


	#include "ElectricField.h"

	//TODO
	//I want to have a trapz function for the integration
	//instead of writing the integration every time6

	//currently a lot of the code is redundant for the parallel/perpendicular functions

	class Distrib2D
	{
	public:
		Distrib2D(double v_min, double v_max, double v_step, double m, double T, int asp_angle, std::string coeff_path) : 
			vf_min(v_min), vf_max(v_max), dv(v_step), mass(m), temp(T), angle(asp_angle), ef(coeff_path)
		{
			nv = std::lround((vf_max - vf_min) / v_step) + 1;
			e_field_ramp_rate = 1.0; //1 mV/m per second

			//populate vf arrays
			vf_par.resize(nv);
			for (int i = 0; i < nv; i++)
			{
				vf_par[i] = vf_min + dv*i;
			}

			vf_perp = vf_par; //copy assignment. same, I think it is redundant, even in the future
		}
		
		std::vector<std::vector<double>> get_f_vf_dist(double t, double dz) const;

		std::vector<double> get_f_vf_par_dist(double t, double dz) const;
		std::vector<double> get_f_vf_perp_dist(double t, double dz) const;

		std::vector<double> calc_zeroth_moment_over_time_par(double t1, double t2, double dt, double dz) const;
		std::vector<double> calc_zeroth_moment_over_time_perp(double t1, double t2, double dt, double dz) const;

		std::vector<double> calc_first_moment_over_time_par(double t1, double t2, double dt, double dz) const;
		std::vector<double> calc_first_moment_over_time_perp(double t1, double t2, double dt, double dz) const;

		std::vector<double> calc_ion_temp_over_time_par(double t1, double t2, double dt, double dz) const;
		std::vector<double> calc_ion_temp_over_time_perp(double t1, double t2, double dt, double dz) const;

		const std::vector<double>& get_vf_par() const { return vf_par; }//same as below
		const std::vector<double>& get_vf_perp() const { return vf_perp; }
		static std::vector<double> calc_times(double t1, double t2, double dt);

	private:
		double vf_max, vf_min, dv, mass, temp, e_field_ramp_rate;
		std::vector<double> vf_par, vf_perp;
		int nv, angle;
		ElectricField ef;

		static constexpr double pi = std::numbers::pi_v<double>;
		static constexpr double kb = 1.380649e-23; //https://physics.nist.gov/cgi-bin/cuu/Value?k
		static constexpr double acceleration = -5.0; //m s^-2

		double electric_field(double t) const; // mV/m

		double maxwell_boltzmann_dist(double v) const;
		double distribution_function(double v, double n, double t) const;

		double boundary_num_density(double t) const;

		double calc_dt(double vi, double vf, double t) const;

		double compute_vi_par(double vf_par, double t, double dz) const;

		double compute_f_vf(double vf_par, double vf_perp, double t, double dz) const;

		double calc_num_density(double t, double dz) const;

		double calc_avg_velocity_par(double t, double dz) const;
		double calc_avg_velocity_perp(double t, double dz) const;

		//TODO check with book
		double calc_ion_temp_par(double t, double dz) const;
		double calc_ion_temp_perp(double t, double dz) const;
		
	};

*/