#pragma once

#include <fstream>
#include <string>
#include <vector>
#include <stdexcept>
#include <ranges>
#include <string_view>
#include <filesystem>
#include <algorithm>
#include <utility>
#include <cmath>
#include <array>
#include <print>

#include "mdarray.h"

#include <Eigen/Dense>

//TODO: use std::mdspan
//Note: only works for the exact files provided
class ElectricField
{
public:
    ElectricField(std::string directory) : e_interp_coeffs{}
    {
        this->dir = directory;
        this->n_files = num_files();
        coeffs.resize(n_angles, n_files, n_cols);
        ion_thermal_speeds.resize(n_angles, std::vector<double>(n_files));
        //read_coeffs();
        read_coeffs_new_fmt();
        compute_interp_coeffs();
    }

    double compute_dist_discrete(int electric_field, int aspect_angle, double x) const;
    double compute_dist(double electric_field, int aspect_angle, double v) const;
    const std::vector<double>& get_ion_thermal_speeds(int aspect_angle) const;
    double compute_integral(double electric_field, int aspect_angle) const; //check
    void print_coeffs() const;

private:
    static constexpr int n_cols = 26, n_angles = 10;
    static constexpr std::string_view delim{"\t\t"};
    static constexpr int order = 6;
    static constexpr int n_e_vals_ptable = 100'000;
    static constexpr double d_en = 1.0/n_e_vals_ptable;

    static constexpr std::array<double, 2*n_cols> leg_a = [] {
        std::array<double, 2*n_cols> a{};
        for (int i = 1; i < 2*(n_cols - 1); i++) a[i] = (2.0*i + 1.0)/(i + 1.0);
        return a;
    }();

    static constexpr std::array<double, 2*n_cols> leg_b = [] {
        std::array<double, 2*n_cols> b{};
        for (int i = 1; i < 2*(n_cols - 1); i++) b[i] = double(i)/(i + 1.0);
        return b;
    }();

    int n_files;
    std::filesystem::path dir;
    std::vector<double> e_field_vals;
    //"coeffs" is indexed as such: coeffs[aspect angle idx][e field idx][i]
    //"ion_thermal_speeds" is indexed as such: ion_thermal_speeds[aspect angle idx][e field idx]
    array3d<double> coeffs;
    std::vector<std::vector<double>> ion_thermal_speeds;

    //indexed as such: e_interp_coeffs[angle idx][coeff idx][i]
    //where coeff idx is the index of the 
    std::array<std::array<std::array<double, order + 1>, n_cols>, n_angles> e_interp_coeffs;

    //indexed as such: ion_thermal_speeds_interp_coeffs[angle idx][i]
    std::array<std::array<double, order + 1>, n_angles> ion_thermal_speeds_interp_coeffs;

    int num_files() const;
    void read_coeffs();
    void read_coeffs_new_fmt();
    int extract_number(const std::filesystem::path& p) const;
    void compute_interp_coeffs();
    double eval_poly(const std::array<double, order + 1>& c, double x) const;
    double eval_legendre_series(const std::array<double, n_cols>& coeffs, double x) const;        
};