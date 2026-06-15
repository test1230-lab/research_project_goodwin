#pragma once

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <stdexcept>
#include <ranges>
#include <string_view>
#include <filesystem>
#include <algorithm>
#include <utility>
#include <cmath>
#include <stdexcept>
#include <array>

#include <Eigen/Dense>

//Note: only works for the exact files provided
class ElectricField
{
public:
    ElectricField(std::string directory) : e_interp_coeffs{}
    {
        this->dir = directory;
        this->n_files = num_files();
        coeffs.resize(n_angles, std::vector<std::vector<double>>(n_files, std::vector<double>(n_cols)));
        ion_thermal_speeds.resize(n_angles, std::vector<double>(n_files));
        ion_thermal_speeds_interp_coeffs.resize(n_angles);
        read_coeffs();
        compute_interp_coeffs();
    }

    double compute_dist_discrete(int electric_field, int aspect_angle, double x) const;
    double compute_dist(double electric_field, int aspect_angle, double v) const;
    double compute_dist2(double electric_field, int aspect_angle, double v) const;
    const std::vector<double>& get_ion_thermal_speeds(int aspect_angle) const;
    double compute_integral(double electric_field, int aspect_angle) const; //check

private:
    static constexpr int n_cols = 26, n_angles = 10;
    static constexpr std::string_view delim{"\t\t"};
    static constexpr int order = 6;
    static constexpr int n_e_vals_ptable = 100'000;
    static constexpr double d_en = 1.0/n_e_vals_ptable;

    int n_files;
    std::filesystem::path dir;
    std::vector<double> e_field_vals;
    //"coeffs" is indexed as such: coeffs[aspect angle idx][e field idx][i]
    //"ion_thermal_speeds" is indexed as such: ion_thermal_speeds[aspect angle idx][e field idx]
    std::vector<std::vector<std::vector<double>>> coeffs;
    std::vector<std::vector<double>> ion_thermal_speeds;

    //indexed as such: e_interp_coeffs[angle idx][coeff idx][i]
    //where coeff idx is the index of the 
    std::array<std::array<std::array<double, order + 1>, n_cols>, n_angles> e_interp_coeffs;

    //indexed as such: ion_thermal_speeds_interp_coeffs[angle idx][i]
    std::vector<std::array<double, order + 1>> ion_thermal_speeds_interp_coeffs;

    int num_files() const;
    void read_coeffs();
    int extract_number(const std::filesystem::path& p) const;
    void compute_interp_coeffs();
    double eval_poly(const std::array<double, order + 1>& c, double x) const;
    double eval_legendre_series(const std::array<double, n_cols>& coeffs, double x) const;

};