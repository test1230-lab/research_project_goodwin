#include "ElectricField.h"

int ElectricField::num_files() const
{
    int n = 0;
    for (const auto& entry : std::filesystem::directory_iterator(dir))
    {
        if (entry.is_regular_file()) 
        {
            n++;
        }
    }

    std::print("{} coefficient files in directory\n", n);

    return n;
}

//chatgpt code
int ElectricField::extract_number(const std::filesystem::path& p) const
{
    const std::string name = p.stem().string();  // "1-E_100"

    const auto pos = name.find('_');
    if (pos == std::string::npos)
    {
        throw std::runtime_error("Invalid filename: " + name);
    }

    return std::stoi(name.substr(pos + 1));
}

/*
void ElectricField::read_coeffs()
{
    std::vector<std::pair<int, std::filesystem::path>> files;

    for (const auto& entry : std::filesystem::directory_iterator(dir))
    {
        if (entry.is_regular_file())
        {
            const int num = extract_number(entry.path());
            e_field_vals.push_back(static_cast<double>(num));
            files.push_back({num, entry.path()});
        }
    }

    std::sort(files.begin(), files.end());
    std::sort(e_field_vals.begin(), e_field_vals.end());

    int n = 0;
    for (const auto& file : files)
    {
        std::filesystem::directory_entry entry(file.second);
        if (entry.is_regular_file()) 
        {
            std::ifstream in(entry.path());
            std::string line;

            //skip line 1,2
            std::getline(in, line);
            std::getline(in, line);

            //get ion thermal speeds
            std::getline(in, line);
            int idx = 0;
            for (auto&& elem : std::views::split(line, delim))
            {
                const double val = std::stod(std::string(elem.begin(), elem.end()));
                ion_thermal_speeds[idx++][n] = val;
            }

            //get coeffs
            for (int i = 0; i < n_cols; i++)
            {
                std::getline(in, line);
                idx = 0;
                for (auto&& elem : std::views::split(line, delim))
                {
                    const double val = std::stod(std::string(elem.begin(), elem.end()));
                    coeffs[idx++, n, i] = val;
                }
            }

            n++;
        }
    }
}
*/


void ElectricField::read_coeffs()
{
    std::vector<std::pair<int, std::filesystem::path>> files;

    for (const auto& entry : std::filesystem::directory_iterator(dir))
    {
        if (entry.is_regular_file())
        {
            const int num = extract_number(entry.path());
            e_field_vals.push_back(static_cast<double>(num));
            files.push_back({num, entry.path()});
        }
    }

    std::sort(files.begin(), files.end());
    std::sort(e_field_vals.begin(), e_field_vals.end());

    int n = 0;
    for (const auto& file : files)
    {
        std::filesystem::directory_entry entry(file.second);
        if (entry.is_regular_file()) 
        {
            std::ifstream in(entry.path());
            std::string line;

            //skip line 1,2
            std::getline(in, line);
            std::getline(in, line);

            //get ion thermal speeds
            std::getline(in, line);
            int idx = 0;
            for (auto&& elem : std::views::split(line, delim))
            {
                const double val = std::stod(std::string(elem.begin(), elem.end()));
                ion_thermal_speeds[idx++][n] = val;
            }

            //get coeffs
            for (int i = 0; i < n_cols; i++)
            {
                std::getline(in, line);
                idx = 0;
                for (auto&& elem : std::views::split(line, delim))
                {
                    const double val = std::stod(std::string(elem.begin(), elem.end()));
                    coeffs[idx++, n, i] = val;
                }
            }

            n++;
        }
    }
}

void ElectricField::read_coeffs_new_fmt()
{
    std::vector<std::pair<int, std::filesystem::path>> files;

    for (const auto& entry : std::filesystem::directory_iterator(dir))
    {
        if (entry.is_regular_file())
        {
            const int num = extract_number(entry.path());
            e_field_vals.push_back(static_cast<double>(num));
            files.push_back({num, entry.path()});
        }
    }

    std::sort(files.begin(), files.end());
    std::sort(e_field_vals.begin(), e_field_vals.end());

    int n = 0;
    for (const auto& file : files)
    {
        std::filesystem::directory_entry entry(file.second);
        if (entry.is_regular_file()) 
        {
            std::ifstream in(entry.path());
            std::print("path:{}\n", entry.path().string());
            std::string line;

            //skip line 1-8 inclusive
            for (int i = 0; i < 8; i++)
            {
                std::getline(in, line);
            }

            //get ion thermal speeds
            std::getline(in, line);
            int idx = 0;
            for (auto&& elem : std::views::split(line, delim))
            {
                const double val = std::stod(std::string(elem.begin(), elem.end()));
                ion_thermal_speeds[idx++][n] = val;
            }

            //get coeffs
            for (int i = 0; i < 2*n_cols - 1; i++)
            {
                std::getline(in, line);
                if (i % 2 == 0)
                {
                    const int coeff_idx = i / 2;
                    idx = 0;
                    for (auto&& elem : std::views::split(line, delim))
                    {
                        const double val = std::stod(std::string(elem.begin(), elem.end()));
                        coeffs[idx++, n, coeff_idx] = val;
                        
                    }
                }
            }
            

            n++;
        }
    }
}

void ElectricField::print_coeffs() const
{
    for (int i = 0; i < n_files; i++)
    {
        std::print("E = {}\n", e_field_vals[i]);
        for (int j = 0; j < n_cols; j++)
        {
            std::print("{}\n", coeffs[0, j, i]);
        }
    }

    std::print("\n\n");
}

double ElectricField::compute_dist_discrete(int electric_field, int aspect_angle, double v) const
{
    const int aspect_angle_idx = aspect_angle / 10;
    const int e_idx = std::max(0, (electric_field) / 10);

    //y is vx/b, where vx is the line-of-sight speed and b is the ion thermal speed).
    const double y = v/ion_thermal_speeds[aspect_angle_idx][e_idx];

    double sum = 0.0;
    for (int i = 0; i < n_cols; i++)
    {
        const double ci = coeffs[aspect_angle_idx, e_idx, i];
        sum += ci*std::legendre(i*2, y/4.0);
    }

    return sum/ion_thermal_speeds[aspect_angle_idx][e_idx];
}

double ElectricField::eval_poly(const std::array<double, order + 1>& c, double x) const
{
    //horners method
    double result = c[order];

    for (int k = order - 1; k >= 0; k--)
    {
        //result = result*x + c[k];
        result = std::fma(result, x, c[k]);
    }

    return result;
}

double ElectricField::compute_dist(double electric_field, int aspect_angle, double v) const
{
    const int angle_idx = aspect_angle / 10;
    electric_field = std::clamp(electric_field, 0.0, 200.0);
    //const double e_norm = (electric_field - 20.0)/90.0 - 1.0;
    const double e_norm = 2.0*((electric_field - e_field_vals[0])/(e_field_vals[n_files - 1] - e_field_vals[0])) - 1.0;

    const double ion_thermal_speed_interp = eval_poly(ion_thermal_speeds_interp_coeffs[angle_idx], e_norm);
    const double y = v/ion_thermal_speed_interp;
    const double x = y/4.0;

    if (std::abs(x) > 1.0)
    {
        return 0.0;
    }

    std::array<double, n_cols> coeffs_;
    for (int i = 0; i < n_cols; i++)
    {
        coeffs_[i] = eval_poly(e_interp_coeffs[angle_idx][i], e_norm);
    }

    return eval_legendre_series(coeffs_, x)/ion_thermal_speed_interp;
}

double ElectricField::eval_legendre_series(const std::array<double, n_cols>& c, double x) const 
{
    double sum = c[0];
    double pm1 = 1.0;   // P_{i-1}
    double p   = x;     // P_i  (i = 1)

    for (int i = 1; i < 2*(n_cols - 1); i++) 
    {
        const double pnext = leg_a[i]*x*p - leg_b[i]*pm1;  //no division
        pm1 = p; p = pnext;
        if (i % 2 == 1)
        {
            sum = std::fma(p, c[(i + 1)/2], sum);
        }
    }

    return sum;
}

const std::vector<double>& ElectricField::get_ion_thermal_speeds(int aspect_angle) const
{
    const int aspect_angle_idx = aspect_angle / 10;
    return ion_thermal_speeds[aspect_angle_idx];
}

void ElectricField::compute_interp_coeffs()
{
    std::vector<double> e_coeffs_at_angle(n_files);
    std::vector<double> norm_e_field_vals(n_files);

    //[0, 200] -> [-1, 1]
    for (int i = 0; i < n_files; i++)
    {
        norm_e_field_vals[i] = 2.0*((e_field_vals[i] - e_field_vals[0])/(e_field_vals[n_files - 1] - e_field_vals[0])) - 1.0;
    }

    //using code from: https://gist.github.com/tabbott/129a32e965fefd202fb9065a4368986f
    auto xvec = Eigen::Map<const Eigen::VectorXd>(norm_e_field_vals.data(), norm_e_field_vals.size());
    Eigen::MatrixXd xs(norm_e_field_vals.size(), order + 1);
    xs.col(0).setOnes();

    for (int i = 1; i <= order; i++) 
    {
        xs.col(i).array() = xs.col(i - 1).array() * xvec.array();
    }

    auto decomposition = xs.householderQr();

    for (int angle_idx = 0; angle_idx < n_angles; angle_idx++)
    {
        for (int coeff_idx = 0; coeff_idx < n_cols; coeff_idx++)
        {
            //remap coeff array to more useful format
            for (int i = 0; i < n_files; i++)
            {
                e_coeffs_at_angle[i] = coeffs[angle_idx, i, coeff_idx];
            }     

            auto ys = Eigen::Map<const Eigen::VectorXd>(e_coeffs_at_angle.data(),
                                                        e_coeffs_at_angle.size());

            auto result_map = Eigen::Map<Eigen::VectorXd>(e_interp_coeffs[angle_idx][coeff_idx].data(),
                                                          e_interp_coeffs[angle_idx][coeff_idx].size());

            result_map = decomposition.solve(ys);
        }
    }
    
    for (int angle_idx = 0; angle_idx < n_angles; angle_idx++)
    {
        auto ys = Eigen::Map<const Eigen::VectorXd>(ion_thermal_speeds[angle_idx].data(),
                                                        ion_thermal_speeds[angle_idx].size());

        auto result_map = Eigen::Map<Eigen::VectorXd>(ion_thermal_speeds_interp_coeffs[angle_idx].data(),
                                                          ion_thermal_speeds_interp_coeffs[angle_idx].size());

        result_map = decomposition.solve(ys);
    }
    
}

double ElectricField::compute_integral(double electric_field, int aspect_angle) const
{
    const int angle_idx = aspect_angle / 10;
    const double e_norm = 2.0*((electric_field - e_field_vals[0])/(e_field_vals[n_files - 1] - e_field_vals[0])) - 1.0;
    return 8.0*eval_poly(e_interp_coeffs[angle_idx][0], e_norm);
}

