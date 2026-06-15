#include <iostream>
#include <vector>
#include <string>
#include <stdexcept>
#include <utility>
#include <tuple>
#include <iomanip>
#include <fstream>
#include <chrono>

#include "distrib2d.h"
#include "gnuplot-iostream.h"

// -----> CLAUDE CODE <-----
// ----------------------------------------------------------------------------
// Minimal 1-D plotter using gnuplot-iostream.
// Empty `savefile` -> interactive window; otherwise a PNG is written.
// ----------------------------------------------------------------------------



void plot1d(const std::vector<double>& x, const std::vector<double>& y,
            const std::string& title  = "", const std::string& xlabel = "x",
            const std::string& ylabel = "y", const std::string& savefile = "")
{
    if (x.size() != y.size())
    {
        throw std::invalid_argument("plot1d: x and y must have equal length");
    }

    Gnuplot gp;

    if (!savefile.empty()) 
    {
        gp << "set terminal pngcairo size 900,600 enhanced\n";
        gp << "set output '" << savefile << "'\n";
    }

    if (!title.empty())
    {
        gp << "set title '" << title << "'\n";
    }

    gp << "set xlabel '" << xlabel << "'\n";
    gp << "set ylabel '" << ylabel << "'\n";
    gp << "set grid\n";
    gp << "plot '-' with lines lt -1 lw 1 notitle\n";
    gp.send1d(std::make_pair(x, y));
}

//NOTE: -----> CLAUDE CODE!!!! <-----
//--------------------------------------------
void heatmap(const std::vector<double>& vpar,                  // x, size Ncol
             const std::vector<double>& vperp,                 // y, size Nrow
             const std::vector<std::vector<double>>& z,        // z[row][col]
             const std::string& title    = "",
             const std::string& xlabel   = "v_{par} [m/s]",
             const std::string& ylabel   = "v_{perp} [m/s]",
             const std::string& clabel   = "",
             const std::string& savefile = "")
{
    if (z.size() != vperp.size())
        throw std::invalid_argument("heatmap: z row count must equal vperp.size()");
    for (const auto& row : z)
        if (row.size() != vpar.size())
            throw std::invalid_argument("heatmap: each z row must have vpar.size() columns");

    Gnuplot gp;

    if (!savefile.empty()) {
        gp << "set terminal pngcairo size 900,600 enhanced\n";
        gp << "set output '" << savefile << "'\n";
    }
    if (!title.empty())
        gp << "set title '" << title << "'\n";

    gp << "set xlabel '" << xlabel << "'\n";
    gp << "set ylabel \"" << ylabel << "\"\n";   // double-quoted: see note
    if (!clabel.empty())
        gp << "set cblabel '" << clabel << "'\n";

    std::vector<std::vector<std::tuple<double, double, double>>> grid;
    grid.reserve(vperp.size());
    for (std::size_t r = 0; r < vperp.size(); ++r) {
        grid.emplace_back();
        grid.back().reserve(vpar.size());
        for (std::size_t c = 0; c < vpar.size(); ++c)
            grid.back().emplace_back(vpar[c], vperp[r], z[r][c]);
    }

    // file2d writes grid to a temp file and returns its (quoted) name.
    // Do NOT also call send2d — file2d both writes and references the data.
    gp << "set autoscale fix\n";
    gp << "set size square\n";
    gp << "plot " << gp.file2d(grid) << " with image notitle\n";
}

//NOTE: -----> CLAUDE CODE!!!! <-----
void heatmap_grid(const std::vector<double>& vpar,
                  const std::vector<double>& vperp,
                  const std::vector<std::vector<std::vector<double>>>& panels,
                  const std::vector<std::string>& titles   = {},
                  const std::string& suptitle = "",
                  const std::string& xlabel   = "v_{par} [m/s]",
                  const std::string& ylabel   = "v_{perp} [m/s]",
                  const std::string& clabel   = "",
                  const std::string& savefile = "")
{
    const std::size_t n = panels.size();
    if (n == 0) throw std::invalid_argument("heatmap_grid: no panels");
    for (const auto& z : panels) {
        if (z.size() != vperp.size())
            throw std::invalid_argument("heatmap_grid: panel rows must equal vperp.size()");
        for (const auto& row : z)
            if (row.size() != vpar.size())
                throw std::invalid_argument("heatmap_grid: panel cols must equal vpar.size()");
    }

    const std::size_t cols = static_cast<std::size_t>(std::ceil(std::sqrt(double(n))));
    const std::size_t rows = (n + cols - 1) / cols;

    // ONE shared color range across all panels, so brightness == magnitude everywhere
    double lo = +1e30, hi = 0.0;
    for (const auto& z : panels)
        for (const auto& row : z)
            for (double v : row) { lo = std::min(lo, v); hi = std::max(hi, v); }

    Gnuplot gp;

    const std::size_t W = 480 * cols + 160, H = 460 * rows + 80;
    if (!savefile.empty()) {
        gp << "set terminal pngcairo size " << W << "," << H << " enhanced font ',13'\n";
        gp << "set output '" << savefile << "'\n";
    } else {
        // swap 'qt' for 'wxt' or 'x11' if that's what your gnuplot build uses
        gp << "set terminal qt size " << W << "," << H << " font ',13'\n";
    }

    gp << "set autoscale fix\n";
    gp << "set cbrange [" << lo << ":" << hi << "]\n";
    gp << "set xtics 3000 font ',11'\n";
    gp << "set ytics 3000 font ',11'\n";
    gp << "set tics scale 0.5\n";
    gp << "set cbtics font ',12'\n";                                  // <-- add
    if (!clabel.empty()) gp << "set cblabel '" << clabel << "'\n";

    // position the single shared colorbar (screen coords), then disable it
    //gp << "set colorbox user origin 0.90,0.5 size 0.025,0.78\n";     // <-- changed
    gp << "unset colorbox\n";

    gp << "set multiplot layout " << rows << "," << cols
       << " margins 0.08,0.88,0.09,0.90 spacing 0.05,0.05";
    if (!suptitle.empty()) gp << " title '" << suptitle << "' font ',23'";
    gp << "\n";

    for (std::size_t p = 0; p < n; ++p) {
        const std::size_t row = p / cols, col = p % cols;

        if (p < titles.size() && !titles[p].empty())
            gp << "set title '" << titles[p] << "' font ',14' offset 0,-0.8\n";
        else
            gp << "unset title\n";

        // labels only on the outer edge
        if (col == 0)        gp << "set ylabel \"" << ylabel << "\"\n";
        else                 gp << "unset ylabel\n";
        if (row == rows - 1) gp << "set xlabel '" << xlabel << "'\n";
        else                 gp << "unset xlabel\n";

        if (p == n - 1) 
        {
            gp << "set colorbox\n";   // draw the shared bar once
            gp << "set colorbox user origin 0.90,0.12 size 0.025,0.78\n";   // full geometry here
        }

        const auto& z = panels[p];
        std::vector<std::vector<std::tuple<double, double, double>>> grid;
        grid.reserve(vperp.size());
        for (std::size_t r = 0; r < vperp.size(); ++r) {
            grid.emplace_back();
            grid.back().reserve(vpar.size());
            for (std::size_t c = 0; c < vpar.size(); ++c)
                grid.back().emplace_back(vpar[c], vperp[r], z[r][c]);
        }
        gp << "plot " << gp.file2d(grid) << " with image notitle\n";
    }

    gp << "unset multiplot\n";
}

//chatgpt code:
// Overlayed 1D plot: same x, multiple y series
void plot1d_overlay(const std::vector<double>& x,
                    const std::vector<std::vector<double>>& ys,
                    const std::vector<std::string>& labels = {},
                    std::string title = "",
                    std::string xlabel = "x",
                    std::string ylabel = "y",
                    std::string savefile = "")
{
    if (ys.empty())
        throw std::invalid_argument("plot1d_overlay: no y-series provided");

    for (const auto& y : ys)
    {
        if (y.size() != x.size())
            throw std::invalid_argument("plot1d_overlay: x and y must have equal length");
    }

    if (!labels.empty() && labels.size() != ys.size())
        throw std::invalid_argument("plot1d_overlay: labels must match number of y-series");

    FILE* gp = popen("gnuplot -persist", "w");
    if (!gp)
        throw std::runtime_error("plot1d_overlay: failed to launch gnuplot");

    if (!savefile.empty())
    {
        std::fprintf(gp, "set terminal pngcairo size 900,600 enhanced\n");
        std::fprintf(gp, "set output '%s'\n", savefile.c_str());
    }

    if (!title.empty())
        std::fprintf(gp, "set title '%s'\n", title.c_str());

    std::fprintf(gp, "set xlabel '%s'\n", xlabel.c_str());
    std::fprintf(gp, "set ylabel '%s'\n", ylabel.c_str());
    std::fprintf(gp, "set grid\n");

    std::fprintf(gp, "plot ");
    for (size_t k = 0; k < ys.size(); ++k)
    {
        if (k > 0) std::fprintf(gp, ", ");
        const char* label = labels.empty() ? "" : labels[k].c_str();
        std::fprintf(gp, "'-' with lines lw 1 title '%s'", label);
    }
    std::fprintf(gp, "\n");

    for (const auto& y : ys)
    {
        for (size_t i = 0; i < x.size(); ++i)
            std::fprintf(gp, "%g %g\n", x[i], y[i]);
        std::fprintf(gp, "e\n");
    }

    if (!savefile.empty())
        std::fprintf(gp, "unset output\n");

    std::fflush(gp);
    pclose(gp);
}

void write_matrix(const std::string& filename, const std::vector<std::vector<double>>& mat)
{
    std::ofstream out(filename);

    if (!out)
    {
        throw std::runtime_error("Failed to open " + filename);
    }

    out << std::setprecision(17);

    for (const auto& row : mat)
    {
        for (int i = 0; i < row.size(); i++)
        {
            if (i != 0)
            {
                out << ' ';
            }

            out << row[i];
        }
        out << '\n';
    }
}


int main()
{
	const double temp = 1000.0;
	const double dz = 100e3;
	const double mass = 2.66e-26;
    
    const double t1 = 0.0;
	const double t2 = 1000.0;
	const double dt = 5;

    const double vmin = -3000.0;
    const double vmax = -vmin;
    const double dv = 5.0;
    

    Distrib2D dist{vmin, vmax, dv, mass, temp, "./coeffs"};

    /*std::vector<double> times = dist.calc_times(t1, t2, dt);
    std::vector<double> ef(times.size());

    for (int i = 0; i < times.size(); i++)
    {
        const double t = t1 + dt*i;
        ef[i] = dist.electric_field(t);
    }

    plot1d(times, ef, "Electric Field Magnitude", "time [s]", "|E| [mV/m]");*/

    std::vector<double> times = dist.calc_times(t1, t2, dt);
    auto start = std::chrono::steady_clock::now();

    
    //Distrib2D::Moments m = dist.get_moments(t1, t2, dt, dz);
    
    std::vector<std::vector<std::vector<double>>> dists(10);
    std::vector<std::string> titles(10);
    #pragma omp parallel for
    for (int i = 0; i < 9; i++)
    {
        const double t = 100.0 + 50.0*i;
        titles[i] =  std::to_string((int)t);
        dists[i] = dist.get_f_vf_dist(t, dz);
    }

    /*std::vector<double> bd(times.size());
    for (int i = 0; i < times.size(); i++)
    {
        const double t = t1 + dt*i;
        bd[i] = dist.boundary_density(t);
    }*/

    auto end = std::chrono::steady_clock::now();
    uint64_t elapsed1 = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << elapsed1/1000.0 << "[s]\n";
    
    /*const std::string param = std::format("for dz={:.1f}km", dz/1000.0);
    plot1d(times, m.density, "Density " + param, "t [s]", "Density [m^-3]");
    plot1d(times, m.avg_par_velocity, "Parallel Veloctiy " + param, "t [s]", "Velocity [m/s]");
    plot1d_overlay(times, {m.ion_temp_par, m.ion_temp_perp},
             {"Ion Temp Parallel " + param, "Ion Temp Perpendicular " + param}, " ", "t [s]", "Temp [K]");

    plot1d(times, bd, "Boundary Density", "time [s]", "m^-3");*/


    for (int i = 0; i < 9; i++)
    {
        write_matrix("./output/" + titles[i] + ".dat", dists[i]);
    }

	return 0;
}

