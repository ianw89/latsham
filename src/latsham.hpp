#pragma once

#include <errno.h>
#include <vector>
#include <math.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>
#include <cstdio>
#include <string>
#include <fstream>
#include <iostream>
#include <boost/math/tools/roots.hpp>
#include <limits>
#include <stdlib.h>

#define SILENT 0
#define VERBOSE 0

using boost::math::tools::eps_tolerance;
using boost::math::tools::toms748_solve;

const std::string GAL_MAGR_DENSITY = "/mount/sirocco1/imw2293/GROUP_CAT/SelfCalGroupFinder/py/parameters/bgs_y3/gal_abs_mag_r_k_best_density_func.dat";
const std::string GAL_COLOR_GMR_DENSITY = "/mount/sirocco1/imw2293/GROUP_CAT/SelfCalGroupFinder/py/parameters/bgs_y3/gal_g-r_density_func.dat";
  
// TODO Let's replace these with Boost functions soon
double qromo(double (*func)(double), double a, double b, double (*choose)(double(*)(double), double, double, int));
double midpnt(double (*func)(double), double a, double b, int n);


// LOGGING
#define LOG_VERBOSE(...)  if (!SILENT && VERBOSE) fprintf(stderr, __VA_ARGS__)
#define LOG_INFO(...)  if (!SILENT) fprintf(stderr, __VA_ARGS__)
#define LOG_PERF(...)  if (!SILENT) fprintf(stderr, __VA_ARGS__)
#define LOG_WARN(...)  if (!SILENT) fprintf(stderr, __VA_ARGS__)
#define LOG_ERROR(...) fprintf(stderr, __VA_ARGS__)

double gsl_spline_eval_extrap(const gsl_spline *spline, const double *x, const double *y, int n, double xq, gsl_interp_accel *acc) {
    // Constant extrapolation below range - this only happens due to floating point inaccuracy issues
    //std::cout << "Extrapolating spline at xq=" << xq << " with range [" << x[0] << ", " << x[n-1] << "]\n";
    if (xq < x[0]) {
        LOG_VERBOSE("xq (%.3f) is below the range. Extrapolating at left edge.\n", xq);
        return gsl_spline_eval(spline, x[0], acc);
    } else if (xq > x[n-1]) {
        LOG_VERBOSE("xq (%.3f) is above the range. Extrapolating at right edge.\n", xq);
        return gsl_spline_eval(spline, x[n-1], acc);
    } else {
        return gsl_spline_eval(spline, xq, acc);
    }
}


#define OMEGA_M 0.315192  // In past, 0.25 was used.
#define OMEGA_L (1 - OMEGA_M)
#define PI 3.141592741
#define RHO_CRIT 2.775E+11
#define SPEED_OF_LIGHT 299792.458 // in km/s
#define c_on_H0 2997.92458 // For H0 = 100 h km/s/Mpc
#define BIG_G 4.304E-9 /* BIG G in units of (km/s)^2*Mpc/M_sol */
#define RT2PI 2.50662827

/*
double func_dr1(double z) {
    return pow(OMEGA_M * (1 + z) * (1 + z) * (1 + z) + OMEGA_L, -0.5);
}

double distance_redshift(double z) {
    double x;
    if (z <= 0)
        return 0;
    x = c_on_H0 * qromo(func_dr1, 0.0, z, midpnt);
    return x;
}
*/

int filesize(const std::string &filename) {
    int size = 0;
    std::string line;
    std::ifstream count_ifs(filename);
    if (!count_ifs.is_open()) {
      std::cerr << "ERROR opening " << filename << std::endl;
      exit(ENOENT);
    }
    while (std::getline(count_ifs, line)) ++size;
    return size;
}

std::ifstream openfile(const std::string &filename) {
    std::ifstream ifs(filename);
    if (!ifs.is_open()) {
      std::cerr << "ERROR opening " << filename << std::endl;
      exit(ENOENT);
    }
    return ifs;
}


class TabulatedDensityFunction {

public:
    TabulatedDensityFunction(const std::string& filename) {
        std::ifstream fp = openfile(filename);
        n = filesize(filename);
        px.resize(n);
        rho.resize(n);
        for (int i = 0; i < n; i++)
            fp >> px[i] >> rho[i];
        fp.close();
        LOG_VERBOSE("Loaded density function with %d points from %s\n", n, filename.c_str());    }

    double getLeftEdge() {
        return px[0];
    }
    double getRightEdge() {
        return px[n - 1];
    }

    int getSize() const {
        return n;
    }
    const std::vector<double>& getX() const {
        return px;
    }
    const std::vector<double>& getDensity() const {
        return rho;
    }

private:
    int n = 0;
    std::vector<double> px; // x values for the density function
    std::vector<double> rho; // Density values at each x. Units are always assumed to be number/(Mpc/h)^3
};


/**
 * A cumulative density function for abundance matching built from a TabulatedDensityFunction. 
 */
class AMCDF {

public:

    AMCDF(const std::string& filename) {

        TabulatedDensityFunction pdf(filename);
        const int cnt = pdf.getSize();
        const std::vector<double> &px = pdf.getX();  // x values for the density function
        const std::vector<double> &rho = pdf.getDensity();  // density at each bin 

        // We calculate the CDF based on the PDF by directly integrating from right to left (end of file backwards) 
        // using the trapezoid rule, so the bin width is (x[i+1] - x[i]) and the density is rho[i] at the bin center.

        if (cnt < 100) {
            LOG_WARN("Density function %s has only %d points, which may be too few for accurate CDF integration.\n", filename.c_str(), cnt);
        }
        // TODO more sanity checks here. Maybe on the density numbers
        // In the group finder the halo mass function is in different units

        cx.resize(cnt);
        logCumulativeDensity.resize(cnt);

        double sum = 0.0;
        //double min = -80.0; // Minimum log cumulative density to avoid issues with log(0)
        //double dmin = 0.01; // This is redundant
        cx[cnt - 1] = px[cnt - 1];
        logCumulativeDensity[cnt - 1] = -80.0;; // assume n(> last point) = 0
        //min += dmin;
        for (int i = cnt - 2; i >= 0; i--) {
            double dx = abs(px[i + 1] - px[i]);  // uneven bins ok, descending or ascending order ok
            double area = 0.5 * (rho[i] + rho[i + 1]) * dx; // trapezoid area
            // Avoid issues with zero density bins. 
            // Small enough to not affect results but big enough to prevent the log(sum) from being the same value in multiple bins.
            // The spline fitting needs this. 
            // TODO Depending on the exact meaning of double on other hardware this might require tweaking.
            if (area < 1e-14) area = 1e-14; 
            sum += area;
            cx[i] = px[i];
            logCumulativeDensity[i] = log(sum);
            //min += dmin; // cannot repeat values for inversion
        }

        
        // Fit spline to the cumulative density function
        acc = gsl_interp_accel_alloc();
        spline = gsl_spline_alloc(gsl_interp_steffen, cnt);
        if (!spline || !acc) {
            LOG_ERROR("Failed to allocate GSL spline for CDF\n");
            exit(1);
        }
        // Reverse order if needed because GSL splines must be monotonically increasing in x
        if (cx[0] > cx[cnt - 1]) {
            std::reverse(cx.begin(), cx.end());
            std::reverse(logCumulativeDensity.begin(), logCumulativeDensity.end());
        }
        int status = gsl_spline_init(spline, cx.data(), logCumulativeDensity.data(), cnt);
        if (status) {
            LOG_ERROR("GSL spline init failed for CDF: %s\n", gsl_strerror(status));
            exit(1);
        }

        
        // Print out the logCumulativeDensity and cx values for debugging
        LOG_VERBOSE("Cumulative density function values for %s:\n", filename.c_str());
        for (int i = 0; i < cnt; i++) {
            LOG_VERBOSE("  x: %.4f, logCumulativeDensity: %.6f\n", cx[i], logCumulativeDensity[i]);
        }       

        
        // Now let's fit another spline for the inverse (reverse if needed)
        cx_inv.resize(cnt);
        logCumulativeDensity_inv.resize(cnt);
        bool reverse = logCumulativeDensity[0] > logCumulativeDensity[cnt - 1]; 
        for (int i = 0; i < cnt; i++) {
            cx_inv[i] = reverse ? cx[cnt - 1 - i] : cx[i];
            logCumulativeDensity_inv[i] = reverse ? logCumulativeDensity[cnt - 1 - i] : logCumulativeDensity[i];
        }

        inv_acc = gsl_interp_accel_alloc();
        inv_spline = gsl_spline_alloc(gsl_interp_steffen, cnt);
        if (!inv_spline || !inv_acc) {
            LOG_ERROR("Failed to allocate GSL spline for CDF inverse\n");
            exit(ENOMEM);
        }
        status = gsl_spline_init(inv_spline, logCumulativeDensity_inv.data(), cx_inv.data(), cnt);
        if (status) {
            LOG_ERROR("GSL spline init failed for CDF inverse: %s\n", gsl_strerror(status));
            exit(ENOMEM);
        }

        // Print out the logCumulativeDensity and cx values for debugging
        LOG_VERBOSE("Inverted cumulative density function values for %s:\n", filename.c_str());
        for (int i = 0; i < cnt; i++) {
            LOG_VERBOSE("  logCumulativeDensity: %.6f  x: %.4f, \n", logCumulativeDensity_inv[i], cx_inv[i]);
        }    

        LOG_VERBOSE("Built cumulative density spline for %s with %d points.\n", filename.c_str(), cnt);
    }

    double getLeftEdge() const {
        return cx[0];
    }
    double getRightEdge() const {
        return cx[cx.size() - 1];
    }

    /**
     * Evaluate the cumulative density function n(>value) at the given value using spline interpolation.
     * 
     * Returns the cumulative density n(>value) in units of number/(Mpc/h)^3. 
     */
    double eval(double value) {
        // Spline works on log CDF, so exponentiate the result to get n(>value)
        return exp(gsl_spline_eval_extrap(spline, cx.data(), logCumulativeDensity.data(), cx.size(), value, acc));
    }

    double eval_inverse(double cumulativeDensity) {
        // Given a cumulative density, return the corresponding property value using the inverse spline
        return gsl_spline_eval_extrap(inv_spline, logCumulativeDensity_inv.data(), cx_inv.data(), cx_inv.size(), log(cumulativeDensity), inv_acc);
    }

//private:
    std::vector<double> cx; // x values for the cumulative density function (same as PDF x values)
    std::vector<double> logCumulativeDensity; // log cumulative number density for each
    gsl_interp_accel* acc;
    gsl_spline* spline;

    std::vector<double> cx_inv; // x values for the cumulative density function (same as PDF x values)
    std::vector<double> logCumulativeDensity_inv; // log cumulative number density for each
    gsl_interp_accel* inv_acc;
    gsl_spline* inv_spline;

};

#define NZBIN 200
class AbundanceMatchingManager {
    public:

    /**
     * Match a provided density (which should be a running cumulative density) 
     * to a AMCDF to acquire a property value from the AMCDF.
     * 
     * Units must be consistent and are usually [1 / (Mpc/h)^3].
     * 
     * Returns the property value at the matched cumulative density.
     */
    virtual double match(double density) = 0;

    /** 
     * Matching for flux-limited mode - tracks the density seperately for overlapping redshift bins.
     * 
     * Using a vmax correction for galaxies that can't make it to the end of the redshift bin.
     */
    double match_in_zbins(double z, double vmax)  {
        /*
        int i, iz;
        double rlo, rhi, dz, dzmin, vv;

        // if first call, get the volume in each dz bin
        if (needs_setup) {
            for (i = 0; i < NZBIN; ++i) {
                zlo[i] = i * 1. / NZBIN;
                zhi[i] = zlo[i] + 0.05;
                if (i == 0)
                    rlo = 0;
                else
                    rlo = distance_redshift(zlo[i]);
                rhi = distance_redshift(zhi[i]);
                volume[i] = 4. / 3. * PI * (rhi * rhi * rhi - rlo * rlo * rlo) * FRAC_AREA;
                vhi[i] = 4. / 3. * PI * rhi * rhi * rhi * FRAC_AREA;
                vlo[i] = 4. / 3. * PI * rlo * rlo * rlo * FRAC_AREA;
            }
            needs_setup = false;
        }

        if (z < 0 || isnan(z) || !isfinite(z)) {
            LOG_ERROR("ERROR: invalid redshift %f\n", z);
            assert(false);
        }

        if (z > 100)
        {
            for (i = 0; i < NZBIN; ++i)
            if (negcnt[i])
                fprintf(stderr, "%d %f %d\n", i, zhi[i] - 0.025, negcnt[i]);
            return 0;
        }

        // Determine what bins this galaxy belong to
        // TODO this can definitely be optimized to not have a loop NZBIN times for each galaxy.
        dzmin = 1;
        for (i = 0; i < NZBIN; ++i) {
            if (z >= zlo[i] && z < zhi[i]) {
            //fprintf(stderr, "Matched z = %f to bin %d\n", z, i);
            if (vmax > vhi[i])
                vv = volume[i];
            else
                vv = vmax - vlo[i];
            if (vv < 0)
                vv = volume[i];
            negcnt[i]++;
            zcnt[i] += 1 / vv;

            if (vv < 0.0)
                LOG_ERROR("vmax = %e.  %e %e %e %e %e %e\n", vmax, vlo[i], vhi[i], zlo[i], zhi[i], z, zcnt[i]);
            }
            dz = fabs(z - (zhi[i] + zlo[i]) / 2);
            if (dz < dzmin) {
            dzmin = dz;
            iz = i;
            }
        }

        return match(zcnt[iz]);
         */
        throw std::runtime_error("match_in_zbins is not implemented yet. Volume-limited mode only for now.");
        //return 0.0; // placeholder for now
    }

    /**
     * Reset SHAM counters for the z-bins. Not needed in volume-limited mode.
     */
    void reset() {
        for (int i = 0; i < NZBIN; ++i)
        zcnt[i] = negcnt[i] = 0;
    }

    private:
      bool needs_setup = true;
      int negcnt[NZBIN];
      double zcnt[NZBIN];
      double volume[NZBIN], zlo[NZBIN], zhi[NZBIN], vhi[NZBIN], vlo[NZBIN];
};
  

class GalaxyMagMatcher : public AbundanceMatchingManager {
public:
    static GalaxyMagMatcher& get() {
        static GalaxyMagMatcher inst;
        return inst;
    }

    double match(double density) {
        if (density > 0.05) {
            LOG_WARN("Warning: Matching density %.3e is too dense for current density functions.\n", density);
        }

        // Best approach is to use use an inverted spline of the CDF to directly get the property
        return mag_cdf->eval_inverse(density); 

        // Boost root finding approach. Slower
        /*
        double min = std::min(mag_cdf->getLeftEdge(), mag_cdf->getRightEdge());
        double max = std::max(mag_cdf->getLeftEdge(), mag_cdf->getRightEdge());

        auto func = [this, density](double value) {
            double matched_dens = mag_cdf->eval(value);
            //std::cout << "Evaluated density at value " << value << " is " << matched_dens << " for target density " << density << std::endl;
            return matched_dens - density;        
        };

        int bits = 18;  // TODO what should it be?  1 / 2^bits tolerance I think is meaning
        eps_tolerance<double> tol(bits);
        boost::uintmax_t max_iter = 100;
        //std::cout << "Matching galaxy density " << density << " to magnitude with bounds [" << min << ", " << max << "]..." << std::endl;

        // std::pair<double, double> result = bracket_and_solve_root(func, guess, factor, rising, tol);
        std::pair<double, double> result = toms748_solve(func, min, max, tol, max_iter);
        //std::cout << "Matched galaxy density " << density << " to magnitude " << result.first << " with function value " << func(result.first) << std::endl;
        return result.first; // root is between result.first and result.second, but they should be very close given the tolerance
        */

        // Old Group Finder appoach. Slower
        //return exp(zbrent(func, min, max, 1.0E-5, density));
    }
    
    //~GalaxyMagMatcher() {
    //    delete mag_cdf;
    //}
private:
    AMCDF *mag_cdf;

    GalaxyMagMatcher() {
        // Load the galaxy magnitude density function and build the spline for abundance matching
        mag_cdf = new AMCDF(GAL_MAGR_DENSITY);
    }


};


