/***********************************************************************
 * "latsham.hpp"
 * This file is a header-only library containing the core interface for
 * latent abundance matching models.
 * 
 * This file is part of the latsham package. 
 * Copyright (c) 2026 Ian Williams under the MIT License.
 ***********************************************************************/

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
#include <algorithm>
#include <sstream>
//#include <boost/math/tools/roots.hpp>
#include <limits>
#include <stdlib.h>
#include <map>
#include <memory>

#define SILENT 0
#define VERBOSE 0
#define PERF 1

#define NORMALIZE_DENSITY_FUNCTIONS 1

//using boost::math::tools::eps_tolerance;
//using boost::math::tools::toms748_solve;

// TODO Let's replace these with Boost functions soon
double qromo(double (*func)(double), double a, double b, double (*choose)(double(*)(double), double, double, int));
double midpnt(double (*func)(double), double a, double b, int n);

/**
 * Get a linearly interpolated value from a set of x and y values.
 *
 * @param x The x values (must be sorted in ascending order).
 * @param y The y values corresponding to the x values.
 * @param xvalue The x value at which to interpolate.
 * @return The interpolated y value.
 */
double lin_interp(const std::vector<double>& x, const std::vector<double>& y, double xvalue) {
    auto it = std::lower_bound(x.begin(), x.end(), xvalue);
    if (it == x.begin()) return y.front();
    if (it == x.end())   return y.back();
    auto i1 = it;
    auto i0 = std::prev(it);
    double t = (xvalue - *i0) / (*i1 - *i0);
    return y[i0 - x.begin()] + t * (y[i1 - x.begin()] - y[i0 - x.begin()]); // Linear interpolation between the two nearest points
}


// LOGGING
#define LOG_VERBOSE(...)  if (!SILENT && VERBOSE) fprintf(stderr, __VA_ARGS__)
#define LOG_INFO(...)  if (!SILENT) fprintf(stderr, __VA_ARGS__)
#define LOG_PERF(...)  if (!SILENT && PERF) fprintf(stderr, __VA_ARGS__)
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

inline double absmagToLogLgal(double mag) {
    return (4.65 - mag) / 2.5; // log10(L/Lsun)
}

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

/**
 * A simple interface for a model that has properties and latent properties that can be accessed by index.
 * The ordering of the properties must match whatever is in the associated LatentModel.
 */
struct Transformable {
    virtual double getProperty(int num) = 0;
    virtual void setProperty(int num, double value) = 0;
    virtual double getLatentProperty(int num) = 0;
    virtual void setLatentProperty(int num, double value) = 0;
};


/**
 * A tabulated density function loaded from a file. 
 * The file should have two columns: x and density(x). 
 * The density is assumed to be in units of number/(Mpc/h)^3. 
 */
class TabulatedDensityFunction {

public:
    TabulatedDensityFunction(const std::string& filename) {
        std::ifstream fp = openfile(filename);
        n = filesize(filename);
        px.resize(n);
        rho.resize(n);
        for (int i = 0; i < n; ++i)
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

        std::vector<double> propValues; // x values (property values) for the CDF
        std::vector<double> logCumulativeDensity; // y values (log cumulative number density) for each x

        // We calculate the CDF based on the PDF by directly integrating from right to left (end of file backwards) 
        // using the trapezoid rule, so the bin width is (x[i+1] - x[i]) and the density is rho[i] at the bin center.

        if (cnt < 100) {
            LOG_WARN("Density function %s has only %d points, which may be too few for accurate CDF integration.\n", filename.c_str(), cnt);
        }
        // TODO more sanity checks here. Maybe on the density numbers
        // In the group finder the halo mass function is in different units

        propValues.resize(cnt);
        logCumulativeDensity.resize(cnt);

        double sum = 0.0;
        //double min = -80.0; // Minimum log cumulative density to avoid issues with log(0)
        //double dmin = 0.01; // This is redundant
        propValues[cnt - 1] = px[cnt - 1];
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
            propValues[i] = px[i];
            logCumulativeDensity[i] = log(sum);
            //min += dmin; // cannot repeat values for inversion
        }

        // Reverse order if needed because GSL splines must be monotonically increasing in x
        if (propValues[0] > propValues[cnt - 1]) {
            std::reverse(propValues.begin(), propValues.end());
            std::reverse(logCumulativeDensity.begin(), logCumulativeDensity.end());
        }

        
        // Print out the logCumulativeDensity and cx values for debugging
        LOG_VERBOSE("Cumulative density function values for %s:\n", filename.c_str());
        for (int i = 0; i < cnt; ++i) {
            LOG_VERBOSE("  x: %.4f, logCumulativeDensity: %.6f\n", propValues[i], logCumulativeDensity[i]);
        }       

        
        // Now let's fit a spline for prop value as a function of log cumulative density (reverse if needed)
        propValues_inv.resize(cnt);
        logCumulativeDensity_inv.resize(cnt);
        bool reverse = logCumulativeDensity[0] > logCumulativeDensity[cnt - 1]; 
        for (int i = 0; i < cnt; ++i) {
            propValues_inv[i] = reverse ? propValues[cnt - 1 - i] : propValues[i];
            logCumulativeDensity_inv[i] = reverse ? logCumulativeDensity[cnt - 1 - i] : logCumulativeDensity[i];
        }

        inv_acc = gsl_interp_accel_alloc();
        inv_spline = gsl_spline_alloc(gsl_interp_steffen, cnt);
        if (!inv_spline || !inv_acc) {
            LOG_ERROR("Failed to allocate GSL spline for CDF inverse\n");
            exit(ENOMEM);
        }
        auto status = gsl_spline_init(inv_spline, logCumulativeDensity_inv.data(), propValues_inv.data(), cnt);
        if (status) {
            LOG_ERROR("GSL spline init failed for CDF inverse: %s\n", gsl_strerror(status));
            exit(ENOMEM);
        }

        // Print out the logCumulativeDensity and cx values for debugging
        LOG_VERBOSE("Inverted cumulative density function values for %s:\n", filename.c_str());
        for (int i = 0; i < cnt; ++i) {
            LOG_VERBOSE("  logCumulativeDensity: %.6f  x: %.4f, \n", logCumulativeDensity_inv[i], propValues_inv[i]);
        }    

        LOG_VERBOSE("Built cumulative density spline for %s with %d points.\n", filename.c_str(), cnt);
    }

    double getLeftEdge() const {
        return this->propValues_inv.front();
    }
    double getRightEdge() const {
        return this->propValues_inv.back();
    }

    /**
     * Get a property value for a given cumulative density value.
     * If NORMALIZE_DENSITY_FUNCTIONS is enabled, it should be a normalized cumulative density value between 0 and 1.
     * Otherwise it should be in physical units matching the tabulated density function file (usually Mpc^-3 h^3).
     */
    double eval_inverse(double cumulativeDensity) {
        // Given a cumulative density, return the corresponding property value using the inverse spline
        //return gsl_spline_eval_extrap(inv_spline, logCumulativeDensity_inv.data(), propValues_inv.data(), propValues_inv.size(), log(cumulativeDensity), inv_acc);

        // Instead let's try linear interpolation for now
        // TODO is this wrong because it's in log space?
        return lin_interp(logCumulativeDensity_inv, propValues_inv, log(cumulativeDensity));
    }

    /**
     * @brief Normalize the cumulative density function so that the maximum cumulative density is 1.0
     */
    void normalize() {
        
        double max_density = exp(logCumulativeDensity_inv.back());
        if (max_density <= 0) {
            LOG_ERROR("Cannot normalize CDF: maximum cumulative density is non-positive (%.6e)\n", max_density);
            return;
        }
        for (double &val : this->logCumulativeDensity_inv) {
            val -= log(max_density);
        }
    }

//private:
    std::vector<double> propValues_inv; // x values (property values) for the inverse CDF
    std::vector<double> logCumulativeDensity_inv;  // log cumulative number density for each, reversed for inversion
    gsl_interp_accel* inv_acc;
    gsl_spline* inv_spline;

};

#define NZBIN 200

/**
 * Performs abundance matching against a cumulative density function. 
 * See AMCDF / TabulatedDensityFunction. See AbundanceMatchingRegistry below for a central place that
 * creates/caches instances of this class, one per density function file.
 */
class AbundanceMatchingManager {
    public:

    AbundanceMatchingManager(const std::string& filename) : cdf(filename) {}

    /**
     * Match a provided density (which should be a running cumulative density) 
     * to a AMCDF to acquire a property value from the AMCDF.
     * 
     * Units must be consistent and are usually [1 / (Mpc/h)^3].
     * 
     * Returns the property value at the matched cumulative density.
     */
    double match(double density) {
        if (NORMALIZE_DENSITY_FUNCTIONS!=1 && density > 0.05) {
            LOG_WARN("Warning: Matching density %.3e is too dense for current density functions.\n", density);
        }

        return cdf.eval_inverse(density);
    }

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

    void normalize() {
        if (this->normalized) return;
        
        this->cdf.normalize();
        // TOOD if we need the zbins version, think about this
        this->normalized = true;
    }

    private:
      AMCDF cdf;

      bool needs_setup = true;
      bool normalized = false;
      int negcnt[NZBIN];
      double zcnt[NZBIN];
      double volume[NZBIN], zlo[NZBIN], zhi[NZBIN], vhi[NZBIN], vlo[NZBIN];
};

/**
 * Central repository of AbundanceMatchingManager instances, keyed by the density function
 * file they were built from. Lazily constructs and caches a manager per unique file so the
 * AMCDF spline construction only ever happens once per file.
 */
class AbundanceMatchingRegistry {
public:
    static AbundanceMatchingRegistry& get() {
        static AbundanceMatchingRegistry inst;
        return inst;
    }

    AbundanceMatchingManager& getMatcher(const std::string& filename) {
        auto it = matchers.find(filename);
        if (it == matchers.end()) {
            it = matchers.emplace(filename, std::make_unique<AbundanceMatchingManager>(filename)).first;
        }
        return *it->second;
    }

private:
    AbundanceMatchingRegistry() = default;
    std::map<std::string, std::unique_ptr<AbundanceMatchingManager>> matchers;
};


// Interface for latent models (PCA/ICA) that can transform data to and from a latent space.
class ILatentModel {
public:
    virtual void forward_transform(Transformable &obj) = 0;
    virtual void inverse_transform(Transformable &obj) = 0;
    virtual void print_model() = 0;
    virtual int getFeatureCount() = 0;
    virtual void load() = 0;
};

// Class that holds a basic model as matrices for transforming data between an original feature space and a latent space.
// The matrix loaded from disk should contain everything needed to transform folded into a single matrix for this to work as intended.
class LatentModel : public ILatentModel {
public:
    LatentModel(std::string model_file) : model_file(model_file) {}

    ~LatentModel() = default;

    // Forward: IC[i] = sum_j (x[j] - mu[j]) * A[j][i]
    void forward_transform(Transformable &obj) {
        if (!loaded) load();
        for (int i = 0; i < NFEAT; ++i) {
            double ic = 0;
            for (int j = 0; j < NFEAT; ++j)
                ic += (obj.getProperty(j) - mu[j]) * A[j][i];
            obj.setLatentProperty(i, ic);
        }
    }

    // Inverse: x[j] = sum_i IC[i] * B[i][j] + mu[j]
    void inverse_transform(Transformable &obj) {
        if (!loaded) load();
        for (int j = 0; j < NFEAT; ++j) {
            double x = mu[j];
            for (int i = 0; i < NFEAT; ++i)
                x += obj.getLatentProperty(i) * B[i][j];
            obj.setProperty(j, x);
        }
    }

    void print_model() {
        if (!loaded) load();
        std::cout << "mean:";
        for (int i = 0; i < NFEAT; ++i) std::cout << " " << mu[i];
        std::cout << "\nforward_matrix (A):\n";
        for (int i = 0; i < NFEAT; ++i) {
            for (int j = 0; j < NFEAT; ++j) std::cout << A[i][j] << " ";
            std::cout << "\n";
        }
        std::cout << "inverse_matrix (B):\n";
        for (int i = 0; i < NFEAT; ++i) {
            for (int j = 0; j < NFEAT; ++j) std::cout << B[i][j] << " ";
            std::cout << "\n";  
        }
    }

    int getFeatureCount() {
        return NFEAT;
    }

private:
    std::string model_file;
    int NFEAT = 0;
    std::vector<double> mu;
    std::vector<std::vector<double>> A; // A[feature][component]
    std::vector<std::vector<double>> B; // B[component][feature]
    bool loaded = false;

    void load() {
        std::ifstream fp(model_file);
        if (!fp) { fprintf(stderr, "Could not open %s\n", model_file.c_str()); exit(1); }
        std::string line;
        int block = 0, j = 0;
        while (std::getline(fp, line)) {
            if (line.find("# n_features:") != std::string::npos) {
                std::istringstream ss(line.substr(13));
                ss >> NFEAT;
                mu.resize(NFEAT);
                A.assign(NFEAT, std::vector<double>(NFEAT));
                B.assign(NFEAT, std::vector<double>(NFEAT));
                continue;
            }
            if (line.empty() || line[0] == '#') continue;
            if (NFEAT <= 0) { LOG_ERROR("n_features must appear before data blocks\n"); exit(1); }
            std::istringstream ss(line);
            switch (block) {
                case 0: for (int i = 0; i < NFEAT; ++i) ss >> mu[i]; break;
                case 1: for (int i = 0; i < NFEAT; ++i) ss >> A[j][i]; ++j; break;
                case 2: for (int i = 0; i < NFEAT; ++i) ss >> B[j][i]; ++j; break;
            }
            if (j == NFEAT) j = 0;
            if (j == 0) ++block;
        }
        loaded = true;
        LOG_VERBOSE("Loaded latent model from %s, NFEAT=%d\n", model_file.c_str(), NFEAT);
    }
};

