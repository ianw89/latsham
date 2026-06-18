/***********************************************************************
 * This file is part of the latsham package. 
 * Copyright (c) 2026 Ian Williams under the MIT License.
 ***********************************************************************/

#pragma once
#include <vector>
#include <algorithm>
#include <string>
#include <stdexcept>
#include <random>
#include "latsham.hpp"

// This file is already cut to central halos only; no subhalos. That's the 'C' in the name.
// It is also cut to M200b > 8E9 Msun
const std::string HALO_CATALOG = "/mount/sirocco1/imw2293/GROUP_CAT/DATA/POPMOCK/smdpl_z0.19717.M8E9.C.h5";
double constexpr BOX_SIZE = 400.0; // Mpc/h, from the simulation specs
double constexpr SIM_VOLUME = BOX_SIZE * BOX_SIZE * BOX_SIZE; // (Mpc/h)^3
double constexpr HALO_DENSITY_INCREMENT = 1.0 / SIM_VOLUME; 
// ('index', 'ID', 'upid', 'M200b', 'Mpeak', 'mvir', 'rvir', 'rs', 'vmax', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'Halfmass_Scale', 'Spin', 'c', 'LOGMHALO', 'NOISE')

const std::string GAL_MAGR_DENSITY = "/mount/sirocco1/imw2293/GROUP_CAT/SelfCalGroupFinder/py/parameters/bgs_y3/gal_abs_mag_r_k_best_density_func.dat";
const std::string GAL_COLOR_GMR_DENSITY = "/mount/sirocco1/imw2293/GROUP_CAT/SelfCalGroupFinder/py/parameters/bgs_y3/gal_g-r_density_func.dat";
const std::string GAL_2P_ICA1_DENSITY = "/mount/sirocco1/imw2293/GROUP_CAT/SelfCalGroupFinder/py/parameters/bgs_y3/gal_2p_ica1_density_func.dat";
const std::string GAL_2P_ICA2_DENSITY = "/mount/sirocco1/imw2293/GROUP_CAT/SelfCalGroupFinder/py/parameters/bgs_y3/gal_2p_ica2_density_func.dat";
  
const std::string HALO_4P_LATENT_MODEL_TEXT_FILE = "/mount/sirocco1/imw2293/GROUP_CAT/SelfCalGroupFinder/py/parameters/bgs_y3/halo_ica_model.txt";
const std::string HALO_2P_LATENT_MODEL_TEXT_FILE = "/mount/sirocco1/imw2293/GROUP_CAT/SelfCalGroupFinder/py/parameters/bgs_y3/halo_2p_model.txt";
const std::string GAL_2P_LATENT_MODEL_TEXT_FILE = "/mount/sirocco1/imw2293/GROUP_CAT/SelfCalGroupFinder/py/parameters/bgs_y3/gal_2p_model.txt";


struct Galaxy : Transformable {
    double abs_mag_r; // absolute magnitude in the r band, k-corrected to z=0.1, with h=1.0
    double color_g_r;   // g-r color (mag)
    double lat1, lat2;

    double getProperty(int num) override {
        switch (num) {
            case 0: return abs_mag_r;
            case 1: return color_g_r;
            default: throw std::out_of_range("Invalid property number");
        }
    }
    void setProperty(int num, double value) override {
        switch (num) {
            case 0: abs_mag_r = value; break;
            case 1: color_g_r = value; break;
            default: throw std::out_of_range("Invalid property number");
        }
    }
    double getLatentProperty(int num) override {
        switch (num) {
            case 0: return lat1;
            case 1: return lat2;
            default: throw std::out_of_range("Invalid latent property number");
        }
    }
    void setLatentProperty(int num, double value) override {
        switch (num) {
            case 0: lat1 = value; break;
            case 1: lat2 = value; break;
            default: throw std::out_of_range("Invalid latent property number");
        }
    }
};

struct Halo : Transformable {
    double x,y,z;
    double logmhalo; // log10 of M200b
    double halfmass_scale;
    double c;
    double spin;
    double lat1, lat2;

    // temp properties used for ranking
    double rank;

    Galaxy galaxy; // The galaxy that occupies this halo after abundance matching

    double getProperty(int num) override {
        switch (num) {
            case 0: return logmhalo;
            case 1: return halfmass_scale;
            default: throw std::out_of_range("Invalid property number");
        }
    }
    void setProperty(int num, double value) override {
        switch (num) {
            case 0: logmhalo = value; break;
            case 1: halfmass_scale = value; break;
            default: throw std::out_of_range("Invalid property number");
        }
    }
    double getLatentProperty(int num) override {
        switch (num) {
            case 0: return lat1;
            case 1: return lat2;
            default: throw std::out_of_range("Invalid latent property number");
        }
    }
    void setLatentProperty(int num, double value) override {
        switch (num) {
            case 0: lat1 = value; break;
            case 1: lat2 = value; break;
            default: throw std::out_of_range("Invalid latent property number");
        }
    }
};

void sortByRank(std::vector<Halo>& halos) {
    std::sort(halos.begin(), halos.end(),
        [](const Halo& a, const Halo& b) { 
            return a.rank > b.rank; 
        });
}

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

class GalaxyColorGMRMatcher : public AbundanceMatchingManager {
public:
    static GalaxyColorGMRMatcher& get() {
        static GalaxyColorGMRMatcher inst;
        return inst;
    }

    double match(double density) {
        if (density > 0.05) {
            LOG_WARN("Warning: Matching density %.3e is too dense for current density functions.\n", density);
        }
        return color_cdf->eval_inverse(density); 
    }

private:
    AMCDF *color_cdf;

    GalaxyColorGMRMatcher() {
        color_cdf = new AMCDF(GAL_COLOR_GMR_DENSITY);
    }
};

class Galaxy2P_ICA1Matcher : public AbundanceMatchingManager {
public:
    static Galaxy2P_ICA1Matcher& get() {
        static Galaxy2P_ICA1Matcher inst;
        return inst;
    }

    double match(double density) {
        if (density > 0.05) {
            LOG_WARN("Warning: Matching density %.3e is too dense for current density functions.\n", density);
        }

        return cdf->eval_inverse(density); 
    }
    
private:
    AMCDF *cdf;

    Galaxy2P_ICA1Matcher() {
        cdf = new AMCDF(GAL_2P_ICA1_DENSITY);
    }
};

class Galaxy2P_ICA2Matcher : public AbundanceMatchingManager {
public:
    static Galaxy2P_ICA2Matcher& get() {
        static Galaxy2P_ICA2Matcher inst;
        return inst;
    }

    double match(double density) {
        if (density > 0.05) {
            LOG_WARN("Warning: Matching density %.3e is too dense for current density functions.\n", density);
        }

        return cdf->eval_inverse(density); 
    }
    
private:
    AMCDF *cdf;

    Galaxy2P_ICA2Matcher() {
        cdf = new AMCDF(GAL_2P_ICA2_DENSITY);
    }
};

class ConnectionModel {
public:
    virtual void setParamsFromList(std::vector<double> params) = 0;
    virtual bool match(std::vector<Halo>& halos) = 0;
    virtual bool writeMocks(std::vector<Halo>& halos) = 0;
    virtual void load() = 0;
    
    std::unique_ptr<LatentModel> galaxyModel;
    std::unique_ptr<LatentModel> haloModel;
};

class ConnectionModel_1 : public ConnectionModel {
public:
    double scatter;

    void setParamsFromList(std::vector<double> params) {
        if (params.size() != 1) {
            LOG_ERROR("Expected 1 parameter for Model1, got %zu\n", params.size());
            throw std::runtime_error("Invalid number of parameters for Model1");
        }
        scatter = params[0];
    }

    bool match(std::vector<Halo>& halos) override {

        // Lognormal scatter in luminosity at fixed halo mass
        std::random_device rd;
        std::mt19937 gen(rd());
        double mean = 0.0; 
        double std_dev = scatter;  
        std::normal_distribution<double> gaussian_dist(mean, std_dev);

        for (Halo &h : halos) {
            h.rank = h.logmhalo + gaussian_dist(gen);
        }
        sortByRank(halos);

        // Now abundance match to galaxy r band abs magnitude using
        double density = 0.0;
        for (Halo& h : halos) {
            density += HALO_DENSITY_INCREMENT; // cumulative density of halos above this mass
            h.galaxy.abs_mag_r = GalaxyMagMatcher::get().match(density);
        }

        return true;
    }    

    bool writeMocks(std::vector<Halo>& halos) override {
        // Write files out in the format that corrfun expects for wp calculation.
        std::string out_filename = "/mount/sirocco1/imw2293/GROUP_CAT/OUTPUT/LATSHAM/mock_M20.dat";
        std::ofstream outputFile(out_filename);

        // Write the central galaxies for this mock
        for (int i=0; i < halos.size(); ++i) {
            if (halos[i].galaxy.abs_mag_r < -21.0 || halos[i].galaxy.abs_mag_r > -20.0)
                continue;

            outputFile << halos[i].x << " " << halos[i].y << " " << halos[i].z << "\n";
        }

        return true;
        // No satellites for this study; we will compare to centrals-only clustering measurements.    
    }
};

class ConnectionModel_2P2P : public ConnectionModel {
public:

    void load() {
        this->galaxyModel = std::make_unique<LatentModel>(GAL_2P_LATENT_MODEL_TEXT_FILE);
        this->haloModel = std::make_unique<LatentModel>(HALO_2P_LATENT_MODEL_TEXT_FILE);
        loaded = true;
    }

    // Best:  [0.88693407 3.33345685] 
    //3.98728008 5.84725974
    // [-0.98, 0.198, 0.824, -0.567]
    // -0.9802, 0.824
    // [3.9827825  0.32284589]
    double p1 = cos(3.9827825);
    double p2 = sin(3.9827825);
    double p3 = cos(0.32284589);
    double p4 = sin(0.32284589);
    
    void setParamsFromList(std::vector<double> params) {
        if (params.size() != 4) {
            LOG_ERROR("Expected 4 parameters for Model_2P2P, got %zu\n", params.size());
            throw std::runtime_error("Invalid number of parameters for Model_2L");
        }
        p1 = params[0];
        p2 = params[1];
        p3 = params[2];
        p4 = params[3];
        LOG_INFO("Set model params: p1=%.3f, p2=%.3f, p3=%.3f, p4=%.3f\n", p1, p2, p3, p4);
    }

    bool match(std::vector<Halo>& halos) override {
        if (!loaded) load();

        // Perf note: I tried storing a lookup of density values, but somehow it was slower than recomputing it.
        // No gains from parallelizing the match(density) call.
        // The sorts take ~250ms each on howdy
        // Not much to improve here. Replacing mock writing with direct call to corrfunc from code is most important.

        //#pragma omp parallel for
        for (size_t i = 0; i < halos.size(); ++i) {
            halos[i].rank = p1 * halos[i].lat1 + p2 * halos[i].lat2;
        }
        sortByRank(halos);

        double density = 0.0;
        for (size_t i = 0; i < halos.size(); ++i) {
            density += HALO_DENSITY_INCREMENT; // cumulative density of halos above this mass
            halos[i].galaxy.setLatentProperty(0, Galaxy2P_ICA1Matcher::get().match(density));
        }

        //#pragma omp parallel for
        for (size_t i = 0; i < halos.size(); ++i) {
            halos[i].rank = p3 * halos[i].lat1 + p4 * halos[i].lat2;
        }
        sortByRank(halos);

        density = 0.0;
        for (size_t i = 0; i < halos.size(); ++i) {
            density += HALO_DENSITY_INCREMENT; // cumulative density of halos above this mass
            halos[i].galaxy.setLatentProperty(1, Galaxy2P_ICA2Matcher::get().match(density));
        }

        #pragma omp parallel for
        for (size_t i = 0; i < halos.size(); ++i) {
            this->galaxyModel->inverse_transform(halos[i].galaxy);
        }

        // Print results for the first 5 halos
        /*for (int i = 0; i < 5 && i < halos.size(); ++i) {
            const Halo& h = halos[i];
            LOG_INFO("Halo %d: logmhalo=%.3f, age=%.3f, color=%.3f, abs_mag_r=%.3f\n", i, 
                h.logmhalo, h.halfmass_scale, h.galaxy.color_g_r, h.galaxy.abs_mag_r);
        }*/

        // Re-abundance properties match onto itself to ensure original propety distributions are respected
        //#pragma omp parallel for
        for (size_t i = 0; i < halos.size(); ++i) {
            halos[i].rank = halos[i].galaxy.getProperty(0); // re-apply the same latent property values to force them to be self-consistent with the mag/color distribution after the inverse transform. This is a hack to deal with the higher-order residual correlation issue.
        }
        sortByRank(halos);
        density = 0.0;
        for (Halo &h : halos) {
            density += HALO_DENSITY_INCREMENT; 
            h.galaxy.setProperty(0, GalaxyMagMatcher::get().match(density));
        }
        
        //#pragma omp parallel for
        for (size_t i = 0; i < halos.size(); ++i) {
            halos[i].rank = halos[i].galaxy.getProperty(1); // re-apply the same latent property values to force them to be self-consistent with the mag/color distribution after the inverse transform. This is a hack to deal with the higher-order residual correlation issue.
        }
        sortByRank(halos);
        density = 0.0;
        for (Halo &h : halos) {
            density += HALO_DENSITY_INCREMENT; 
            h.galaxy.setProperty(1, GalaxyColorGMRMatcher::get().match(density));
        }

        return true;
    }    

    bool writeMocks(std::vector<Halo>& halos) override {
        // Write files out in the format that corrfun expects for wp calculation.
        std::string out_basename = "/mount/sirocco1/imw2293/GROUP_CAT/OUTPUT/LATSHAM/mock_%s_M%d.dat";

        std::vector<int> magbins = {-17, -18, -19, -20, -21, -22};
        std::vector<std::string> colors = {"red", "blue"};

        int successes = 0;
        #pragma omp parallel for collapse(2) reduction(+:successes)
        for (int i = 0; i < magbins.size(); ++i) {
            for (int j = 0; j < colors.size(); ++j) {
                auto magbin = magbins[i];
                auto color = colors[j];
                std::string out_filename = out_basename;
                out_filename = out_filename.replace(out_filename.find("%s"), 2, color);
                out_filename = out_filename.replace(out_filename.find("%d"), 2, std::to_string(abs(magbin)));
                LOG_VERBOSE("Writing mock for magbin %d, color %s to %s\n", magbin, color.c_str(), out_filename.c_str());
                successes += writeMockForBin(halos, magbin, color, out_filename);
            }
        }
        LOG_INFO("Successfully wrote %d mocks out of %d\n", successes, (int)(magbins.size() * colors.size()));
        return successes == (magbins.size() * colors.size());
    }

private:
    bool loaded = false;

    bool writeMockForBin(const std::vector<Halo>& halos, int magbin, const std::string& color, const std::string& out_filename) {
        std::ofstream outputFile(out_filename);
        double color_cut = 0.76; 
        int count = 0;

        // No satellites for this study; we will compare to centrals-only clustering measurements.    

        for (const Halo& h : halos) {
            if (h.galaxy.abs_mag_r < (magbin - 1) || h.galaxy.abs_mag_r > magbin)
                continue;
            if (color == "red" && h.galaxy.color_g_r < color_cut)
                continue;
            if (color == "blue" && h.galaxy.color_g_r >= color_cut)
                continue;

            ++count;
            outputFile << h.x << " " << h.y << " " << h.z << "\n";
        }
        LOG_VERBOSE("Wrote %d galaxies for magbin %d, color %s\n", count, magbin, color.c_str());
        return count > 0;
    }
};