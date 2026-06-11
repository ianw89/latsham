/***********************************************************************
 * This file is part of the latsham package. 
 * Copyright (c) 2026 Ian Williams under the MIT License.
 ***********************************************************************/

#pragma once

#include "latsham.hpp"

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

    // temp property used for ranking
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