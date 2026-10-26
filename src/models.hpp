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
#include <regex>
#include <set>
#include <map>
#include <filesystem>
#include <iostream>
#include "latsham.hpp"

// This file is already cut to central halos only; no subhalos. That's the 'C' in the name.
const std::string HALO_CATALOG = "/mount/sirocco1/imw2293/GROUP_CAT/DATA/POPMOCK/smdpl_z0.19717.M1E11.C.MagBin4.h5";
double constexpr BOX_SIZE = 400.0; // Mpc/h, from the simulation specs
double constexpr SIM_VOLUME = BOX_SIZE * BOX_SIZE * BOX_SIZE; // (Mpc/h)^3
double constexpr HALO_DENSITY_INCREMENT = 1.0 / SIM_VOLUME; 
// ('index', 'ID', 'upid', 'M200b', 'Mpeak', 'mvir', 'rvir', 'rs', 'vmax', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'Halfmass_Scale', 'Spin', 'c', 'LOGMHALO', 'NOISE')

const std::string GAL_MAGR_DENSITY = "/mount/sirocco1/imw2293/GROUP_CAT/SelfCalGroupFinder/py/parameters/bgs_y3/gal_abs_mag_r_k_best_density_func.dat";
const std::string GAL_COLOR_GMR_DENSITY = "/mount/sirocco1/imw2293/GROUP_CAT/SelfCalGroupFinder/py/parameters/bgs_y3/gal_g-r_density_func.dat";

const std::string GAL_MAGBIN4_MAGR_DENSITY = "/mount/sirocco1/imw2293/GROUP_CAT/SelfCalGroupFinder/py/parameters/bgs_y3/gal_magbin4_abs_mag_r_k_best_density_func.dat";
const std::string GAL_MAGBIN4_C9050_DENSITY = "/mount/sirocco1/imw2293/GROUP_CAT/SelfCalGroupFinder/py/parameters/bgs_y3/gal_magbin4_c9050_density_func.dat";
const std::string GAL_MAGBIN4_COLOR_GMR_DENSITY = "/mount/sirocco1/imw2293/GROUP_CAT/SelfCalGroupFinder/py/parameters/bgs_y3/gal_magbin4_g-r_density_func.dat";
const std::string GAL_MAGBIN4_DN4000_DENSITY = "/mount/sirocco1/imw2293/GROUP_CAT/SelfCalGroupFinder/py/parameters/bgs_y3/gal_magbin4_dn4000_model_density_func.dat";
const std::string GAL_MAGBIN4_LOGHALPHA_DENSITY = "/mount/sirocco1/imw2293/GROUP_CAT/SelfCalGroupFinder/py/parameters/bgs_y3/gal_magbin4_log_halpha_ew_density_func.dat";
const std::string GAL_MAGBIN4_LOGMSTAR_DENSITY = "/mount/sirocco1/imw2293/GROUP_CAT/SelfCalGroupFinder/py/parameters/bgs_y3/gal_magbin4_logmstar_density_func.dat";
const std::string GAL_MAGBIN4_HALFLIGHT_DENSITY = "/mount/sirocco1/imw2293/GROUP_CAT/SelfCalGroupFinder/py/parameters/bgs_y3/gal_magbin4_shape_r_kpc_density_func.dat";

std::string galDensityFilename(const std::string& latent_model_name, int component_index) {
    return "/mount/sirocco1/imw2293/GROUP_CAT/SelfCalGroupFinder/py/parameters/bgs_y3/gal_" + latent_model_name + "_pca" + std::to_string(component_index + 1) + "_density_func.dat";
}

const std::string GAL_2P_ICA1_DENSITY = "/mount/sirocco1/imw2293/GROUP_CAT/SelfCalGroupFinder/py/parameters/bgs_y3/gal_2p_ica1_density_func.dat";
const std::string GAL_2P_ICA2_DENSITY = "/mount/sirocco1/imw2293/GROUP_CAT/SelfCalGroupFinder/py/parameters/bgs_y3/gal_2p_ica2_density_func.dat";
  
const std::string HALO_4P_LATENT_MODEL_TEXT_FILE = "/mount/sirocco1/imw2293/GROUP_CAT/SelfCalGroupFinder/py/parameters/bgs_y3/halo_ica_model.txt";
const std::string HALO_2P_LATENT_MODEL_TEXT_FILE = "/mount/sirocco1/imw2293/GROUP_CAT/SelfCalGroupFinder/py/parameters/bgs_y3/halo_2p_model.txt";
const std::string GAL_2P_LATENT_MODEL_TEXT_FILE = "/mount/sirocco1/imw2293/GROUP_CAT/SelfCalGroupFinder/py/parameters/bgs_y3/gal_2p_model.txt";

const std::string GAL_PCA_7P_M4_MODEL_TEXT_FILE = "/mount/sirocco1/imw2293/GROUP_CAT/SelfCalGroupFinder/py/parameters/bgs_y3/gal_pca_7p_m4_model.txt";
//const std::string GAL_PCA_7P_M4_Q_MODEL_TEXT_FILE = "/mount/sirocco1/imw2293/GROUP_CAT/SelfCalGroupFinder/py/parameters/bgs_y3/gal_pca_7p_m4_q_model.txt";
//const std::string GAL_PCA_7P_M4_SF_MODEL_TEXT_FILE = "/mount/sirocco1/imw2293/GROUP_CAT/SelfCalGroupFinder/py/parameters/bgs_y3/gal_pca_7p_m4_sf_model.txt";

struct Galaxy2P : Transformable {
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

struct Galaxy7P : Transformable {
    double abs_mag_r; // absolute magnitude in the r band, k-corrected to z=0.1, with h=1.0
    double logmstar;
    double color_g_r;   // g-r color (mag)
    double dn4000;
    double halflight_radius;
    double c9050;
    double loghalpha;
    double lat1, lat2, lat3, lat4, lat5, lat6, lat7;

    double getProperty(int num) override {
        switch (num) {
            case 0: return abs_mag_r;
            case 1: return logmstar;
            case 2: return color_g_r;
            case 3: return dn4000;
            case 4: return halflight_radius;
            case 5: return c9050;
            case 6: return loghalpha;
            default: throw std::out_of_range("Invalid property number");
        }
    }
    void setProperty(int num, double value) override {
        switch (num) {
            case 0: abs_mag_r = value; break;
            case 1: logmstar = value; break;
            case 2: color_g_r = value; break;
            case 3: dn4000 = value; break;
            case 4: halflight_radius = value; break;
            case 5: c9050 = value; break;
            case 6: loghalpha = value; break;
            default: throw std::out_of_range("Invalid property number");
        }
    }
    double getLatentProperty(int num) override {
        switch (num) {
            case 0: return lat1;
            case 1: return lat2;
            case 2: return lat3;
            case 3: return lat4;
            case 4: return lat5;
            case 5: return lat6;
            case 6: return lat7;
            default: throw std::out_of_range("Invalid latent property number");
        }
    }
    void setLatentProperty(int num, double value) override {
        switch (num) {
            case 0: lat1 = value; break;
            case 1: lat2 = value; break;
            case 2: lat3 = value; break;
            case 3: lat4 = value; break;
            case 4: lat5 = value; break;
            case 5: lat6 = value; break;
            case 6: lat7 = value; break;
            default: throw std::out_of_range("Invalid latent property number");
        }
    }
};

template<typename T>
struct Halo : Transformable {
    double x,y,z;
    double logmhalo; // log10 of M200b
    double halfmass_scale;
    double c;
    double spin;
    double lat1, lat2, lat3, lat4;

    // temp properties used for ranking
    double rank;

    T galaxy; // The galaxy that occupies this halo after abundance matching

    double getProperty(int num) override {
        switch (num) {
            case 0: return logmhalo;
            case 1: return halfmass_scale;
            case 2: return c;
            case 3: return spin;
            default: throw std::out_of_range("Invalid property number");
        }
    }
    void setProperty(int num, double value) override {
        switch (num) {
            case 0: logmhalo = value; break;
            case 1: halfmass_scale = value; break;
            case 2: c = value; break;
            case 3: spin = value; break;
            default: throw std::out_of_range("Invalid property number");
        }
    }
    double getLatentProperty(int num) override {
        switch (num) {
            case 0: return lat1;
            case 1: return lat2;
            case 2: return lat3;
            case 3: return lat4;
            default: throw std::out_of_range("Invalid latent property number");
        }
    }
    void setLatentProperty(int num, double value) override {
        switch (num) {
            case 0: lat1 = value; break;
            case 1: lat2 = value; break;
            case 2: lat3 = value; break;
            case 3: lat4 = value; break;
            default: throw std::out_of_range("Invalid latent property number");
        }
    }
};

template<typename T>
void sortByRank(std::vector<Halo<T>>& halos) {
    std::sort(halos.begin(), halos.end(),
        [](const Halo<T>& a, const Halo<T>& b) { 
            return a.rank > b.rank; 
        });
}



/** Hmm, this multi-model approach won't work because I need to already have assigned galaxy properties to figure out which model ot use! It's circular...
class QSFLatentModel : public ILatentModel {
public:

    QSFLatentModel(const std::string &quiescentModelFile, const std::string &starFormingModelFile) {
        this->quiescentModel = std::make_unique<LatentModel>(quiescentModelFile);
        this->starFormingModel = std::make_unique<LatentModel>(starFormingModelFile);

        // Ensure they have the same feature count
        if (this->quiescentModel->getFeatureCount() != this->starFormingModel->getFeatureCount()) {
            LOG_ERROR("Quiescent and star-forming galaxy models must have the same number of features\n");
            exit(1);
        }
    }

    void forward_transform(Transformable &obj) override {
        // TODO split based on if the galaxy is quiescent or star-forming and apply the corresponding model
        // Issue is that it is just a Transformable here, not a Galaxy7P
    }

    void inverse_transform(Transformable &obj) override {
        // TODO split based on if the galaxy is quiescent or star-forming and apply the corresponding model
    }

    void print_model() override {
        LOG_INFO("Quiescent Model:\n");
        quiescentModel->print_model();
        LOG_INFO("Star Forming Model:\n");
        starFormingModel->print_model();
    }

    ~QSFLatentModel() = default;


private:
    std::unique_ptr<LatentModel> quiescentModel;
    std::unique_ptr<LatentModel> starFormingModel;
};
*/


// Previously there was a separate singleton class per density function (GalaxyMagMatcher,
// GalaxyColorGMRMatcher, Galaxy2P_ICA1Matcher, Galaxy2P_ICA2Matcher) each with an identical
// match() implementation, differing only in which density function file they loaded. That
// common behavior now lives directly in AbundanceMatchingManager (see latsham.hpp), and
// AbundanceMatchingRegistry::get().getMatcher(filename) lazily creates/caches one instance
// per file. Call sites below just fetch the matcher for the density file they need.

template<typename T>
class ConnectionModel {
public:
    using GalaxyModel = T;
    virtual void setParamsFromList(std::vector<double> params) = 0;
    virtual bool match(std::vector<Halo<GalaxyModel>>& halos) = 0;
    virtual bool writeMocks(std::vector<Halo<GalaxyModel>>& halos) = 0;
    virtual void load() = 0;
    
    std::unique_ptr<ILatentModel> galaxyModel;
    std::unique_ptr<ILatentModel> haloModel;
};

class ConnectionModel_1 : public ConnectionModel<Galaxy2P> {
public:
    double scatter;

    void setParamsFromList(std::vector<double> params) {
        if (params.size() != 1) {
            LOG_ERROR("Expected 1 parameter for Model1, got %zu\n", params.size());
            throw std::runtime_error("Invalid number of parameters for Model1");
        }
        scatter = params[0];
    }

    bool match(std::vector<Halo<Galaxy2P>>& halos) override {

        // Lognormal scatter in luminosity at fixed halo mass
        std::random_device rd;
        std::mt19937 gen(rd());
        double mean = 0.0; 
        double std_dev = scatter;  
        std::normal_distribution<double> gaussian_dist(mean, std_dev);

        for (Halo<Galaxy2P> &h : halos) {
            h.rank = h.logmhalo + gaussian_dist(gen);
        }
        sortByRank<Galaxy2P>(halos);

        // Now abundance match to galaxy r band abs magnitude using
        double density = 0.0;
        for (Halo<Galaxy2P>& h : halos) {
            density += HALO_DENSITY_INCREMENT; // cumulative density of halos above this mass
            h.galaxy.abs_mag_r = AbundanceMatchingRegistry::get().getMatcher(GAL_MAGR_DENSITY).match(density);
        }

        return true;
    }    

    bool writeMocks(std::vector<Halo<Galaxy2P>>& halos) override {
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

class ConnectionModel_2P2P : public ConnectionModel<Galaxy2P> {
public:

    void load() {
        this->galaxyModel = std::make_unique<LatentModel>(GAL_2P_LATENT_MODEL_TEXT_FILE);
        this->haloModel = std::make_unique<LatentModel>(HALO_2P_LATENT_MODEL_TEXT_FILE);
        loaded = true;
    }

    double p1 = 1.0;
    double p2 = 0.0;
    double p3 = 0.0;
    double p4 = 1.0;

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

    bool match(std::vector<Halo<Galaxy2P>>& halos) override {
        if (!loaded) load();

        // Perf note: I tried storing a lookup of density values, but somehow it was slower than recomputing it.
        // No gains from parallelizing the match(density) call.
        // The sorts take ~250ms each on howdy
        // Not much to improve here. Replacing mock writing with direct call to corrfunc from code is most important.

        //#pragma omp parallel for
        for (size_t i = 0; i < halos.size(); ++i) {
            halos[i].rank = p1 * halos[i].lat1 + p2 * halos[i].lat2;
        }
        sortByRank<Galaxy2P>(halos);

        double density = 0.0;
        for (size_t i = 0; i < halos.size(); ++i) {
            density += HALO_DENSITY_INCREMENT; // cumulative density of halos above this mass
            halos[i].galaxy.setLatentProperty(0, AbundanceMatchingRegistry::get().getMatcher(GAL_2P_ICA1_DENSITY).match(density));
        }

        //#pragma omp parallel for
        for (size_t i = 0; i < halos.size(); ++i) {
            halos[i].rank = p3 * halos[i].lat1 + p4 * halos[i].lat2;
        }
        sortByRank<Galaxy2P>(halos);

        density = 0.0;
        for (size_t i = 0; i < halos.size(); ++i) {
            density += HALO_DENSITY_INCREMENT; // cumulative density of halos above this mass
            halos[i].galaxy.setLatentProperty(1, AbundanceMatchingRegistry::get().getMatcher(GAL_2P_ICA2_DENSITY).match(density));
        }

        #pragma omp parallel for
        for (size_t i = 0; i < halos.size(); ++i) {
            this->galaxyModel->inverse_transform(halos[i].galaxy);
        }

        // Re-abundance properties match onto itself to ensure original propety distributions are respected
        // This is a trick to deal with the higher-order residual correlations in the latent space,
        // which can result in unrealistic property values.

        //#pragma omp parallel for
        for (size_t i = 0; i < halos.size(); ++i) {
            halos[i].rank = halos[i].galaxy.getProperty(1); 
        }
        sortByRank<Galaxy2P>(halos);
        density = 0.0;
        for (Halo<Galaxy2P> &h : halos) {
            density += HALO_DENSITY_INCREMENT; 
            h.galaxy.setProperty(1, AbundanceMatchingRegistry::get().getMatcher(GAL_COLOR_GMR_DENSITY).match(density));
        }

        //#pragma omp parallel for
        for (size_t i = 0; i < halos.size(); ++i) {
            // Minus sign for magnitudes!!
            halos[i].rank = - halos[i].galaxy.getProperty(0); 
        }
        sortByRank<Galaxy2P>(halos);
        density = 0.0;
        for (Halo<Galaxy2P> &h : halos) {
            density += HALO_DENSITY_INCREMENT; 
            h.galaxy.setProperty(0, AbundanceMatchingRegistry::get().getMatcher(GAL_MAGR_DENSITY).match(density));
        }

        for (size_t i = 0; i < 5 && i < halos.size(); ++i) {
            const Halo<Galaxy2P>& h = halos[i];
            LOG_INFO("Halo %d: logmhalo=%.3f, age=%.3f, color=%.3f, abs_mag_r=%.3f\n", i, 
                h.logmhalo, h.halfmass_scale, h.galaxy.color_g_r, h.galaxy.abs_mag_r);
        }

        return true;
    }    

    bool writeMocks(std::vector<Halo<Galaxy2P>>& halos) override {
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

    bool writeMockForBin(const std::vector<Halo<Galaxy2P>>& halos, int magbin, const std::string& color, const std::string& out_filename) {
        std::ofstream outputFile(out_filename);
        double color_cut = 0.76; 
        int count = 0;

        // No satellites for this study; we will compare to centrals-only clustering measurements.    

        for (const Halo<Galaxy2P>& h : halos) {
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

class ConnectionModel_2P4P : public ConnectionModel<Galaxy2P> {
public:

    void load() {
        this->galaxyModel = std::make_unique<LatentModel>(GAL_2P_LATENT_MODEL_TEXT_FILE);
        this->haloModel = std::make_unique<LatentModel>(HALO_4P_LATENT_MODEL_TEXT_FILE);

        const std::string data_folder = "/mount/sirocco1/imw2293/GROUP_CAT/latsham/data/dr2_wp_1/";
        const std::regex pattern(R"(wp_mag([\d\.-]+to[\d\.-]+)_gmr([\d\.-]+to[\d\.-]+)\.dat)");
        std::set<double, std::greater<double>> magbins_set;
        std::map<std::string, std::set<double>> gmrbins_set;

        for (const auto &entry : std::filesystem::directory_iterator(data_folder)) {
            if (!entry.is_regular_file())
                continue;

            std::string fname = entry.path().filename().string();

            std::smatch match;
            if (!std::regex_match(fname, match, pattern))
                continue;

            std::string mag_range = match[1];
            std::string gr_range  = match[2];

            size_t pos1 = mag_range.find("to");
            double mag_min = std::stod(mag_range.substr(0, pos1));
            double mag_max = std::stod(mag_range.substr(pos1 + 2));
            magbins_set.insert(mag_min);
            magbins_set.insert(mag_max);

            size_t pos2 = gr_range.find("to");
            double gr_min = std::stod(gr_range.substr(0, pos2));
            double gr_max = std::stod(gr_range.substr(pos2 + 2));
            gmrbins_set[mag_range].insert(gr_min);
            gmrbins_set[mag_range].insert(gr_max);
        }

        if (magbins_set.empty()) {
            LOG_ERROR("No magnitude bins found in data folder %s\n", data_folder.c_str());
            throw std::runtime_error("No magnitude bins found");
        }

        magbins.assign(magbins_set.begin(), magbins_set.end());

        gmrbins.clear();

        for (auto it = gmrbins_set.begin(); it != gmrbins_set.end(); ++it)
        {
            gmrbins.emplace_back(it->second.begin(), it->second.end());
        }

        total_bins = 0;
        for (size_t i = 0; i < magbins.size() - 1; ++i) {
            LOG_VERBOSE("Mag bin %zu: [%.2f, %.2f]\n", i, magbins[i], magbins[i + 1]);
            if (i < gmrbins.size()) {
                LOG_VERBOSE("  GMR bins: ");
                for (double gmr : gmrbins[i]) {
                    LOG_VERBOSE("%.3f ", gmr);
                }
                LOG_VERBOSE("\n");
                total_bins += gmrbins[i].size() - 1; // number of bins is number of edges - 1
            }
        }

        loaded = true;
    }

    // Initial guess from loadings for LOGMHALO and AGE
    /*double g1_h1 = 0.999;
    double g1_h2 = -0.137;
    double g1_h3 = 0.26;
    double g1_h4 = -0.008;
    double g2_h1 = -0.111;
    double g2_h2 = 0.422;
    double g2_h3 = 1.183;
    double g2_h4 = -0.168;*/
    double g1_h1 = 0.8;
    double g1_h2 = 0.05;
    double g1_h3 = 0.15;
    double g1_h4 = 0.0;
    double g2_h1 = 0.0;
    double g2_h2 = 0.2;
    double g2_h3 = 0.8;
    double g2_h4 = 0.0;

    void setParamsFromList(std::vector<double> params) {
        if (params.size() != 8) {
            LOG_ERROR("Expected 8 parameters for Model_2P4P, got %zu\n", params.size());
            throw std::runtime_error("Invalid number of parameters for Model_2P4P");
        }
        g1_h1 = params[0];
        g1_h2 = params[1];
        g1_h3 = params[2];
        g1_h4 = params[3];
        g2_h1 = params[4];
        g2_h2 = params[5];
        g2_h3 = params[6];
        g2_h4 = params[7];
        LOG_INFO("Set model params: g1_h1=%.3f, g1_h2=%.3f, g1_h3=%.3f, g1_h4=%.3f, g2_h1=%.3f, g2_h2=%.3f, g2_h3=%.3f, g2_h4=%.3f\n", g1_h1, g1_h2, g1_h3, g1_h4, g2_h1, g2_h2, g2_h3, g2_h4);
    }

    bool match(std::vector<Halo<Galaxy2P>>& halos) override {
        if (!loaded) load();

        // Perf note: I tried storing a lookup of density values, but somehow it was slower than recomputing it.
        // No gains from parallelizing the match(density) call.
        // The sorts take ~250ms each on howdy
        // Not much to improve here. Replacing mock writing with direct call to corrfunc from code is most important.

        //#pragma omp parallel for
        for (size_t i = 0; i < halos.size(); ++i) {
            halos[i].rank = g1_h1 * halos[i].lat1 + g1_h2 * halos[i].lat2 + g1_h3 * halos[i].lat3 + g1_h4 * halos[i].lat4;
        }
        sortByRank<Galaxy2P>(halos);

        double density = 0.0;
        for (size_t i = 0; i < halos.size(); ++i) {
            density += HALO_DENSITY_INCREMENT; // cumulative density of halos above this mass
            halos[i].galaxy.setLatentProperty(0, AbundanceMatchingRegistry::get().getMatcher(GAL_2P_ICA1_DENSITY).match(density));
        }

        //#pragma omp parallel for
        for (size_t i = 0; i < halos.size(); ++i) {
            halos[i].rank = g2_h1 * halos[i].lat1 + g2_h2 * halos[i].lat2 + g2_h3 * halos[i].lat3 + g2_h4 * halos[i].lat4;
        }
        sortByRank<Galaxy2P>(halos);

        density = 0.0;
        for (size_t i = 0; i < halos.size(); ++i) {
            density += HALO_DENSITY_INCREMENT; // cumulative density of halos above this mass
            halos[i].galaxy.setLatentProperty(1, AbundanceMatchingRegistry::get().getMatcher(GAL_2P_ICA2_DENSITY).match(density));
        }

        // Get the original galaxy property values from the assigned latent ones
        #pragma omp parallel for
        for (size_t i = 0; i < halos.size(); ++i) {
            this->galaxyModel->inverse_transform(halos[i].galaxy);
        }

        // **************************************************************************************************
        // Re-abundance properties match onto itself to ensure original propety distributions are respected
        // This is a trick to deal with the higher-order residual correlations in the latent space,
        // which can result in unrealistic property values.
        // **************************************************************************************************

        //#pragma omp parallel for
        for (size_t i = 0; i < halos.size(); ++i) {
            halos[i].rank = halos[i].galaxy.getProperty(1); 
        }
        sortByRank<Galaxy2P>(halos);
        density = 0.0;
        for (Halo<Galaxy2P> &h : halos) {
            density += HALO_DENSITY_INCREMENT; 
            h.galaxy.setProperty(1, AbundanceMatchingRegistry::get().getMatcher(GAL_COLOR_GMR_DENSITY).match(density));
        }

        //#pragma omp parallel for
        for (size_t i = 0; i < halos.size(); ++i) {
            // Minus sign for magnitudes!!
            halos[i].rank = - halos[i].galaxy.getProperty(0); 
        }
        sortByRank<Galaxy2P>(halos);
        density = 0.0;
        for (Halo<Galaxy2P> &h : halos) {
            density += HALO_DENSITY_INCREMENT; 
            h.galaxy.setProperty(0, AbundanceMatchingRegistry::get().getMatcher(GAL_MAGR_DENSITY).match(density));
        }

        for (size_t i = 0; i < 5 && i < halos.size(); ++i) {
            const Halo<Galaxy2P>& h = halos[i];
            LOG_INFO("Halo %d: logmhalo=%.3f, age=%.3f, color=%.3f, abs_mag_r=%.3f\n", i, 
                h.logmhalo, h.halfmass_scale, h.galaxy.color_g_r, h.galaxy.abs_mag_r);
        }

        return true;
    }    

    bool writeMocks(std::vector<Halo<Galaxy2P>>& halos) override {
        // Write files out in the format that corrfun expects for wp calculation.
        int successes = 0;

        //#pragma omp parallel for collapse(2) reduction(+:successes)
        #pragma omp parallel for reduction(+:successes) // Only parallizes outer loop, but 10 at a time is fine...
        for (size_t i = 0; i < this->magbins.size() - 1; ++i) {
            for (size_t j = 0; j < this->gmrbins[i].size() - 1; ++j) {
                successes += this->writeMockForBin(halos, i, j);
            }
        }
        
        LOG_INFO("Successfully wrote %d mocks out of %zu\n", successes, this->total_bins);
        return successes == this->total_bins;
    }

private:
    bool loaded = false;
    std::vector<double> magbins;
    std::vector<std::vector<double>> gmrbins;
    size_t total_bins;

    bool writeMockForBin(const std::vector<Halo<Galaxy2P>>& halos, size_t i, size_t j) {
        //double color_cut = 0.76; 
        int count = 0;

        auto mag_left = magbins[i];
        auto mag_right = magbins[i + 1];
        auto gmr_left = gmrbins[i][j];
        auto gmr_right = gmrbins[i][j + 1];
        std::ostringstream oss;
        oss << "/mount/sirocco1/imw2293/GROUP_CAT/OUTPUT/LATSHAM/mock_"
        << std::fixed << std::setprecision(4) << mag_left << "_" << mag_right
        << "_" << gmr_left << "_" << gmr_right << ".dat";
        std::ofstream outputFile(oss.str());

        // No satellites for this study; we will compare to centrals-only clustering measurements
        for (const Halo<Galaxy2P>& h : halos) {
            if (h.galaxy.abs_mag_r > mag_left || h.galaxy.abs_mag_r < mag_right) // They go like -16, -17, -18...
                continue;
            if (h.galaxy.color_g_r < gmr_left || h.galaxy.color_g_r > gmr_right)
                continue;

            ++count;
            outputFile << std::setprecision(8) << h.x << " " << h.y << " " << h.z << "\n";
        }
        LOG_VERBOSE("Mock building: wrote %d galaxies for magbin [%.2f, %.2f], gmrbin [%.3f, %.3f]\n", count, mag_left, mag_right, gmr_left, gmr_right);
        return count > 0;
    }
};


struct BiasInfo {
    double magbin_brighter;
    double magbin_fainter;
    double magbin_mean;
    bool quiescent;
    std::string third_property_name;
    double third_property_lower;
    double third_property_upper;
    double third_property_mean;
    double b; // overall bias
    double b_err_low;
    double b_err_high;
    double b_p; // property bias
    double b_p_err_low;
    double b_p_err_high;
    double b_m; // magnitude bias
};

/**
 * This connection model uses:
 * 1) The DR2 v0.3 bias measurements
 * 2) The galaxy 7-prop PCA model from the 4th mag bin (~ -21 mag)
 * 3) The generic halo 4-prop ICA model
 */
class ConnectionModel_DR2v0_3_7P4P : public ConnectionModel<Galaxy7P> {
public:

    void load() {
        this->galaxyModel = std::make_unique<LatentModel>(GAL_PCA_7P_M4_MODEL_TEXT_FILE); // right now this is derived only from ~ -21 mag galaxies as an experiment
        this->haloModel = std::make_unique<LatentModel>(HALO_4P_LATENT_MODEL_TEXT_FILE);
        this->galaxyModel->load();
        this->haloModel->load();

        // Ensure the number of properties match
        if (this->galaxyModel->getFeatureCount() != 7) {
            LOG_ERROR("Expected 7 features for galaxy model, got %d\n", this->galaxyModel->getFeatureCount());
            throw std::runtime_error("Invalid number of features for galaxy model");
        }
        if (this->haloModel->getFeatureCount() != 4) {
            LOG_ERROR("Expected 4 features for halo model, got %d\n", this->haloModel->getFeatureCount());
            throw std::runtime_error("Invalid number of features for halo model");
        }

        const std::string data_folder = "/mount/sirocco1/imw2293/GROUP_CAT/latsham/data/dr2_wp_2/";
        const std::string path = data_folder + "biases_BGS_DR2_5-30Mpc_v0.3.txt";

        std::ifstream inputFile(path);
        if (!inputFile.is_open()) {
            LOG_ERROR("Failed to open bias data file: %s\n", path.c_str());
            throw std::runtime_error("Failed to open bias data file");
        }

        // Skip the header line (column names); without this the first >> would try
        // to parse "magbin_fainter" as a double, fail immediately, and the while loop
        // below would never execute, silently loading 0 bias measurements.
        std::string header_line;
        std::getline(inputFile, header_line);

        BiasInfo bias;
        double quiescent_flag;
        // Note: the file's column order is magbin_fainter then magbin_brighter.
        // quiescent is stored in the file as "1.0"/"0.0", not "1"/"0", so it must be
        // read into a double first (bool>>'s integer-only parser would read just "1"
        // from "1.0", leave ".0" in the stream, and desync every field after it).
        while (inputFile >> bias.magbin_fainter >> bias.magbin_brighter >> bias.magbin_mean
                         >> quiescent_flag >> bias.third_property_name
                         >> bias.third_property_lower >> bias.third_property_upper >> bias.third_property_mean
                         >> bias.b >> bias.b_err_low >> bias.b_err_high
                         >> bias.b_p >> bias.b_p_err_low >> bias.b_p_err_high
                         >> bias.b_m) {
            bias.quiescent = quiescent_flag != 0.0;

            // Only keep the right magnitude bin for now while we focuson on the ~ -21 mag bin
            if (bias.magbin_brighter >= -21.1 && bias.magbin_fainter <= -20.70)
                measured_biases.push_back(bias);
        }
        inputFile.close();

        // Print off how many bias measurements we kept
        LOG_INFO("Loaded %zu bias measurements for magbin ~ -21 from %s\n", measured_biases.size(), path.c_str());

        loaded = true;
    }

    double parameters[7][4] = {
        {0.0, 0.5, -0.5, 0.0}, // red and old, diffuse
        {1.0, 0.0, 0.0, 0.0},   // bright, massive
        {-1.0, 0.0, 0.0, 0.0},   // dim and large
        {-1.0, 0.0, 0.0, 0.0},   // compact, massive, dim
        {0.0, 0.0, -1.0, 0.0},   // diffuse, young, starforming
        {0.0, 0.0, -1.0, 0.0}, // old with halpha 
        {0.0, 0.0, 0.0, 1.0} // low mass edge on dusty star-forming
    };

    void setParamsFromList(std::vector<double> params) {
        if (params.size() != 7*4) {
            LOG_ERROR("Expected 28 parameters for ConnectionModel_DR2v0_3_7P4P, got %zu\n", params.size());
            throw std::runtime_error("Invalid number of parameters for ConnectionModel_DR2v0_3_7P4P");
        }
        for (size_t i = 0; i < 7; ++i) {
            for (size_t j = 0; j < 4; ++j) {
                parameters[i][j] = params[i*4 + j];
            }
        }
        LOG_INFO("Set model params for ConnectionModel_DR2v0_3_7P4P\n");
    }

    bool match(std::vector<Halo<Galaxy7P>>& halos) override {
        if (!loaded) load();

        // Perf note: I tried storing a lookup of density values, but somehow it was slower than recomputing it.
        // No gains from parallelizing the match(density) call.
        // The sorts take ~250ms each on howdy
        // Not much to improve here. Replacing mock writing with direct call to corrfunc from code is most important.

        // SET GALAXY PCs
        for (size_t propNum = 0; propNum < 7; ++propNum) {

            //#pragma omp parallel for
            for (size_t i = 0; i < halos.size(); ++i) {
                halos[i].rank = this->parameters[propNum][0] * halos[i].lat1 + this->parameters[propNum][1] * halos[i].lat2 + this->parameters[propNum][2] * halos[i].lat3 + this->parameters[propNum][3] * halos[i].lat4;
            }
            sortByRank<Galaxy7P>(halos);

            // TODO think about if HALO_DENSITY_INCREMENT is still right for the downsampled model...
            double density = 0.0;
            auto matcher =  AbundanceMatchingRegistry::get().getMatcher(galDensityFilename("7p", propNum));
            for (size_t i = 0; i < halos.size(); ++i) {
                density += HALO_DENSITY_INCREMENT; // cumulative density of halos above this mass
                halos[i].galaxy.setLatentProperty(propNum, matcher.match(density));
            }
        }

        // Get the original galaxy property values from the assigned latent ones
        #pragma omp parallel for
        for (size_t i = 0; i < halos.size(); ++i) {
            this->galaxyModel->inverse_transform(halos[i].galaxy);
        }

        // **************************************************************************************************
        // Re-abundance properties match onto itself to ensure original propety distributions are respected
        // This is a trick to deal with the higher-order residual correlations in the latent space,
        // which can result in unrealistic property values.
        // **************************************************************************************************

        for (size_t i = 0; i < halos.size(); ++i) {
            halos[i].rank = halos[i].galaxy.color_g_r;
        }
        sortByRank(halos);
        double density = 0.0;
        for (Halo<Galaxy7P> &h : halos) {
            density += HALO_DENSITY_INCREMENT; 
            h.galaxy.color_g_r = AbundanceMatchingRegistry::get().getMatcher(GAL_MAGBIN4_COLOR_GMR_DENSITY).match(density);
        }

        for (size_t i = 0; i < halos.size(); ++i) {
            // Minus sign for magnitudes!!
            halos[i].rank = - halos[i].galaxy.abs_mag_r;
        }
        sortByRank(halos);
        density = 0.0;
        for (Halo<Galaxy7P> &h : halos) {
            density += HALO_DENSITY_INCREMENT; 
            h.galaxy.setProperty(0, AbundanceMatchingRegistry::get().getMatcher(GAL_MAGBIN4_MAGR_DENSITY).match(density));
        }

        for (size_t i = 0; i < halos.size(); ++i) {
            halos[i].rank = halos[i].galaxy.c9050;
        }
        sortByRank(halos);
        density = 0.0;
        for (Halo<Galaxy7P> &h : halos) {
            density += HALO_DENSITY_INCREMENT; 
            h.galaxy.c9050 = AbundanceMatchingRegistry::get().getMatcher(GAL_MAGBIN4_C9050_DENSITY).match(density);
        }

        for (size_t i = 0; i < halos.size(); ++i) {
            halos[i].rank = halos[i].galaxy.dn4000;
        }
        sortByRank(halos);
        density = 0.0;
        for (Halo<Galaxy7P> &h : halos) {
            density += HALO_DENSITY_INCREMENT; 
            h.galaxy.dn4000 = AbundanceMatchingRegistry::get().getMatcher(GAL_MAGBIN4_DN4000_DENSITY).match(density);
        }

        for (size_t i = 0; i < halos.size(); ++i) {
            halos[i].rank = halos[i].galaxy.halflight_radius;
        }
        sortByRank(halos);
        density = 0.0;
        for (Halo<Galaxy7P> &h : halos) {
            density += HALO_DENSITY_INCREMENT; 
            h.galaxy.halflight_radius = AbundanceMatchingRegistry::get().getMatcher(GAL_MAGBIN4_HALFLIGHT_DENSITY).match(density);
        } // TODO why is this negative

        for (size_t i = 0; i < halos.size(); ++i) {
            halos[i].rank = halos[i].galaxy.loghalpha;
        }
        sortByRank(halos);
        density = 0.0;
        for (Halo<Galaxy7P> &h : halos) {
            density += HALO_DENSITY_INCREMENT; 
            h.galaxy.loghalpha = AbundanceMatchingRegistry::get().getMatcher(GAL_MAGBIN4_LOGHALPHA_DENSITY).match(density);
        } // TODO why is this negative

        for (size_t i = 0; i < halos.size(); ++i) {
            halos[i].rank = halos[i].galaxy.logmstar;
        }
        sortByRank(halos);
        density = 0.0;
        for (Halo<Galaxy7P> &h : halos) {
            density += HALO_DENSITY_INCREMENT; 
            h.galaxy.logmstar = AbundanceMatchingRegistry::get().getMatcher(GAL_MAGBIN4_LOGMSTAR_DENSITY).match(density);
        }
        
        for (size_t i = 0; i < 5 && i < halos.size(); ++i) {
            const Halo<Galaxy7P>& h = halos[i];
            LOG_INFO("Halo %d: logmhalo=%.3f, a1/2=%.3f, color=%.3f, abs_mag_r=%.3f, c9050=%.3f, dn4000=%.3f, halflight_radius=%.3f, loghalpha=%.3f, logmstar=%.3f\n", i, 
                h.logmhalo, h.halfmass_scale, h.galaxy.color_g_r, h.galaxy.abs_mag_r,
                h.galaxy.c9050, h.galaxy.dn4000, h.galaxy.halflight_radius,
                h.galaxy.loghalpha, h.galaxy.logmstar);
        }

        return true;
    }    

    bool writeMocks(std::vector<Halo<Galaxy7P>>& halos) override {
        // Write files out in the format that corrfun expects for wp calculation.
        int successes = 0;
        size_t total_bins = this->measured_biases.size();

        if (total_bins <= 0) {
            LOG_ERROR("No bias measurements loaded, cannot write mocks\n");
            return false;
        }

        #pragma omp parallel for reduction(+:successes)
        for (size_t i = 0; i < this->measured_biases.size() - 1; ++i) {
            successes += this->writeMockForBin(halos, i);
        }
        
        LOG_INFO("Successfully wrote %d mocks out of %zu\n", successes, total_bins);
        return successes == total_bins; 
    }

private:
    bool loaded = false;
    std::vector<BiasInfo> measured_biases;

    bool writeMockForBin(std::vector<Halo<Galaxy7P>>& halos, size_t bin) {
        int count = 0;
        BiasInfo &info = this->measured_biases[bin];

        std::ostringstream oss;
        oss << "/mount/sirocco1/imw2293/GROUP_CAT/OUTPUT/LATSHAM/mock_" << bin;
        std::ofstream outputFile(oss.str());

        // Write the halos whose galaxy (central) fits the criteria defined by this bias bin.
        for (const Halo<Galaxy7P> &h : halos) {
            if (h.galaxy.abs_mag_r < info.magbin_brighter || h.galaxy.abs_mag_r > info.magbin_fainter)
                continue;
            // TODO need updated galaxy model for the other properties

            ++count;
            outputFile << std::setprecision(8) << h.x << " " << h.y << " " << h.z << "\n";
        }

        LOG_VERBOSE("Mock building: wrote %d galaxies for bias bin %zu\n", count, bin);
        return count > 0;
    }
};