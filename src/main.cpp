/***********************************************************************
 * This file is part of the latsham package. 
 * Copyright (c) 2026 Ian Williams under the MIT License.
 ***********************************************************************/

#include <cstdio>
#include <unistd.h>
#include <argp.h>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <chrono>
#include <algorithm>
#include <string>
#include <map>
#include <stdexcept>
#include <highfive/H5File.hpp>
#include <hdf5.h>
#include <random>
#include "latsham.hpp"
#include "models.hpp"

// Message passing protocol via a pipe to python wrapper
#define MSG_REQUEST 0
#define MSG_FSAT 1
#define MSG_LHMR 2
#define MSG_LSAT 3
#define MSG_HOD 4
#define MSG_HODFIT 5
#define MSG_COMPLETED 6
#define MSG_ABORTED 7
#define TYPE_FLOAT 0
#define TYPE_DOUBLE 1

static FILE* MSG_PIPE = nullptr;

// This file is already cut to central halos only; no subhalos. That's the 'C' in the name.
// It is also cut to M200b > 8E9 Msun
const std::string HALO_CATALOG = "/mount/sirocco1/imw2293/GROUP_CAT/DATA/POPMOCK/smdpl_z0.19717.M8E9.C.h5";
double constexpr BOX_SIZE = 400.0; // Mpc/h, from the simulation specs
double constexpr SIM_VOLUME = BOX_SIZE * BOX_SIZE * BOX_SIZE; // (Mpc/h)^3
double constexpr HALO_DENSITY_INCREMENT = 1.0 / SIM_VOLUME; 
// ('index', 'ID', 'upid', 'M200b', 'Mpeak', 'mvir', 'rvir', 'rs', 'vmax', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'Halfmass_Scale', 'Spin', 'c', 'LOGMHALO', 'NOISE')

static std::string to_lower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), ::tolower);
    return s;
}

class Timer {
public:
    void begin() {
        start = std::chrono::high_resolution_clock::now();
    }
    void endAndLog(const std::string& message) {
        auto now = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(now - start).count();
        LOG_INFO("%s: %lld ms\n", message.c_str(), duration);
    }
private:
    std::chrono::high_resolution_clock::time_point start;
};


// Reads one field from a compound dataset by case-insensitive name match.
// HDF5 will convert float->double automatically.
static std::vector<double> read_field(hid_t ds, const std::string& field_lower,
                                      const std::map<std::string, std::string>& name_map,
                                      hsize_t n) {
    auto it = name_map.find(field_lower);
    if (it == name_map.end())
        throw std::runtime_error("Column not found (case-insensitive): " + field_lower);

    hid_t mem_type = H5Tcreate(H5T_COMPOUND, sizeof(double));
    H5Tinsert(mem_type, it->second.c_str(), 0, H5T_NATIVE_DOUBLE);

    std::vector<double> buf(n);
    herr_t err = H5Dread(ds, mem_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf.data());
    H5Tclose(mem_type);

    if (err < 0)
        throw std::runtime_error("H5Dread failed for column: " + it->second);
    return buf;
}

std::vector<Halo> read_halo_catalog(const std::string& path) {
    HighFive::File file(path, HighFive::File::ReadOnly);
    auto dataset = file.getDataSet("halos/table");
    hsize_t n = dataset.getDimensions()[0];

    hid_t ds_id = dataset.getId();
    hid_t file_type = H5Dget_type(ds_id);
    int nmembers = H5Tget_nmembers(file_type);
    std::map<std::string, std::string> name_map;
    for (int i = 0; i < nmembers; ++i) {
        char* raw = H5Tget_member_name(file_type, i);
        std::string name(raw);
        H5free_memory(raw);
        name_map[to_lower(name)] = name;
    }
    H5Tclose(file_type);

    auto x_v        = read_field(ds_id, "x",              name_map, n);
    auto y_v        = read_field(ds_id, "y",              name_map, n);
    auto z_v        = read_field(ds_id, "z",              name_map, n);
    //auto vx_v       = read_field(ds_id, "vx",             name_map, n);
    //auto vy_v       = read_field(ds_id, "vy",             name_map, n);
    //auto vz_v       = read_field(ds_id, "vz",             name_map, n);
    auto logm_v     = read_field(ds_id, "logmhalo",       name_map, n);
    auto halfm_v    = read_field(ds_id, "halfmass_scale", name_map, n);
    auto c_v        = read_field(ds_id, "c",              name_map, n);
    auto spin_v     = read_field(ds_id, "spin",           name_map, n);

    std::vector<Halo> halos(n);
    for (hsize_t i = 0; i < n; ++i) {
        halos[i].x              = x_v[i];
        halos[i].y              = y_v[i];
        halos[i].z              = z_v[i];
        //halos[i].vx             = vx_v[i];
        //halos[i].vy             = vy_v[i];
        //halos[i].vz             = vz_v[i];
        halos[i].logmhalo       = logm_v[i];
        halos[i].halfmass_scale = halfm_v[i];
        halos[i].c              = c_v[i];
        halos[i].spin           = spin_v[i];
    }

    LOG_INFO("Read %zu halos from catalog %s\n", halos.size(), path.c_str());
    return halos;
}

static char doc[] = "latsham: Latent model abundance matching of a halos catalog"; 

static char args_doc[] = "";

static struct argp_option options[] = {
  {"pipe",         'P', "PIPEID",                             0,  "Specify a pipe ID communication with emcee manager", 1},
  { 0 }
};

/* Parse a single option. */
static error_t parse_opt (int key, char *arg, struct argp_state *state) {
    int pipe_id;
    switch (key) {
        case 'P':
            if (arg == NULL) {
                LOG_ERROR("Pipe ID must be specified with -P option.\n");
                exit(EPERM);
            }
            pipe_id = atoi(arg);
            if (pipe_id < 0) {
                LOG_ERROR("Invalid pipe ID: %d\n", pipe_id);
                exit(EPERM);
            }
            MSG_PIPE = fdopen(pipe_id, "w");
            if (!MSG_PIPE) 
                perror("fdopen");
            break;   
        default:
            return ARGP_ERR_UNKNOWN;
    }
    return 0;
}
static struct argp argp = { options, parse_opt, args_doc, doc };

void send_completed_msg() {
    if (MSG_PIPE != NULL) {
      //LOG_INFO("Group finding completed, sending message and awaiting next request...\n");
      uint8_t resp_msg_type = MSG_COMPLETED;
      uint8_t resp_data_type = TYPE_FLOAT;
      uint32_t resp_count = 0;

      fwrite(&resp_msg_type, 1, 1, MSG_PIPE);
      fwrite(&resp_data_type, 1, 1, MSG_PIPE);
      fwrite(&resp_count, sizeof(uint32_t), 1, MSG_PIPE);
      fflush(MSG_PIPE);
    }
}

void sortByRank(std::vector<Halo>& halos) {
    std::sort(halos.begin(), halos.end(),
        [](const Halo& a, const Halo& b) { 
            return a.rank > b.rank; 
        });
}

class ConnectionModel {
public:
    virtual void setParamsFromList(std::vector<double> params) = 0;
    virtual void match(std::vector<Halo>& halos) = 0;
    virtual void writeMocks(std::vector<Halo>& halos) = 0;
};

class ConnectionModel_1 : ConnectionModel {
public:
    double scatter;

    void setParamsFromList(std::vector<double> params) {
        if (params.size() != 1) {
            LOG_ERROR("Expected 1 parameter for Model1, got %zu\n", params.size());
            throw std::runtime_error("Invalid number of parameters for Model1");
        }
        scatter = params[0];
    }

    void match(std::vector<Halo>& halos) override {

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
    }    

    void writeMocks(std::vector<Halo>& halos) override {
        // Write files out in the format that corrfun expects for wp calculation.
        std::string out_filename = "/mount/sirocco1/imw2293/GROUP_CAT/OUTPUT/LATSHAM/mock_M20.dat";
        std::ofstream outputFile(out_filename);

        // Write the central galaxies for this mock
        for (int i=0; i < halos.size(); ++i) {
            if (halos[i].galaxy.abs_mag_r < -21.0 || halos[i].galaxy.abs_mag_r > -20.0)
                continue;

            outputFile << halos[i].x << " " << halos[i].y << " " << halos[i].z << "\n";
        }
        // No satellites for this study; we will compare to centrals-only clustering measurements.    
    }
};

class ConnectionModel_2P2P : ConnectionModel {
public:

    void load() {
        this->galaxyModel = std::make_unique<LatentModel>(GAL_2P_LATENT_MODEL_TEXT_FILE);
        this->haloModel = std::make_unique<LatentModel>(HALO_2P_LATENT_MODEL_TEXT_FILE);
    }

    std::unique_ptr<LatentModel> galaxyModel;
    std::unique_ptr<LatentModel> haloModel;
    double p1 = 1.0;
    double p2 = 0.0;
    double p3 = 0.0; 
    double p4 = 1.0;
    
    void setParamsFromList(std::vector<double> params) {
        if (params.size() != 4) {
            LOG_ERROR("Expected 4 parameters for Model_2L, got %zu\n", params.size());
            throw std::runtime_error("Invalid number of parameters for Model_2L");
        }
        p1 = params[0];
        p2 = params[1];
        p3 = params[2];
        p4 = params[3];
        LOG_INFO("Set model params: p1=%.3f, p2=%.3f, p3=%.3f, p4=%.3f\n", p1, p2, p3, p4);
    }

    void match(std::vector<Halo>& halos) override {

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

        // TODO re-abundance match onto itself (write tests to force this first)
        /* 
        for (Halo &h : halos) {
            h.rank = h.galaxy.getLatentProperty(0); // re-apply the same latent property values to force them to be self-consistent with the mag/color distribution after the inverse transform. This is a hack to deal with the higher-order residual correlation issue.
        }
        sortByRank(halos);

        density = 0.0;
        for (Halo &h : halos) {
            density += HALO_DENSITY_INCREMENT; 
            h.galaxy.setLatentProperty(0, GalaxyMagMatcher::get().match(density));
        }
        */
        
    }    

    void writeMocks(std::vector<Halo>& halos) override {
        // Write files out in the format that corrfun expects for wp calculation.
        std::string out_basename = "/mount/sirocco1/imw2293/GROUP_CAT/OUTPUT/LATSHAM/mock_%s_M%d.dat";

        std::vector<int> magbins = {-17, -18, -19, -20, -21, -22};
        std::vector<std::string> colors = {"red", "blue"};

        #pragma omp parallel for collapse(2)
        for (int i = 0; i < magbins.size(); ++i) {
            for (int j = 0; j < colors.size(); ++j) {
                auto magbin = magbins[i];
                auto color = colors[j];
                std::string out_filename = out_basename;
                out_filename = out_filename.replace(out_filename.find("%s"), 2, color);
                out_filename = out_filename.replace(out_filename.find("%d"), 2, std::to_string(abs(magbin)));
                LOG_VERBOSE("Writing mock for magbin %d, color %s to %s\n", magbin, color.c_str(), out_filename.c_str());
                writeMockForBin(halos, magbin, color, out_filename);
            }
        }

    }

private:
    void writeMockForBin(const std::vector<Halo>& halos, int magbin, const std::string& color, const std::string& out_filename) {
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
    }
};


std::vector<double> await_request() {
  uint8_t msg_type, data_type;
  uint32_t count;
  std::vector<double> params;

  // Read header: 1 byte msg_type, 1 byte data_type, 4 bytes count (little-endian)
  // TODO end gracefully when stdin is closed
  //if (feof(stdin)) {
  //  LOG_INFO("End of input stream reached, exiting...\n");
  //  return false; // No more requests
  //}
  uint8_t header[6];
  size_t n = fread(header, 1, 6, stdin);
  if (n != 6) {
    LOG_ERROR("Failed to read MSG_REQUEST header (got %zu bytes)\n", n);
    throw std::runtime_error("Failed to read MSG_REQUEST header");
  }
  msg_type = header[0];
  data_type = header[1];
  memcpy(&count, header + 2, sizeof(uint32_t));
  if (msg_type != MSG_REQUEST) {
    LOG_ERROR("Expected MSG_REQUEST, got %d\n", msg_type);
    throw std::runtime_error("Expected MSG_REQUEST");
  }
  if (data_type != TYPE_DOUBLE) {
    LOG_ERROR("Expected TYPE_DOUBLE, got %d\n", data_type);
    throw std::runtime_error("Expected TYPE_DOUBLE");
  }

  // Read payload: count doubles
  params.resize(count);
  size_t n_payload = fread(params.data(), sizeof(double), count, stdin);
  if (n_payload != count) {
    fprintf(stderr, "Failed to read MSG_REQUEST payload (got %zu doubles)\n", n_payload);
    throw std::runtime_error("Failed to read MSG_REQUEST payload");
  }

  return params;
}


int main(int argc, char **argv) {

    LOG_INFO("Welcome to Latent-SHAM\n");

    argp_parse (&argp, argc, argv, 0, 0, nullptr);

    srand(5981);

    Timer timer;
    timer.begin();
    std::vector<Halo> halos = read_halo_catalog(HALO_CATALOG);
    timer.endAndLog("Time to read halo catalog");

    // Remove halos below a mass threshold for now
    double logmhalo_cut = 11.0;
    halos.erase(std::remove_if(halos.begin(), halos.end(),
        [logmhalo_cut](const Halo& h) { return h.logmhalo < logmhalo_cut; }), halos.end());

    //ConnectionModel_1 model;
    //model.scatter = 0.83; 

    ConnectionModel_2P2P model = ConnectionModel_2P2P();
    model.load();

    // Apply halo model to get latent parameters
    int k = 0;
    for (Halo &h : halos) {
        model.haloModel->forward_transform(h);
    }

    timer.begin();
    model.match(halos);
    timer.endAndLog("Time for initial abundance matching");
    timer.begin();
    model.writeMocks(halos);
    timer.endAndLog("Time to write initial mocks");

    std::vector<double> params;
    timer.begin();
    while (MSG_PIPE != NULL) {

        try {
            params = await_request();
            model.setParamsFromList(params);
            timer.endAndLog("Time to get new parameters from pipe");
        } catch (const std::exception& e) {
            LOG_ERROR("Error while awaiting request: %s\n", e.what());
            exit(1);
        }

        timer.begin();
        model.match(halos);
        timer.endAndLog("Time to match with new parameters");

        timer.begin();
        model.writeMocks(halos);
        timer.endAndLog("Time to write mocks with new parameters");

        timer.begin();
        send_completed_msg();
    }

    return 0;
}