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
// ('index', 'ID', 'upid', 'M200b', 'Mpeak', 'mvir', 'rvir', 'rs', 'vmax', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'Halfmass_Scale', 'Spin', 'c', 'LOGMHALO', 'NOISE')

struct halo {
    // Halo Properties
    double x,y,z;
    double vx,vy,vz;
    double logmhalo; // log10 of M200b
    double halfmass_scale;
    double c;
    double spin;

    // temp property used for ranking
    double rank;

    // Assigned Galaxy Properties
    double abs_mag_r; // absolute magnitude in the r band, k-corrected to z=0.1, with h=1.0
    double color_g_r;   // g-r color (mag)
};

static std::string to_lower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), ::tolower);
    return s;
}

void run_corrfunc() {
    //  cmd = f"{corrfunc_path}/wp {boxsize} mock_{color}_M{m}.dat a {self.caldata.rpbinsfile} {pimax} {nthreads} > wp_mock_{color}_M{m}.dat 2> wp_stderr.txt"
    const std::string path = "/mount/sirocco1/tinker/src/Corrfunc/bin/wp";

}


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

std::vector<halo> read_halo_catalog(const std::string& path) {
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
    auto vx_v       = read_field(ds_id, "vx",             name_map, n);
    auto vy_v       = read_field(ds_id, "vy",             name_map, n);
    auto vz_v       = read_field(ds_id, "vz",             name_map, n);
    auto logm_v     = read_field(ds_id, "logmhalo",       name_map, n);
    auto halfm_v    = read_field(ds_id, "halfmass_scale", name_map, n);
    auto c_v        = read_field(ds_id, "c",              name_map, n);
    auto spin_v     = read_field(ds_id, "spin",           name_map, n);

    std::vector<halo> halos(n);
    for (hsize_t i = 0; i < n; ++i) {
        halos[i].x              = x_v[i];
        halos[i].y              = y_v[i];
        halos[i].z              = z_v[i];
        halos[i].vx             = vx_v[i];
        halos[i].vy             = vy_v[i];
        halos[i].vz             = vz_v[i];
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


struct Model1Params {
    double scatter;
};


void do_matching(Model1Params params, std::vector<halo>& halos) {

    auto t0 = std::chrono::high_resolution_clock::now();

    // Lognormal scatter in luminosity at fixed halo mass
    std::random_device rd;
    std::mt19937 gen(rd());
    double mean = 0.0; 
    double std_dev = params.scatter;  
    std::normal_distribution<double> gaussian_dist(mean, std_dev);

    for (halo &h : halos) {
        h.rank = h.logmhalo + gaussian_dist(gen);
    }

    std::sort(halos.begin(), halos.end(),
        [](const halo& a, const halo& b) { 
            return a.rank > b.rank; 
        });
    
    auto t1 = std::chrono::high_resolution_clock::now();

    // Now abundance match to galaxy r band abs magnitude using
    double density = 0.0;
    for (halo& h : halos) {
        density += 1.0 / SIM_VOLUME; // cumulative density of halos above this mass
        h.abs_mag_r = GalaxyMagMatcher::get().match(density);
    }

    auto t2 = std::chrono::high_resolution_clock::now();

    // Write files out in the format that corrfun expects for wp calculation.
    std::string out_filename = "/mount/sirocco1/imw2293/GROUP_CAT/OUTPUT/LATSHAM/mock_M20.dat";
    std::ofstream outputFile(out_filename);

    // Loop through halos and populate with galaxies
    // TODO different lum/color bin mocks
    // For now lets do -20 to -21 mag only

    // Write the central galaxies for this mock
    for (int i=0; i < halos.size(); ++i) {
        if (halos[i].abs_mag_r < -21.0 || halos[i].abs_mag_r > -20.0)
            continue;

        // Corrfunc only needs x,y,z locations. The rest was stuff group finder wrote; not sure why.
        outputFile << halos[i].x << " " << halos[i].y << " " << halos[i].z << "\n";  //<< " "
                   //<< halos[i].vx << " " << halos[i].vy << " " << halos[i].vz << " "
                   //<< 0 << " " 
                   //<< halos[i].logmhalo << "\n";
    }
    // No satellites for this study; we will compare to centrals-only clustering measurements.

    auto t3 = std::chrono::high_resolution_clock::now();

    // Done. Send pipe message that's were done and await instructions
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

    auto duration1 = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
    auto duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    auto duration3 = std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2).count();

    LOG_INFO("do_matching() completed for scatter = %.4f. Perf: sort = %ld ms, match = %ld ms, write = %ld ms\n",
             params.scatter, duration1, duration2, duration3);
}

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

    LOG_INFO("Welcome to Latent-SHAM");

    argp_parse (&argp, argc, argv, 0, 0, nullptr);

    srand(5981);

    auto t0 = std::chrono::high_resolution_clock::now();

    std::vector<halo> halos = read_halo_catalog(HALO_CATALOG);

    auto t1 = std::chrono::high_resolution_clock::now();

    // Remove halos below a mass threshold for now
    double logmhalo_cut = 11.0;
    halos.erase(std::remove_if(halos.begin(), halos.end(),
        [logmhalo_cut](const halo& h) { return h.logmhalo < logmhalo_cut; }), halos.end());

    auto t2 = std::chrono::high_resolution_clock::now();

    Model1Params model_params;
    model_params.scatter = 0.83; 
    do_matching(model_params, halos);
    
    bool run = MSG_PIPE != NULL; 
    while (run) {
        std::vector<double> params;

        // Read parameters from pipe
        try {
            params = await_request();
            model_params.scatter = params[0];
        } catch (const std::exception& e) {
            LOG_ERROR("Error while awaiting request: %s\n", e.what());
            break;
        }

        do_matching(model_params, halos);
    }

    // Print perf summary
    auto duration_read = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
    auto duration_filter = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    std::cerr << "Performance summary:\n";
    std::cerr << "  Read halo catalog: " << duration_read << " ms\n";
    std::cerr << "  Filter halos: " << duration_filter << " ms\n";


    return 0;
}