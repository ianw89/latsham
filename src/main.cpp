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
#include "latsham.hpp"


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

    std::cout << "Read " << n << " halos." << std::endl;
    return halos;
}

int main() {
    std::cout << "Welcome to Latent-SHAM" << std::endl;
    srand(5981);

    std::vector<halo> halos = read_halo_catalog(HALO_CATALOG);

    // Sort by a function of the halo properties
    std::sort(halos.begin(), halos.end(),
        [](const halo& a, const halo& b) { 
            return a.logmhalo > b.logmhalo; 
        });



    // Now abundance match to galaxy r band abs magnitude using
    double density = 0.0;
    for (halo& h : halos) {
        density += 1.0 / SIM_VOLUME; // cumulative density of halos above this mass
        h.abs_mag_r = GalaxyMagMatcher::get().match(density);
    }
    
    // Print the first few halo masses and galaxy assigned
    for (int i = 0; i < 5 && i < halos.size(); ++i) {
        std::cout << "Halo " << i << ": logmhalo = " << halos[i].logmhalo << ", abs_mag_r = " << halos[i].abs_mag_r << std::endl;
    }

    return 0;
}