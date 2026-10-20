import numpy as np
from math import log, exp, sqrt, cos, sin, pi
import os
import subprocess as sp
import time

from caldata import CalData

def spherical_to_linear_parameters(hyperparams):
    if len(hyperparams) == 6:
        theta1a, theta1b, theta1c, theta2a, theta2b, theta2c = hyperparams
        return[cos(theta1a), sin(theta1a)*cos(theta1b), sin(theta1a)*sin(theta1b)*cos(theta1c), sin(theta1a)*sin(theta1b)*sin(theta1c),
                    cos(theta2a), sin(theta2a)*cos(theta2b), sin(theta2a)*sin(theta2b)*cos(theta2c), sin(theta2a)*sin(theta2b)*sin(theta2c)]
    if len(hyperparams) == 2:
        theta1, theta2 = hyperparams
        return [cos(theta1), sin(theta1), cos(theta2), sin(theta2)]
    else:        
        raise ValueError("Hyperparams must be either 2 or 6 elements long.")
    

def project_to_hypersphere(vec, radius=1.0):
    norm = sqrt(sum(v*v for v in vec))
    if norm == 0:
        raise ValueError("Cannot project the zero vector.")
    return [radius * v / norm for v in vec]


def linear_to_spherical_parameters(linear_params):
    if len(linear_params) == 8:
        x1, y1, z1, w1, x2, y2, z2, w2 = linear_params

        vec1 = project_to_hypersphere([x1, y1, z1, w1])
        vec2 = project_to_hypersphere([x2, y2, z2, w2])

        x1, y1, z1, w1 = vec1
        x2, y2, z2, w2 = vec2

        theta1a = np.arccos(x1)
        theta1b = np.arctan2(np.sqrt(z1*z1 + w1*w1), y1)
        theta1c = np.arctan2(w1, z1)

        theta2a = np.arccos(x2)
        theta2b = np.arctan2(np.sqrt(z2*z2 + w2*w2), y2)
        theta2c = np.arctan2(w2, z2)

        # Make them all positive and in the range [0, 2*pi) (some are really in [0, pi] and some are in [0, 2*pi])
        theta1a = theta1a % (2 * np.pi)
        theta1b = theta1b % (2 * np.pi)
        theta1c = theta1c % (2 * np.pi)
        theta2a = theta2a % (2 * np.pi)
        theta2b = theta2b % (2 * np.pi)
        theta2c = theta2c % (2 * np.pi)

        return [
            theta1a, theta1b, theta1c,
            theta2a, theta2b, theta2c
        ]
    if len(linear_params) == 4:
        x1, y1, x2, y2 = linear_params

        # Project the 2D vectors onto the unit circle
        vec1 = project_to_hypersphere([x1, y1])
        vec2 = project_to_hypersphere([x2, y2])

        x1, y1 = vec1
        x2, y2 = vec2

        theta1 = np.arccos(x1)
        theta2 = np.arccos(x2)

        return [theta1, theta2]
    else:
        raise ValueError("Linear params must be either 4 or 8 elements long.")
    

assert spherical_to_linear_parameters(linear_to_spherical_parameters([1, 0, 0, 0, 1, 0, 0, 0])) == [1, 0, 0, 0, 1, 0, 0, 0]


angles = np.random.rand(6)
angles[0] *= np.pi
angles[1] *= np.pi
angles[2] *= 2*np.pi
angles[3] *= np.pi
angles[4] *= np.pi
angles[5] *= 2*np.pi

lin = spherical_to_linear_parameters(angles)
angles2 = linear_to_spherical_parameters(lin)
lin2 = spherical_to_linear_parameters(angles2)

assert(np.max(np.abs(np.array(lin) - np.array(lin2))) < 1e-10)



def get_mock_wp(path):
    output_mock = np.loadtxt(path)
    wp_m = output_mock[:,4]
    return wp_m

def calc_wp_from_mocks(caldata: CalData):
    t1 = time.time()
    
    # Define the path to the Corrfunc binary
    corrfunc_path = "/export/sirocco1/tinker/src/Corrfunc/bin/"
    if not os.path.exists(os.path.join(corrfunc_path, "wp")):
        corrfunc_path = "/mount/sirocco1/tinker/src/Corrfunc/bin/"

    output_folder = "/mount/sirocco1/imw2293/GROUP_CAT/OUTPUT/LATSHAM/"
    nthreads = os.cpu_count()
    pimax = 40
    boxsize = 400 # Mpc/h

    # Force a sync to ensure the mock files from the C++ process are fully written to disk
    # before we try to read them with corrfunc. This prevents a potential race condition.
    sp.run(["sync"], check=True)

    if not os.path.exists(f"{corrfunc_path}/wp"):
        print(f"Corrfunc wp executable not found at {corrfunc_path}/wp. Skipping wp calculation.")
        return

    # Call corrfunc compiled executable to compute wp on the mocks from latsham
    for i in range(len(caldata.magbins) - 1):
        for j in range(len(caldata.gmrbins[i]) - 1):
            mag_min = caldata.magbins[i]
            mag_max = caldata.magbins[i + 1]
            gmr_min = caldata.gmrbins[i][j]
            gmr_max = caldata.gmrbins[i][j + 1]
            mock_fname = f"mock_{mag_min:.4f}_{mag_max:.4f}_{gmr_min:.4f}_{gmr_max:.4f}.dat"
            wp_output_fname = f"wp_mock_{mag_min:.4f}_{mag_max:.4f}_{gmr_min:.4f}_{gmr_max:.4f}.dat"
            mock_path = os.path.join(output_folder, mock_fname)
            wp_output_path = os.path.join(output_folder, wp_output_fname)

            cmd = f"{corrfunc_path}/wp {boxsize} {mock_path} a {caldata.rpbinsfile} {pimax} {nthreads} > {wp_output_path} 2> wp_stderr.txt"
            result = sp.run(cmd, cwd=output_folder, shell=True, check=True)
            if result.returncode != 0:
                print(f"Error running corrfunc!")

    t2 = time.time()

    print(f"Calculated all wp in {t2 - t1:.2f} seconds.")

def get_mock_err(volume, wp_err, wp_m):
    vfac = (volume / 400.0**3)**.5 # factor by which to multiply errors
    efac = 0.1 # TODO tweak eventually
    wp_m_err = vfac*wp_err + efac*wp_m
    return wp_m_err

def chisqr_direct_wp(caldata: CalData):
    IDX_START = 3 # Let's skip first 3 bins due to noise right now
    
    chivec = []
    for i in range(len(caldata.magbins) - 1):
        for j in range(len(caldata.gmrbins[i]) - 1):
            mag_min = caldata.magbins[i]
            mag_max = caldata.magbins[i + 1]
            gmr_min = caldata.gmrbins[i][j]
            gmr_max = caldata.gmrbins[i][j + 1]
                    
            # Read in the wp data from the data
            wp, wp_err, r = caldata.get_wp_new(i, j)
            r = r[IDX_START:]
            wp = wp[IDX_START:]
            wp_err = wp_err[IDX_START:]

            # mock
            wp_m = get_mock_wp(f"/mount/sirocco1/imw2293/GROUP_CAT/OUTPUT/LATSHAM/wp_mock_{mag_min:.4f}_{mag_max:.4f}_{gmr_min:.4f}_{gmr_max:.4f}.dat")
            wp_m = wp_m[IDX_START:]
            # If anything in nan, immediatly make these params very bad chisqr
            if np.isnan(wp_m).any():
                print(f"NaN found in wp for magbin {mag_min}-{mag_max}, color {gmr_min}-{gmr_max}. Assigning infinite chi^2.")
                return np.inf
            wp_m_err = get_mock_err(caldata.measurement_volumes[i], wp_err, wp_m)

            sigma = wp_err**2 + wp_m_err**2
            chivec.append(np.sum(((wp - wp_m) / sigma)**2))

    chi2 = np.sum(chivec)
    print(f"Chi^2: {chi2:.2f}")
    return chi2

def chisqr_on_offsets(caldata: CalData):

    chivec = []
    for i in range(len(caldata.magbins) - 1):
        for j in range(len(caldata.gmrbins[i]) - 1):
            # data
            # TODO copy over 
            data_offset, data_offset_err = caldata.get_offset(i, j)

            # mock
            mag_min = caldata.magbins[i]
            mag_max = caldata.magbins[i + 1]
            gmr_min = caldata.gmrbins[i][j]
            gmr_max = caldata.gmrbins[i][j + 1]
            wp_m = get_mock_wp(f"/mount/sirocco1/imw2293/GROUP_CAT/OUTPUT/LATSHAM/wp_mock_{mag_min:.4f}_{mag_max:.4f}_{gmr_min:.4f}_{gmr_max:.4f}.dat")
            
            # TODO what about mock cov matrix?
            mock_offset, mock_offset_err = caldata.compute_wp_offset(wp_m)