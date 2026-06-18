import numpy as np
import sys

from pyutils import *
from dataloc import *
from footprintmanager import FootprintManager


def read_wp_file(fname):
    data = np.loadtxt(fname, skiprows=0, dtype='float')
    wp = data[:,1]
    wp_err = data[:,2]
    radius = data[:,0]
    return wp,wp_err,radius

class CalData:
    """
    Class to specify the calibration data for the SelfCalGroupFinder.
    """
    def __init__(self, paramsfolder, magbins: np.ndarray, color_separation: np.ndarray, magcut):
        if len(color_separation) != len(magbins) - 1:
            raise ValueError("color_separation must have one less element than magbins")
        
        self.paramsfolder = paramsfolder
        self.rpbinsfile = os.path.join(self.paramsfolder, 'wp_rbins.dat')
        self.magbins = magbins # absolute magnitude bin definitions (edges)
        self.color_separation = color_separation # boolean array indicating whether for each L bin the wp should be red/blue seperate (True), or all together (False)
        self.magcut = magcut # apparent magnitude cut to use when calculating the volumes. Negelct .04 difference for N vs S
        self.zmaxes = np.array([get_max_observable_z(m, self.magcut) for m in self.magbins[:-1]])
        fsky = FootprintManager().get_footprint("Y1", min_passes=1) / DEGREES_ON_SPHERE # Measurements were made on 1-pass footprint
        self.measurement_volumes = np.array([get_volume_at_z(z, fsky) for z in self.zmaxes])
        self.bincount = len(self.magbins) - 1

        self._wp_red_cache = {}
        self._wp_blue_cache = {}
        self._wp_all_cache = {}

    def get_wp_red(self, mag: int):
        mag = abs(mag)
        if mag not in self._wp_red_cache:
            fname = os.path.join(self.paramsfolder, f'wp_red_M{mag:d}.dat')
            self._wp_red_cache[mag] = read_wp_file(fname)
        return self._wp_red_cache[mag]
    
    def get_wp_blue(self, mag: int):
        mag = abs(mag)
        if mag not in self._wp_blue_cache:
            fname = os.path.join(self.paramsfolder, f'wp_blue_M{mag:d}.dat')
            self._wp_blue_cache[mag] = read_wp_file(fname)
        return self._wp_blue_cache[mag]
    
    def get_wp_all(self, mag: int):
        mag = abs(mag)
        if mag not in self._wp_all_cache:
            fname = os.path.join(self.paramsfolder, f'wp_all_M{mag:d}.dat')
            self._wp_all_cache[mag] = read_wp_file(fname)
        return self._wp_all_cache[mag]
    
    def mag_to_idx(self, mag: float):
        m = abs(mag)
        return np.asarray(self.magbins == -m).nonzero()[0][0]
    
    @staticmethod
    def BGS_Y1_4bin(magcut: float):
        return CalData(PARAMS_BGSY1_FOLDER, np.array([-19, -20, -21, -22, -23]), np.array([True, True, True, True]), magcut)
    
    @staticmethod
    def BGS_Y1_6bin(magcut: float):
        return CalData(PARAMS_BGSY1_FOLDER, np.array([-17, -18, -19, -20, -21, -22, -23]), np.array([True, True, True, True, True, True]), magcut)

    def __str__(self):
        return f"CalibrationData(\nmagbins={self.magbins}\nzmaxes={self.zmaxes}\nmagcut={self.magcut}\nparamsfolder={self.paramsfolder}\nbinsfile={self.rpbinsfile}\n)"
    
