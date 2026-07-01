import numpy as np
import sys
from pyutils import *
from dataloc import *
from footprintmanager import FootprintManager
import re 

def read_wp_file(fname):
    data = np.loadtxt(fname, skiprows=0, dtype='float')
    wp = data[:,1]
    wp_err = data[:,2]
    radius = data[:,0]
    return wp,wp_err,radius


class CalData:
    """
    Class to specify the calibration data for the a galaxy halo connection model.
    """
    def __init__(self, paramsfolder, magbins: np.ndarray, gmrbins: np.ndarray, magcut: float, fsky: float):
        self.paramsfolder = paramsfolder
        self.rpbinsfile = os.path.join(self.paramsfolder, 'wp_rbins.dat')
        self.magbins = magbins # absolute magnitude bin definitions (edges)
        self.gmrbins = gmrbins # g-r color bin definitions (edges)

        if self.magbins is None:
            self.set_bins_from_folder()

        self.magcut = magcut # apparent magnitude cut to use when calculating the volumes. Negelct .04 difference for N vs S
        self.zmaxes = np.array([get_max_observable_z(m, self.magcut) for m in self.magbins[:-1]])
        self.measurement_volumes = np.array([get_volume_at_z(z, fsky) for z in self.zmaxes])
        self.bincount = len(self.magbins) - 1

        self._wp_red_cache = {}
        self._wp_blue_cache = {}
        self._wp_all_cache = {}

    def set_bins_from_folder(self):
        # Auto read files from folder. Like: wp_mag-22.6473to-21.9413_gmr0.9778to1.0755.data
        filename_pattern = re.compile(
            r"wp"
            r"_mag(?P<mag_range>[\d\.-]+to[\d\.-]+)"  #  magnitude range, needs to handle negative sign too
            r"_gmr(?P<gr_range>[\d\.-]+to[\d\.-]+)"  #  g-r range
            r"\.dat$"
        )

        # Read all files in the folder and extract magbins and gmrbins
        magbins_set = set()
        gmrbins_set = {} # Dictionary of sets for each magbin. Key is like '-22.6473to-21.9413'
        for fname in os.listdir(self.paramsfolder):
            match = filename_pattern.match(fname)
            if match:
                mag_range = match.group("mag_range")
                gr_range = match.group("gr_range")

                if mag_range:
                    mag_min, mag_max = map(float, mag_range.split("to"))
                    magbins_set.add(mag_min)
                    magbins_set.add(mag_max)

                if gr_range:
                    gr_min, gr_max = map(float, gr_range.split("to"))
                    if mag_range not in gmrbins_set:
                        gmrbins_set[mag_range] = set()
                    gmrbins_set[mag_range].add(gr_min)
                    gmrbins_set[mag_range].add(gr_max)

        # Sort and convert to numpy arrays
        self.magbins = np.array(sorted(magbins_set, reverse=True))
        self.gmrbins = []
        for mag_range in sorted(gmrbins_set.keys(), reverse=True):
            self.gmrbins.append(np.array(sorted(gmrbins_set[mag_range])))

        self.gmrbins.reverse()  # Reverse the list of gmrbins to match the order of magbins
        
        # Print off entire bin structure for verification
        print("Detected magbins:", self.magbins)
        for arr in self.gmrbins:
            print("Detected gmrbins for magbin:", arr)  

    # Old methods for 1st test compatibility
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
    
    def get_wp_new(self, mag_idx: float, gmr_idx: float):
        mag_min = self.magbins[mag_idx]
        mag_max = self.magbins[mag_idx + 1]
        gmr_min = self.gmrbins[mag_idx][gmr_idx]
        gmr_max = self.gmrbins[mag_idx][gmr_idx + 1]
        data = np.loadtxt(os.path.join(self.paramsfolder, f'wp_mag{mag_max:.4f}to{mag_min:.4f}_gmr{gmr_min:.4f}to{gmr_max:.4f}.dat'))
        r = data[:,0]
        wp = data[:,1]
        wp_err = data[:,2]

        wp_err = np.where(np.isnan(wp_err), abs(0.05*wp), wp_err) # TODO

        return wp, wp_err, r
    
    def mag_to_idx(self, mag: float):
        m = abs(mag)
        return np.asarray(self.magbins == -m).nonzero()[0][0]
    
    @staticmethod
    def BGS_Y1_4bin(magcut: float):
        fsky = FootprintManager().get_footprint("Y1", min_passes=1) / DEGREES_ON_SPHERE # Measurements were made on 1-pass footprint
        return CalData(PARAMS_BGSY1_FOLDER, np.array([-19, -20, -21, -22, -23]), None, magcut, fsky)
    
    @staticmethod
    def BGS_Y1_6bin(magcut: float):
        fsky = FootprintManager().get_footprint("Y1", min_passes=1) / DEGREES_ON_SPHERE # Measurements were made on 1-pass footprint
        return CalData(PARAMS_BGSY1_FOLDER, np.array([-17, -18, -19, -20, -21, -22, -23]), None, magcut, fsky)

    @staticmethod
    def BGS_Y3_Test1(magcut: float):
        fsky = FootprintManager().get_footprint("Y3", min_passes=1) / DEGREES_ON_SPHERE # Measurements were made on 1-pass footprint
        return CalData(
            '/mount/sirocco1/imw2293/GROUP_CAT/latsham/data/dr2_wp_1/', 
            None,
            None,
            magcut,
            fsky)

            #np.array([-22.6473, -21.9413, -21.488,  -21.091,  -20.7092, -20.3206, -19.9002, -19.4147, -18.8037, -17.9408, -16.4536]), 
