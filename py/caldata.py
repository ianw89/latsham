
import numpy as np
import re 
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import subprocess
import os
import pandas as pd

class BiasData:
    """
    Data and utilities for projected clustering bias measurements.

    Built from a .npy numpy array file that should contain data, 1 row per sample, with the following fields:
    1. magbin_fainter - fainter bound of the magnitude bin
    2. magbin_brighter - brighter bound of the magnitude bin
    3. magbin_mean - mean magnitude of the sample
    4. quiescent - bool indicating if the sample is quiescent or false for star-forming
    5. third_property_name - optional third property name (nan if not applicable)
    6. third_propety_lower - optional lower bound of the third property (nan if not applicable)
    7. third_propety_upper - optional upper bound of the third property (nan if not applicable)
    8. third_property_mean - optional mean value of the sample for the third property (nan if not applicable)
    9. b - measured overall bias for the magnitude bin
    10. b_err_low - lower error bound for the measured bias
    11. b_err_high - upper error bound for the measured bias
    12. b_p - measured property bias
    13. b_p_err_low - lower error bound for the measured property bias
    14. b_p_err_high - upper error bound for the measured property bias
    15. b_m - magnitude bias used in calculating b_p from b
    """

    def __init__(self, datafolder: str):
        self.datafolder = datafolder
        self.data = self.load_biases()

    @staticmethod
    def BGS_DR2_5_30Mpc_v0_3():
        return BiasData('/mount/sirocco1/imw2293/GROUP_CAT/latsham/data/dr2_wp_2/')

    def load_biases(self):
        fn = 'biases_BGS_DR2_5-30Mpc_v0.3.npy'
        path = os.path.join(self.datafolder, fn)
        if not os.path.exists(path):
            raise FileNotFoundError(f"Bias data file not found: {path}. SCP it from NERSC with this command:\n scp ianw89@saul.nersc.gov:/global/cfs/cdirs/desi/users/ianw89/clustering712/DA2/LSS/loa-v1/LSScats/v1.1/v0.3/export_bias/{fn} /mount/sirocco1/imw2293/GROUP_CAT/latsham/data/dr2_wp_2/")

        loaded = np.load(path)
        return loaded

    def write_data_as_txt(self):
        if (self.data is not None):
            df = pd.DataFrame(self.data)
            txt_path = os.path.join(self.datafolder, 'biases_BGS_DR2_5-30Mpc_v0.3.txt')
            df.to_csv(txt_path, sep=' ', index=False)
        else:
            print("No data available to write.")

    
    def get(self, mag: float, quiescent: bool, propname: str = None, propvalue: float = None):
        """
        Get the bias for a given magnitude and quiescent/star-forming status, optionally with a third property.
        """
        # Filter by magnitude bin
        mag_mask = (self.data['magbin_fainter'] >= mag) & (mag >= self.data['magbin_brighter'])
        filtered_data = self.data[mag_mask]

        # Filter by quiescent/star-forming
        filtered_data = filtered_data[filtered_data['quiescent'] == quiescent]

        # If a third property is specified, filter by that as well
        if propname and propvalue is not None:
            filtered_data = filtered_data[filtered_data['third_property_name'] == propname]

            filtered_data = filtered_data[
                (filtered_data['third_property_lower'] <= propvalue) &
                (propvalue <= filtered_data['third_property_upper'])
            ]
        # If not, filter for entries where the third property is NaN
        else:
            filtered_data = filtered_data[np.isnan(filtered_data['third_property_name'])]

        if len(filtered_data) == 0:
            raise ValueError("No matching bias data found for the specified parameters.")

        # There should only be one match
        if len(filtered_data) > 1:
            raise ValueError("Multiple matching bias data entries found. Please check the input parameters.")

        match = filtered_data[0]

        # If requested mag is away from the mag mean of this bias measurement by more than 0.1 mag, raise a warning
        mag_mean = match['magbin_mean']
        if abs(mag - mag_mean) > 0.1:
            print(f"Warning: Requested magnitude {mag} is away from the mean magnitude {mag_mean} of the bias measurement by more than 0.1 mag.")

        # If requeseted third property value is away from the mean of this bias measurement by more than 0.1, raise a warning
        if propname and propvalue is not None:
            measured_prop_mean = match['third_property_mean']
            if not np.isclose(propvalue, measured_prop_mean, rtol=0.05):
                print(f"Warning: Requested third property value {propvalue} is away from the mean {measured_prop_mean} of the bias measurement by more than 5%.")

        # Assuming there's only one matching entry after filtering
        return match['b'], match['b_err_low'], match['b_err_high']
