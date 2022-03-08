"""Compare the WODEN outputs to hyperbeam. Need to run the ctest and
copy the text outputs to the directory, something like

cp ../../build/cmake_testing/FEE_primary_beam_cuda/*double*.txt .
"""

from os import environ
##Change hyperbeam env variable to point to interpolated one
environ['MWA_BEAM_FILE'] = environ['MWA_FEE_HDF5_INTERP']

import mwa_hyperbeam
import numpy as np
import erfa

zenith_delays = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
off_zenith1_delays = [0, 4, 8, 12, 0, 4, 8, 12, 0, 4, 8, 12, 0, 4, 8, 12]
off_zenith2_delays = [0, 2, 4, 8, 2, 4, 8, 12, 4, 8, 12, 16, 8, 12, 16, 20]

all_delays = [zenith_delays, off_zenith1_delays, off_zenith2_delays]
all_delay_names = ["zenith", "offzen1", "offzen2"]

MWA_LAT_RAD = -0.46606083776035967
DD2R = np.pi / 180.0

def reorder_jones(jones):

    reordered_jones = np.empty(jones.shape, dtype=complex)
    reordered_jones[0] =  -jones[3]
    reordered_jones[1] =  jones[2]
    reordered_jones[2] =  -jones[1]
    reordered_jones[3] =  jones[0]

    return reordered_jones

def rotate_jones(jones, para_ang):
    rot_jones = np.empty(jones.shape, dtype=complex)
    sinrot = np.sin(para_ang + np.pi/2.0)
    cosrot = np.cos(para_ang + np.pi/2.0)

    rot_jones[0] = -jones[1]*sinrot + jones[0]*cosrot
    rot_jones[1] = jones[1]*cosrot + jones[0]*sinrot
    rot_jones[2] = -jones[3]*sinrot + jones[2]*cosrot
    rot_jones[3] = jones[3]*cosrot + jones[2]*sinrot

    return rot_jones

def run_frequency_set(beam, azs, zas, freqs, freq_tag):
    """With the give hyperbeam, run all freqs and az,za coords for each
    set of delays"""

    for delays, delay_name in zip(all_delays, all_delay_names):
        with open(f"hyperbeam_multi_{delay_name}_{freq_tag}_rot.txt", 'w+') as outfile:
            for freq in freqs:
                # print(freq)
                all_jones = beam.calc_jones_array(azs, zas, freq, delays, gains, True)

                for ind in np.arange(len(azs)):

                    el = np.pi/2 - zas[ind];

                    ha, dec = para_ang = erfa.ae2hd(azs[ind], el, MWA_LAT_RAD)
                    para_angle = erfa.hd2pa(ha, dec, MWA_LAT_RAD)

                    reordered_jones = reorder_jones(all_jones[ind,:])
                    rotated_jones = rotate_jones(reordered_jones, para_angle)

                    gx = rotated_jones[0]
                    Dx = rotated_jones[1]
                    Dy = rotated_jones[2]
                    gy = rotated_jones[3]

                    outfile.write(f"{azs[ind]:.16f} {zas[ind]:.16f} "
                                  f"{np.real(gx):.16f} {np.imag(gx):.16f} {np.real(Dx):.16f} {np.imag(Dx):.16f} "
                                  f"{np.real(Dy):.16f} {np.imag(Dy):.16f} {np.real(gy):.16f} {np.imag(gy):.16f}\n")

if __name__ == '__main__':

    gains = [1]*16

    beam = mwa_hyperbeam.FEEBeam()

    azs = [270.0*DD2R, 270.0*DD2R, 0.0*DD2R, 90.0*DD2R, 90.0*DD2R]
    zas = [60.0*DD2R, 30.0*DD2R, 0.0*DD2R, 30.0*DD2R, 60.0*DD2R]

    freqs1 = np.arange(167e+6, 167e+6 + 32*80e+3, 80e+3)
    freqs2 = np.arange(167e+6, 167e+6 + 12*1.28e+6, 1.28e+6)
    freqs3 = np.arange(190e+6, 190e+6 + 12*320e+3, 320e+3)

    run_frequency_set(beam, azs, zas, freqs1, "freqs1")
    run_frequency_set(beam, azs, zas, freqs2, "freqs2")
    run_frequency_set(beam, azs, zas, freqs3, "freqs3")

    # print(environ['MWA_FEE_HDF5_INTERP'])
