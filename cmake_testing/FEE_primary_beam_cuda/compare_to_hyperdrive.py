"""Compare the WODEN outputs to hyperbeam. Need to run the ctest and
copy the text outputs to the directory, something like

cp ../../build/cmake_testing/FEE_primary_beam_cuda/*double*.txt .
"""

import mwa_hyperbeam
import numpy as np
import erfa


zenith_delays = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
off_zenith1_delays = [0, 4, 8, 12, 0, 4, 8, 12, 0, 4, 8, 12, 0, 4, 8, 12]
off_zenith2_delays = [0, 2, 4, 8, 2, 4, 8, 12, 4, 8, 12, 16, 8, 12, 16, 20]

all_delays = [zenith_delays, off_zenith1_delays, off_zenith2_delays]
all_delay_names = ["zenith", "offzen1", "offzen2"]

freqs = [100e+6, 150e+6, 200e+6]

MWA_LAT_RAD = -0.46606083776035967

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

if __name__ == '__main__':
    beam = mwa_hyperbeam.FEEBeam()

    gains = [1]*16

    for freq in freqs:
        for delays, delay_name in zip(all_delays, all_delay_names):

            # woden_data = np.loadtxt(f"{delay_name}_{int(freq/1e+6)}_double.txt")
            #
            # azs = woden_data[:,0]
            # zas = woden_data[:,1]

            # with open(f"az-za_coords.txt", 'w+') as outfile:
            #     for ind in np.arange(len(azs)):
            #         outfile.write(f"{azs[ind]:.16f} {zas[ind]:.16f}\n")

            azs, zas = np.loadtxt("az-za_coords.txt", unpack=True)
            all_jones = beam.calc_jones_array(azs, zas, freq, delays, gains, True)

            # with open(f"hyperbeam_{delay_name}_{int(freq/1e+6)}.txt", 'w+') as outfile:
            #     for ind in np.arange(len(azs)):
            #         reordered_jones = reorder_jones(all_jones[ind,:])
            #
            #         gx = reordered_jones[0]
            #         Dx = reordered_jones[1]
            #         Dy = reordered_jones[2]
            #         gy = reordered_jones[3]
            #
            #         outfile.write(f"{azs[ind]:.16f} {zas[ind]:.16f} "
            #                       f"{np.real(gx):.16f} {np.imag(gx):.16f} {np.real(Dx):.16f} {np.imag(Dx):.16f} "
            #                       f"{np.real(Dy):.16f} {np.imag(Dy):.16f} {np.real(gy):.16f} {np.imag(gy):.16f}\n")

            with open(f"hyperbeam_{delay_name}_{int(freq/1e+6)}_rot.txt", 'w+') as outfile:
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
