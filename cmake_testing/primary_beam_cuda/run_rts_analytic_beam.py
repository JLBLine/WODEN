import numpy as np
from os import environ
from subprocess import call, check_output

# BEAM_PATH = environ['MWA_FEE_HDF5']

D2R = np.pi / 180.0

coords = np.load('az-za_coords.npz')

azs = coords['az_arr']
zas = coords['za_arr']
nside = coords['nside']

# print(azs)

azs /= D2R
zas /= D2R


gx_res = []
Dx_res = []
Dy_res = []
gy_res = []


for ind, az, za in zip(range(len(azs)), azs, zas):

    print(f"Doing coord {ind} of {len(azs)}")

    returned_text = check_output(f"./rts_mwa_beam $MWA_FEE_HDF5 {az:.8f} {za:.8f} 200e+6 0", shell=True, universal_newlines=True)
    # print(returned_text)

    lines = returned_text.split('\n')

    chunks = lines[3].split()
    gx_re = float(chunks[0].strip("['j"))
    Dx_re = float(chunks[4].strip("['j"))

    chunks = lines[4].split()
    Dy_re = float(chunks[0].strip("['j"))
    gy_re = float(chunks[4].strip("['j"))

    # print(gx_re, Dx_re, Dy_re, gy_re)
    gx_res.append(gx_re)
    Dx_res.append(Dx_re)
    Dy_res.append(Dy_re)
    gy_res.append(gy_re)


save_array = np.zeros((len(gx_res), 11))

save_array[:,0] = azs*D2R
save_array[:,1] = zas*D2R
save_array[:,2] = gx_res
save_array[:,4] = Dx_res
save_array[:,6] = Dy_res
save_array[:,8] = gy_res
save_array[:,10] = 150e+6

np.savetxt(f"RTS_analytic_beam_gains_azza{nside*nside}.txt", save_array)
