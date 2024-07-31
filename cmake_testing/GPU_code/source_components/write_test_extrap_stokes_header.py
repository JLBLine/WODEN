import numpy as np
np.random.seed(34)

num_powers = 6
num_curves = 6
num_lists = 6
num_pol_frac = 6
num_extrap_freqs = 100

##=======POWER_LAW==============================================================

ref_stokesI = np.random.uniform(1e-3, 10.0, num_powers)
ref_linpol = np.random.uniform(1e-3, 10.0, num_powers)
ref_stokesV = np.random.uniform(1e-3, 10.0, num_powers)

stokesI_power_SIs = np.random.uniform(-1.5, 0.5, num_powers)
linpol_power_SIs = np.random.uniform(-1.5, 0.5, num_powers)
stokesV_power_SIs = np.random.uniform(-1.5, 0.5, num_powers)

extrap_freqs = np.linspace(50e+6, 300e+6, num_extrap_freqs)


##=======CURVED_POWER_LAW=======================================================

##values pulled from model fits to LoBES data
stokesI_qs = np.array([-1.1540648022374427,-1.2848486571584865,-8.442091275232713 ,-7.984679649768583 ,1.4366766327372646 ,1.439801523394266])
stokesI_curve_SIs = np.array([3.9805869816599606, 1.7473922337782646,-0.7152342578625103,-4.432750841164616,-0.910260545384106,-0.5190814949517933])

##just reorder them for the other polarisations
linpol_qs = stokesI_qs[np.array([2,5,3,4,1,0])]
linpol_curve_SIs = stokesI_curve_SIs[np.array([2,5,3,4,1,0])]
stokesV_qs = stokesI_qs[np.array([5,2,3,0,1,4])]
stokesV_curve_SIs = stokesI_curve_SIs[np.array([5,2,3,0,1,4])]

##========POLARISATION FRACTION=================================================

linpol_pol_fracs = np.random.uniform(-1.0, 1.0, num_pol_frac) 
stokesV_pol_fracs = np.random.uniform(-1.0, 1.0, num_pol_frac)

##========RM and all that=======================================================

intr_pol_angle = np.random.uniform(0, 2*np.pi, num_powers + num_pol_frac + num_curves) 
##these are probably big RMs, but it's just a test, so chill
# rms = np.random.uniform(10, 100, num_powers + num_pol_frac + num_curves)
# rms = np.random.uniform(10, 30, num_powers + num_pol_frac + num_curves)
rms = np.random.uniform(0,1, num_powers + num_pol_frac + num_curves)

##=======LIST_STUFF=============================================================

num_list_values = np.random.uniform(5, 32, num_lists).astype(int)
total_vals = int(np.sum(num_list_values))

list_start_indexes = np.empty(num_lists, dtype=int)

list_freqs = np.empty(total_vals)
list_stokesI = np.empty(total_vals)
list_stokesQ = np.empty(total_vals)
list_stokesU = np.empty(total_vals)
list_stokesV = np.empty(total_vals)

start_index = 0

for ind, num_list_vals in enumerate(num_list_values):

    list_freqs[start_index:start_index+num_list_vals] = np.linspace(50e+6, 300e+6, num_list_vals)

    list_stokesI[start_index:start_index+num_list_vals] = np.random.uniform(-1,10,num_list_vals)
    list_stokesQ[start_index:start_index+num_list_vals] = np.random.uniform(-1,10,num_list_vals)
    list_stokesU[start_index:start_index+num_list_vals] = np.random.uniform(-1,10,num_list_vals)
    list_stokesV[start_index:start_index+num_list_vals] = np.random.uniform(-1,10,num_list_vals)

    list_start_indexes[ind] = start_index

    start_index += num_list_vals

def write_header_file(filename, arrays, array_names, precisions):
    with open(filename, 'w') as outfile:
        outfile.write('#include "woden_precision_defs.h"\n\n')
        outfile.write(f"int num_powers = {num_powers};\n")
        outfile.write(f"int num_curves = {num_curves};\n")
        outfile.write(f"int num_lists = {num_lists};\n")
        outfile.write(f"int num_pol_frac = {num_pol_frac};\n")
        outfile.write(f"int num_extrap_freqs = {num_extrap_freqs};\n")

        for array, array_name, precision in zip(arrays, array_names, precisions):

            arr_string = f"{precision} {array_name}[] = \u007b"
            for ind in range(len(array)-1):
                if precision == "user_precision_t":
                    arr_string += f'{array[ind]:.16f},'
                elif precision == "double":
                    arr_string += f'{array[ind]:.16e},'
                elif precision == "int":
                    arr_string += f'{array[ind]:d},'

            if precision == "user_precision_t":
                arr_string += f'{array[-1]:.16f}}};\n'
            elif precision == "double":
                arr_string += f'{array[-1]:.16e}}};\n'
            elif precision == "int":
                arr_string += f'{array[-1]:d}}};\n'

            outfile.write(arr_string)

    return


if __name__ == '__main__':

    filename = "test_extrap_stokes.h"
    arrays = [ref_stokesI, ref_linpol, ref_stokesV,
              extrap_freqs, stokesI_power_SIs, linpol_power_SIs, stokesV_power_SIs,
              stokesI_qs, stokesI_curve_SIs, linpol_qs, linpol_curve_SIs,
              stokesV_qs, stokesV_curve_SIs,
              num_list_values, list_start_indexes, list_freqs, list_stokesI,
              list_stokesQ, list_stokesU, list_stokesV,
              linpol_pol_fracs, stokesV_pol_fracs, intr_pol_angle, rms]
    array_names = ["ref_stokesI", "ref_linpol", "ref_stokesV",
                   "extrap_freqs", "stokesI_power_SIs", "linpol_power_SIs", "stokesV_power_SIs",
                   "stokesI_qs", "stokesI_curve_SIs", "linpol_qs", "linpol_curve_SIs",
                   "stokesV_qs", "stokesV_curve_SIs",
                   "num_list_values", "list_start_indexes", "list_freqs", "list_stokesI",
                   "list_stokesQ", "list_stokesU", "list_stokesV",
                   "linpol_pol_fracs", "stokesV_pol_fracs", "intr_pol_angle", "rms"]
    precisions = ["user_precision_t", "user_precision_t", "user_precision_t",
                   "double", "user_precision_t", "user_precision_t", "user_precision_t",
                   "user_precision_t", "user_precision_t", "user_precision_t", "user_precision_t",
                   "user_precision_t", "user_precision_t",
                   "int", "int", "double", "user_precision_t",
                   "user_precision_t", "user_precision_t", "user_precision_t",
                   "user_precision_t", "user_precision_t",
                   "user_precision_t", "user_precision_t"]
    
    write_header_file(filename, arrays, array_names, precisions)
