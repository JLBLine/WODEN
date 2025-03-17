import numpy as np
np.random.seed(34)


def make_list_arrays(num_lists):
    
    num_list_values = np.random.uniform(5, 32, num_lists).astype(int)
    total_vals = int(np.sum(num_list_values))

    list_start_indexes = np.empty(num_lists, dtype=int)

    list_freqs = np.empty(total_vals)
    list_stokes = np.empty(total_vals)

    start_index = 0

    for ind, num_list_vals in enumerate(num_list_values):

        list_freqs[start_index:start_index+num_list_vals] = np.linspace(50e+6, 300e+6, num_list_vals)
        list_stokes[start_index:start_index+num_list_vals] = np.random.uniform(-1,10,num_list_vals)
        list_start_indexes[ind] = start_index
        start_index += num_list_vals

    return num_list_values, list_start_indexes, list_freqs, list_stokes

num_powers = 5
num_curves = 5
num_lists = 5
num_pol_frac = 3
num_pol_power = 3
num_pol_curv = 3
num_pol_list = 3

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
linpol_qs = stokesI_qs[np.array([4,1,0])]
linpol_curve_SIs = stokesI_curve_SIs[np.array([4,1,0])]
stokesV_qs = stokesI_qs[np.array([5,2,3])]
stokesV_curve_SIs = stokesI_curve_SIs[np.array([5,2,3])]

##========POLARISATION FRACTION=================================================

linpol_pol_fracs = np.random.uniform(-1.0, 1.0, num_pol_frac) 
stokesV_pol_fracs = np.random.uniform(-1.0, 1.0, num_pol_frac)

##========RM and all that=======================================================

intr_pol_angle = np.random.uniform(0, 2*np.pi, num_pol_power + num_pol_frac + num_pol_curv + num_pol_list) 
rms = np.random.uniform(0,1, num_pol_power + num_pol_frac + num_pol_curv + num_pol_list)

##=======LIST_STUFF=============================================================

num_list_values, list_start_indexes, list_freqs, list_stokesI = make_list_arrays(num_lists)

stokesQ_num_list_values, stokesQ_list_start_indexes, stokesQ_list_ref_freqs, stokesQ_list_ref_flux = make_list_arrays(num_pol_list)
stokesU_num_list_values, stokesU_list_start_indexes, stokesU_list_ref_freqs, stokesU_list_ref_flux = make_list_arrays(num_pol_list)
stokesV_num_list_values, stokesV_list_start_indexes, stokesV_list_ref_freqs, stokesV_list_ref_flux = make_list_arrays(num_pol_list)
linpol_p_num_list_values, linpol_p_list_start_indexes, linpol_p_list_ref_freqs, linpol_p_list_ref_flux = make_list_arrays(num_pol_list)



def write_header_file(filename, arrays, array_names, precisions):
    with open(filename, 'w') as outfile:
        outfile.write('#include "woden_precision_defs.h"\n\n')
        outfile.write(f"int num_powers = {num_powers};\n")
        outfile.write(f"int num_curves = {num_curves};\n")
        outfile.write(f"int num_lists = {num_lists};\n")
        outfile.write(f"int num_pol_frac = {num_pol_frac};\n")
        outfile.write(f"int num_extrap_freqs = {num_extrap_freqs};\n")
        outfile.write(f"int num_pol_power = {num_pol_power};\n")
        outfile.write(f"int num_pol_curv = {num_pol_curv};\n")
        outfile.write(f"int num_pol_list = {num_pol_list};\n")

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
              linpol_pol_fracs, stokesV_pol_fracs, intr_pol_angle, rms,
              stokesQ_num_list_values, stokesQ_list_start_indexes, stokesQ_list_ref_freqs, stokesQ_list_ref_flux,
              stokesU_num_list_values, stokesU_list_start_indexes, stokesU_list_ref_freqs, stokesU_list_ref_flux,
              stokesV_num_list_values, stokesV_list_start_indexes, stokesV_list_ref_freqs, stokesV_list_ref_flux,
              linpol_p_num_list_values, linpol_p_list_start_indexes, linpol_p_list_ref_freqs, linpol_p_list_ref_flux]
    array_names = ["ref_stokesI", "ref_linpol", "ref_stokesV",
                   "extrap_freqs", "stokesI_power_SIs", "linpol_power_SIs", "stokesV_power_SIs",
                   "stokesI_qs", "stokesI_curve_SIs", "linpol_qs", "linpol_curve_SIs",
                   "stokesV_qs", "stokesV_curve_SIs",
                   "num_list_values", "list_start_indexes", "list_freqs", "list_stokesI",
                   "linpol_pol_fracs", "stokesV_pol_fracs", "intr_pol_angle", "rms",
                   "stokesQ_num_list_values", "stokesQ_list_start_indexes", "stokesQ_list_ref_freqs", "stokesQ_list_ref_flux",
                   "stokesU_num_list_values", "stokesU_list_start_indexes", "stokesU_list_ref_freqs", "stokesU_list_ref_flux",
                   "stokesV_num_list_values", "stokesV_list_start_indexes", "stokesV_list_ref_freqs", "stokesV_list_ref_flux",
                   "linpol_p_num_list_values", "linpol_p_list_start_indexes", "linpol_p_list_ref_freqs", "linpol_p_list_ref_flux"]
    precisions = ["user_precision_t", "user_precision_t", "user_precision_t",
                   "double", "user_precision_t", "user_precision_t", "user_precision_t",
                   "user_precision_t", "user_precision_t", "user_precision_t", "user_precision_t",
                   "user_precision_t", "user_precision_t",
                   "int", "int", "double", "user_precision_t",
                   "user_precision_t", "user_precision_t",
                   "user_precision_t", "user_precision_t",
                   "int", "int", "double", "user_precision_t",
                   "int", "int", "double", "user_precision_t",
                   "int", "int", "double", "user_precision_t",
                   "int", "int", "double", "user_precision_t"]
    
    write_header_file(filename, arrays, array_names, precisions)
