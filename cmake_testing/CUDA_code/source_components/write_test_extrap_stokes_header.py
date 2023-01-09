import numpy as np
np.random.seed(4576)

num_powers = 6
num_curves = 6
num_lists = 6
num_extrap_freqs = 25

##=======POWER_LAW==============================================================

ref_stokesI = np.random.uniform(1e-3, 100.0, num_powers)
ref_stokesQ = np.random.uniform(1e-3, 100.0, num_powers)
ref_stokesU = np.random.uniform(1e-3, 100.0, num_powers)
ref_stokesV = np.random.uniform(1e-3, 100.0, num_powers)

extrap_freqs = np.linspace(50e+6, 300e+6, num_extrap_freqs)
ref_freqs = np.array([50e+6, 100e+6, 150e+6, 200e+6, 250e+6, 300e+6])
ref_power_SIs = np.random.uniform(-1.5, 0.5, num_powers)

##=======CURVED_POWER_LAW=======================================================

##want to have peaks between 100 and 200MHz so we can visually see things
##in this test, so define curvature and peak freq, and calculate SI from that
peak_freqs = np.random.uniform(100, 200, num_curves)*1e+6
ref_qs = np.random.uniform(-2, 0.5, num_curves)
ref_curve_SIs = -2*ref_qs*np.log(peak_freqs)


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

def write_header_file(filename, arrays, array_names, precisions, formats):
    with open(filename, 'w') as outfile:
        outfile.write('#include "woden_precision_defs.h"\n\n')
        outfile.write(f"int num_powers = {num_powers};\n")
        outfile.write(f"int num_curves = {num_curves};\n")
        outfile.write(f"int num_lists = {num_lists};\n")
        outfile.write(f"int num_extrap_freqs = {num_extrap_freqs};\n")

        for array, array_name, precision, format in zip(arrays, array_names, precisions, formats):

            arr_string = f"{precision} {array_name}[] = \u007b"
            for ind in range(len(array)-1):
                if format == 'f':
                    arr_string += f'{array[ind]:.16f},'
                elif format == "e":
                    arr_string += f'{array[ind]:.16e},'
                elif format == "d":
                    arr_string += f'{array[ind]:d},'

            if format == 'f':
                arr_string += f'{array[-1]:.16f}}};\n'
            elif format == 'e':
                arr_string += f'{array[-1]:.16e}}};\n'
            elif format == 'd':
                arr_string += f'{array[-1]:d}}};\n'

            outfile.write(arr_string)

    return


if __name__ == '__main__':

    filename = "test_extrap_stokes.h"
    arrays = [ref_stokesI, ref_stokesQ, ref_stokesU, ref_stokesV,
              extrap_freqs, ref_freqs, ref_power_SIs, ref_curve_SIs, ref_qs,
              num_list_values, list_start_indexes, list_freqs, list_stokesI,
              list_stokesQ, list_stokesU, list_stokesV]
    array_names = ["ref_stokesI", "ref_stokesQ", "ref_stokesU", "ref_stokesV",
              "extrap_freqs", "ref_freqs", "ref_power_SIs", "ref_curve_SIs", "ref_qs",
              "num_list_values", "list_start_indexes", "list_freqs", "list_stokesI",
              "list_stokesQ", "list_stokesU", "list_stokesV"]
    precisions = ["user_precision_t", "user_precision_t", "user_precision_t",
                  "user_precision_t", "double", "double", "user_precision_t",
                  "user_precision_t", "user_precision_t", "int", "int", "double",
                  "user_precision_t", "user_precision_t", "user_precision_t",
                  "user_precision_t"]
    formats = ["f", "f", "f", "f", "e", "e", "f", "f", "f", "d", "d", "e",
               "f", "f", "f", "f"]

    write_header_file(filename, arrays, array_names, precisions, formats)
