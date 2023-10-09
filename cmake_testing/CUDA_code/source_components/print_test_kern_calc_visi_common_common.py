import numpy as np
np.random.seed(3249085)

num_curves = 10
num_lists = 5

##=======CURVED_POWER_LAW=======================================================

ref_stokesI = np.random.uniform(1e-3, 100.0, num_curves)
ref_stokesQ = np.zeros(num_curves)
ref_stokesU = np.zeros(num_curves)
ref_stokesV = np.zeros(num_curves)

ref_freqs = np.linspace(50e+6, 300e+6, num_curves)

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

curve_inds = range(num_curves, 2*num_curves)
list_inds = range(2*num_curves, 2*num_curves + num_lists)

start_index = 0

for ind, num_list_vals in enumerate(num_list_values):

    list_freqs[start_index:start_index+num_list_vals] = np.linspace(50e+6, 300e+6, num_list_vals)

    list_stokesI[start_index:start_index+num_list_vals] = np.random.uniform(-5,10,num_list_vals)
    list_stokesQ[start_index:start_index+num_list_vals] = 0;
    list_stokesU[start_index:start_index+num_list_vals] = 0;
    list_stokesV[start_index:start_index+num_list_vals] = 0;

    list_start_indexes[ind] = start_index

    start_index += num_list_vals

def print_header_file(arrays, array_names, precisions, formats):

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

        print(arr_string)

    print(f"int total_num_flux_entires = {int(np.sum(num_list_values)):d};")

    return


if __name__ == '__main__':

    arrays = [ref_stokesI, ref_stokesQ, ref_stokesU, ref_stokesV,
              ref_freqs, ref_curve_SIs, ref_qs,
              num_list_values, list_start_indexes, list_freqs, list_stokesI,
              list_stokesQ, list_stokesU, list_stokesV, curve_inds, list_inds]
    array_names = ["k_ref_stokesI", "k_ref_stokesQ", "k_ref_stokesU", "k_ref_stokesV",
                  "k_ref_freqs", "k_ref_curve_SIs", "k_ref_qs",
                  "k_num_list_values", "k_list_start_indexes", "k_list_freqs", "k_list_stokesI",
                  "k_list_stokesQ", "k_list_stokesU", "k_list_stokesV",
                  "k_curve_comp_inds", "k_list_comp_inds"]
    precisions = ["user_precision_t", "user_precision_t", "user_precision_t",
                  "user_precision_t", "double",                   "user_precision_t", "user_precision_t", "int", "int", "double",
                  "user_precision_t", "user_precision_t", "user_precision_t",
                  "user_precision_t", "int", "int"]
    formats = ["f", "f", "f", "f", "e", "f", "f", "d", "d", "e",
               "f", "f", "f", "f","d", "d"]

    print_header_file(arrays, array_names, precisions, formats)
