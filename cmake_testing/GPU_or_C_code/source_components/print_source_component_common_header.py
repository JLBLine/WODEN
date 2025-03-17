import numpy as np
np.random.seed(5476)

num_powers = 3
num_curves = 3
num_lists = 3

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

##=======POWER_LAW==============================================================

ref_stokesI = np.random.uniform(1e-3, 100.0, num_powers)

ref_freqs = np.linspace(50e+6, 300e+6, num_powers)
ref_power_SIs = np.random.uniform(-1.5, 0.5, num_powers)

##=======CURVED_POWER_LAW=======================================================

##want to have peaks between 100 and 200MHz so we can visually see things
##in this test, so define curvature and peak freq, and calculate SI from that
peak_freqs = np.random.uniform(100, 200, num_curves)*1e+6
ref_qs = np.random.uniform(-2, 0.5, num_curves)
ref_curve_SIs = np.random.uniform(-1, 1, num_curves)


##=======LIST_STUFF=============================================================

num_list_values, list_start_indexes, list_freqs, list_stokesI = make_list_arrays(num_lists)

##polarisation stuff============================================================

n_stokesV_pol_frac = 2
n_stokesV_power = 2
n_stokesV_curve = 2
n_stokesV_list = 3
n_linpol_pol_frac = 2
n_linpol_power = 2
n_linpol_curve = 2
n_linpol_list = 2
n_linpol_p_list = 1

n_linpol_angles = n_linpol_pol_frac + n_linpol_power + n_linpol_curve + n_linpol_p_list

stokesV_pol_fracs = np.random.uniform(-1.0,  1.0, n_stokesV_pol_frac)
stokesV_pol_frac_comp_inds = np.array([0, 4])

ref_stokesV = np.random.uniform(1e-3, 10.0, n_stokesV_power)
stokesV_power_SIs = np.random.uniform(-1, 1, n_linpol_curve)
stokesV_power_comp_inds = np.array([1, 5])

stokesV_qs = np.random.uniform(-2, 0.5, n_stokesV_curve)
stokesV_curve_SIs = -2*stokesV_qs*np.log(peak_freqs[:n_stokesV_curve])
stokesV_curve_comp_inds = np.array([2, 6])

stokesV_num_list_values, stokesV_list_start_indexes, stokesV_list_ref_freqs, stokesV_list_ref_flux = make_list_arrays(n_stokesV_list)
stokesV_list_comp_inds = np.array([3, 7, 8])

linpol_pol_fracs = np.random.uniform(-1.0,  1.0, n_linpol_pol_frac)
linpol_pol_frac_comp_inds = np.array([0, 4])

ref_linpol = np.random.uniform(1e-3, 10.0, n_linpol_power)
linpol_power_SIs = np.random.uniform(-1.5, 0.5, n_linpol_power)
linpol_power_comp_inds = np.array([1, 5])

linpol_qs = np.random.uniform(-2, 0.5, n_linpol_curve)
linpol_curve_SIs = np.random.uniform(-1, 1, n_linpol_curve)
# linpol_curve_SIs = -2*linpol_qs*np.log(peak_freqs)
linpol_curve_comp_inds = np.array([2, 6])

stokesQ_num_list_values, stokesQ_list_start_indexes, stokesQ_list_ref_freqs, stokesQ_list_ref_flux = make_list_arrays(n_linpol_list)
stokesU_num_list_values, stokesU_list_start_indexes, stokesU_list_ref_freqs, stokesU_list_ref_flux = make_list_arrays(n_linpol_list)
stokesQ_list_comp_inds = np.array([3, 7])
stokesU_list_comp_inds = np.array([3, 7])

linpol_p_num_list_values, linpol_p_list_start_indexes, linpol_p_list_ref_freqs, linpol_p_list_ref_flux = make_list_arrays(n_linpol_p_list)
linpol_p_list_comp_inds = np.array([8])


intr_pol_angle = np.random.uniform(0, 2*np.pi, n_linpol_angles)
rms = np.random.uniform(0, 80, n_linpol_angles)
linpol_angle_inds = np.concatenate([linpol_pol_frac_comp_inds, linpol_power_comp_inds, linpol_curve_comp_inds, linpol_p_list_comp_inds])

def print_header_file(arrays, array_names, precisions):

    print(f"int num_powers = {num_powers};")
    print(f"int num_curves = {num_curves};")
    print(f"int num_lists = {num_lists};")
    print(f"int n_stokesV_pol_frac = {n_stokesV_pol_frac};")
    print(f"int n_stokesV_power = {n_stokesV_power};")
    print(f"int n_stokesV_curve = {n_stokesV_curve};")
    print(f"int n_stokesV_list = {n_stokesV_list};")
    print(f"int n_linpol_pol_frac = {n_linpol_pol_frac};")
    print(f"int n_linpol_power = {n_linpol_power};")
    print(f"int n_linpol_curve = {n_linpol_curve};")
    print(f"int n_linpol_list = {n_linpol_list};")
    print(f"int n_linpol_p_list = {n_linpol_p_list};")
    print(f"int n_stokesV_list_flux_entries = {int(stokesV_num_list_values.sum())};")
    print(f"int n_stokesQ_list_flux_entries = {int(stokesQ_num_list_values.sum())};")
    print(f"int n_stokesU_list_flux_entries = {int(stokesU_num_list_values.sum())};")
    print(f"int n_linpol_p_list_flux_entries = {int(linpol_p_num_list_values.sum())};")

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

        print(arr_string)

    return


if __name__ == '__main__':

    # arrays = [ref_stokesI, ref_stokesQ, ref_stokesU, ref_stokesV,
    #           ref_freqs, ref_power_SIs, ref_curve_SIs, ref_qs,
    #           num_list_values, list_start_indexes, list_freqs, list_stokesI,
    #           list_stokesQ, list_stokesU, list_stokesV]
    # array_names = ["ref_stokesI", "ref_stokesQ", "ref_stokesU", "ref_stokesV",
    #           "ref_freqs", "ref_power_SIs", "ref_curve_SIs", "ref_qs",
    #           "num_list_values", "list_start_indexes", "list_freqs", "list_stokesI",
    #           "list_stokesQ", "list_stokesU", "list_stokesV"]
    # precisions = ["user_precision_t", "user_precision_t", "user_precision_t",
    #               "user_precision_t", "double", "user_precision_t",
    #               "user_precision_t", "user_precision_t", "int", "int", "double",
    #               "user_precision_t", "user_precision_t", "user_precision_t",
    #               "user_precision_t"]
    # formats = ["f", "f", "f", "f", "e", "f", "f", "f", "d", "d", "e",
    #            "f", "f", "f", "f"]
    
    arrays = [ref_stokesI, ref_power_SIs,
              ref_freqs, ref_curve_SIs, ref_qs,
              num_list_values, list_start_indexes, list_freqs, list_stokesI,
              stokesV_pol_fracs, ref_stokesV, stokesV_power_SIs, stokesV_qs, stokesV_curve_SIs,
              linpol_pol_fracs, ref_linpol, linpol_power_SIs, linpol_qs, linpol_curve_SIs,
              intr_pol_angle, rms, stokesV_pol_frac_comp_inds, stokesV_power_comp_inds, stokesV_curve_comp_inds,
              linpol_pol_frac_comp_inds, linpol_power_comp_inds, linpol_curve_comp_inds, linpol_angle_inds,
              stokesV_list_ref_flux, stokesV_list_ref_freqs,
              stokesV_num_list_values, stokesV_list_start_indexes, stokesV_list_comp_inds,
              stokesQ_list_ref_flux, stokesQ_list_ref_freqs,
              stokesQ_num_list_values, stokesQ_list_start_indexes, stokesQ_list_comp_inds,
              stokesU_list_ref_flux, stokesU_list_ref_freqs,
              stokesU_num_list_values, stokesU_list_start_indexes, stokesU_list_comp_inds,
              linpol_p_list_ref_flux, linpol_p_list_ref_freqs,
              linpol_p_num_list_values, linpol_p_list_start_indexes, linpol_p_list_comp_inds]
    
    
    
    array_names = ["ref_stokesI",  "ref_power_SIs", "ref_freqs", "ref_curve_SIs", "ref_qs",
                  "num_list_values", "list_start_indexes", "list_freqs", "list_stokesI",
                  "stokesV_pol_fracs", "ref_stokesV", "stokesV_power_SIs", "stokesV_qs", "stokesV_curve_SIs",
                  "linpol_pol_fracs", "ref_linpol", "linpol_power_SIs", "linpol_qs", "linpol_curve_SIs",
                  "intr_pol_angle", "rms", "stokesV_pol_frac_comp_inds", "stokesV_power_comp_inds", "stokesV_curve_comp_inds",
                  "linpol_pol_frac_comp_inds", "linpol_power_comp_inds", "linpol_curve_comp_inds", "linpol_angle_inds",
                  "stokesV_list_ref_flux", "stokesV_list_ref_freqs",
                  "stokesV_num_list_values", "stokesV_list_start_indexes", "stokesV_list_comp_inds",
                  "stokesQ_list_ref_flux", "stokesQ_list_ref_freqs",
                  "stokesQ_num_list_values", "stokesQ_list_start_indexes", "stokesQ_list_comp_inds",
                  "stokesU_list_ref_flux", "stokesU_list_ref_freqs",
                  "stokesU_num_list_values", "stokesU_list_start_indexes", "stokesU_list_comp_inds",
                  "linpol_p_list_ref_flux", "linpol_p_list_ref_freqs",
                  "linpol_p_num_list_values", "linpol_p_list_start_indexes", "linpol_p_list_comp_inds"]
    
    precisions = ["user_precision_t", "user_precision_t", "double", "user_precision_t", "user_precision_t",
                  "int", "int", "double", "user_precision_t",
                  "user_precision_t", "user_precision_t", "user_precision_t", "user_precision_t",
                  "user_precision_t", "user_precision_t", "user_precision_t", "user_precision_t",
                  "user_precision_t", "user_precision_t", "user_precision_t", "user_precision_t",
                  "int", "int","int", "int","int", "int","int",
                  "user_precision_t", "double", "int", "int","int",
                  "user_precision_t", "double", "int", "int","int",
                  "user_precision_t", "double", "int", "int","int",
                  "user_precision_t", "double", "int", "int","int"]

    print_header_file(arrays, array_names, precisions)
