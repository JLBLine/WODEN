import matplotlib.pyplot as plt

##This is terrible practice, so sue me
##A bunch of variables are defined in here
from write_test_extrap_stokes_header import *

def power_law(extrap_freqs, ref_freq, SI, ref_stokesI, ref_stokesQ, ref_stokesU, ref_stokesV):

    flux_ratio = (extrap_freqs / ref_freq)**SI

    extrap_stokesI = ref_stokesI*flux_ratio
    extrap_stokesQ = ref_stokesQ*flux_ratio
    extrap_stokesU = ref_stokesU*flux_ratio
    extrap_stokesV = ref_stokesV*flux_ratio

    return extrap_stokesI, extrap_stokesQ, extrap_stokesU, extrap_stokesV

def curved_power_law(extrap_freqs, ref_freq, SI, q, ref_stokesI, ref_stokesQ, ref_stokesU, ref_stokesV):

    si_ratio = (extrap_freqs / ref_freq)**SI
    
    extrap_freq_use = extrap_freqs
    ref_freq_use = ref_freq

    # exp_ratio = np.exp(q*np.log(extrap_freq_use)**2) / np.exp(q*np.log(ref_freq_use)**2)
    
    exp_bit = np.exp(q*np.log(extrap_freqs / ref_freq)**2)

    extrap_stokesI = ref_stokesI*si_ratio*exp_bit
    extrap_stokesQ = ref_stokesQ*si_ratio*exp_bit
    extrap_stokesU = ref_stokesU*si_ratio*exp_bit
    extrap_stokesV = ref_stokesV*si_ratio*exp_bit

    return extrap_stokesI, extrap_stokesQ, extrap_stokesU, extrap_stokesV


def extrapolate_list_flux(ref_freqs, ref_fluxes, desired_freq) -> float:

    ##Happen to be extrapolating to a reference frequency
    if desired_freq in ref_freqs:
        extrap_flux = ref_fluxes[np.where(ref_freqs == desired_freq)][0]
    else:

        freq_diffs = ref_freqs - desired_freq

        low_ind_1 = -1.0

        low_val_1 = 1e16
        # low_val_2 = 1e16

        for ind in np.arange(len(ref_freqs)):
            abs_diff = abs(freq_diffs[ind])

            if abs_diff < low_val_1:
                low_val_1 = abs_diff
                low_ind_1 = ind

        # print("initial", low_ind_1)

        ##Closest frequency is the lowest
        if low_ind_1 == 0:
            low_ind_2 = 1
        ##Closest frequency is the highest
        elif low_ind_1 == len(ref_freqs) - 1:
            low_ind_2 = low_ind_1 - 1
        ##otherwise, choose either above or below
        else:
            ##closest freq is higher than desired
            if ref_freqs[low_ind_1] > desired_freq:
                low_ind_2 = low_ind_1 - 1
            else:
                low_ind_2 = low_ind_1 + 1

        if ref_fluxes[low_ind_1] <= 0 or ref_fluxes[low_ind_2] <= 0:

            gradient = (ref_fluxes[low_ind_2] - ref_fluxes[low_ind_1]) / (ref_freqs[low_ind_2] - ref_freqs[low_ind_1])
            extrap_flux =  ref_fluxes[low_ind_1] + gradient*(desired_freq - ref_freqs[low_ind_1])

        else:

            flux1 = np.log10(ref_fluxes[low_ind_1])
            flux2 = np.log10(ref_fluxes[low_ind_2])
            freq1 = np.log10(ref_freqs[low_ind_1])
            freq2 = np.log10(ref_freqs[low_ind_2])
            desired_freq = np.log10(desired_freq)

            gradient = (flux2 - flux1) / (freq2 - freq1)
            extrap_flux =  flux1 + gradient*(desired_freq - freq1)

            extrap_flux = 10**extrap_flux


    return extrap_flux

if __name__ == '__main__':
    ##read in the data output by ctest - obviously need to run ctest first!


    ctest_data = np.loadtxt("../../../build/cmake_testing/GPU_or_C_code/source_components/test_extrap_stokes.txt")

    fig, axs = plt.subplots(2, 2, figsize=(8, 8))

    axs = axs.flatten()
    
    num_plots = 4

    for power in range(num_plots):

        expec_I, expec_Q, expec_U, expec_V = power_law(extrap_freqs,
                                          200e+6, stokesI_power_SIs[power],
                                          ref_stokesI[power], ref_linpol[power],
                                          ref_linpol[power], ref_stokesV[power])

        axs[power].loglog(extrap_freqs/1e+6, expec_I, label='Expected I')

        low = power*num_extrap_freqs
        high = low + num_extrap_freqs

        extrap_I = ctest_data[low:high,0]
        extrap_Q = ctest_data[low:high,1]
        extrap_U = ctest_data[low:high,2]
        extrap_V = ctest_data[low:high,3]

        axs[power].loglog(extrap_freqs/1e+6, extrap_I, 'C1o', mfc='none', label='Calculated I')

        axs[power].set_xlabel("Freq (MHz)")
        axs[power].set_ylabel("Flux density (Jy)")
        axs[power].legend()

    fig.tight_layout()
    fig.savefig("test_extrap_power_laws.png", bbox_inches='tight')


    fig, axs = plt.subplots(2, 2, figsize=(8, 8))

    axs = axs.flatten()

    for curve in range(num_plots):


        expec_I, expec_Q, expec_U, expec_V = curved_power_law(extrap_freqs,
                                          200e+6, stokesI_curve_SIs[curve], stokesI_qs[curve],
                                          ref_stokesI[curve], ref_linpol[curve],
                                          ref_linpol[curve], ref_stokesV[curve])

        axs[curve].loglog(extrap_freqs/1e+6, expec_I, label='Expected I')

        low = curve*num_extrap_freqs + num_powers*num_extrap_freqs
        high = low + num_extrap_freqs

        extrap_I = ctest_data[low:high,0]
        extrap_Q = ctest_data[low:high,1]
        extrap_U = ctest_data[low:high,2]
        extrap_V = ctest_data[low:high,3]

        # print(extrap_I)

        axs[curve].loglog(extrap_freqs/1e+6, extrap_I, 'C1o', mfc='none', label='Calculated I')


        axs[curve].set_xlabel("Freq (MHz)")
        axs[curve].set_ylabel("Flux density (Jy)")

        axs[curve].legend()

    fig.tight_layout()
    fig.savefig("test_extrap_curve_power_laws.png", bbox_inches='tight')


    fig, axs = plt.subplots(2, 2, figsize=(8, 8))

    axs = axs.flatten()

    for list_ind, num_list_vals in enumerate(num_list_values[:num_plots]):

        axs[list_ind].set_xscale('log')
        axs[list_ind].set_yscale('symlog')

        start_index = list_start_indexes[list_ind]

        l_freqs = list_freqs[start_index:start_index+num_list_vals]

        l_stokesI = list_stokesI[start_index:start_index+num_list_vals]

        axs[list_ind].plot(l_freqs/1e+6, l_stokesI, '-ok', ms=4, mfc='none',label='List entries')

        python_extrap_I = np.empty(num_extrap_freqs)

        for i in range(num_extrap_freqs):

            python_extrap_I[i] = extrapolate_list_flux(l_freqs, l_stokesI, extrap_freqs[i])

        axs[list_ind].plot(extrap_freqs/1e+6, python_extrap_I, 'C1x', ms=4, mfc='none', label='python extrap')

        low = list_ind*num_extrap_freqs + num_powers*num_extrap_freqs + num_curves*num_extrap_freqs
        high = low + num_extrap_freqs

        extrap_I = ctest_data[low:high,0]
        
        axs[list_ind].plot(extrap_freqs/1e+6, extrap_I, 'cs', ms=6, mfc='none', label='GPU extrap')

        axs[list_ind].set_xlabel("Freq (MHz)")
        axs[list_ind].set_ylabel("Flux density (Jy)")
    axs[0].legend(loc='upper left')

    fig.tight_layout()
    fig.savefig("test_extrap_list_laws.png", bbox_inches='tight')
