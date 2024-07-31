import matplotlib.pyplot as plt

##This is terrible practice, so sue me
##A bunch of variables are defined in here
from write_test_extrap_stokes_header import *

if __name__ == '__main__':
    ##read in the data output by ctest - obviously need to run ctest first!


    ctest_data = np.loadtxt("../../../build/cmake_testing/GPU_code/source_components/test_extrap_stokes.txt")

    

    for comp in range(15):
        
        fig, axs = plt.subplots(2, 2, figsize=(8, 8))
        
        low = comp*num_extrap_freqs
        high = low + num_extrap_freqs
        
        extrap_I = ctest_data[low:high,0]
        extrap_Q = ctest_data[low:high,1]
        extrap_U = ctest_data[low:high,2]
        extrap_V = ctest_data[low:high,3]
        
        expec_I = ctest_data[low:high,4]
        expec_Q = ctest_data[low:high,5]
        expec_U = ctest_data[low:high,6]
        expec_V = ctest_data[low:high,7]

        axs[0,0].loglog(extrap_freqs/1e+6, expec_I, label='CPU Stokes I')
        axs[0,0].loglog(extrap_freqs/1e+6, extrap_I, 'C1o', mfc='none', label='GPU Stokes I')
        
        axs[0,1].plot(extrap_freqs/1e+6, expec_Q, label='CPU Stokes Q')
        axs[0,1].plot(extrap_freqs/1e+6, extrap_Q, 'C1o', mfc='none', label='GPU Stokes Q')
        
        axs[1,0].plot(extrap_freqs/1e+6, expec_U, label='CPU Stokes U')
        axs[1,0].plot(extrap_freqs/1e+6, extrap_U, 'C1o', mfc='none', label='GPU Stokes U')
        
        
        if np.all(expec_V < 0):
            axs[1,1].semilogx(extrap_freqs/1e+6, expec_V, label='CPU Stokes V')
            axs[1,1].semilogx(extrap_freqs/1e+6, extrap_V, 'C1o', mfc='none', label='GPU Stokes V')
        else:
            axs[1,1].loglog(extrap_freqs/1e+6, expec_V, label='CPU Stokes V')
            axs[1,1].loglog(extrap_freqs/1e+6, extrap_V, 'C1o', mfc='none', label='GPU Stokes V')
        
        axs[0,1].set_xscale('log')
        axs[1,1].set_xscale('log')


        for ax in axs.flatten():
            ax.set_xlabel("Freq (MHz)")
            ax.set_ylabel("Flux density (Jy)")
            ax.legend()

        fig.tight_layout()
        fig.savefig(f"eg_fluxes_comp{comp:02d}.png", bbox_inches='tight')


    