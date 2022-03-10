import matplotlib.pyplot as plt
import numpy as np
from subprocess import call
from shamfi.shamfi_plotting import add_colourbar


def load_data(filename):
    azs, zas, gx_re, gx_im, Dx_re, Dx_im, Dy_re, Dy_im, gy_re, gy_im, freqs = np.loadtxt(filename, unpack=True)

    gx = gx_re + gx_im*1j
    Dx = Dx_re + Dx_im*1j
    Dy = Dy_re + Dy_im*1j
    gy = gy_re + gy_im*1j

    return azs, zas, gx, Dx, Dy, gy, freqs


def plot_jones(azs, zas, gx, Dx, Dy, gy, freqs, filename):
    outname = filename.split('/')[-1].split('_')[-1]

    num_comps = int(outname.strip('az.txt'))
    print(num_comps)

    num_freqs = 2
    num_times = 2

    for time in range(num_times):
        for freq_ind, freq in enumerate([100e+6, 200e+6]):
            low = time*num_freqs*num_comps + freq_ind*num_comps
            high = low + num_comps

            fig = plt.figure(figsize=(12, 12))

            ax1 = fig.add_subplot(2,2,1, projection='polar')
            ax2 = fig.add_subplot(2,2,2, projection='polar')
            ax3 = fig.add_subplot(2,2,3, projection='polar')
            ax4 = fig.add_subplot(2,2,4, projection='polar')
            # ax5 = fig.add_subplot(4,2,5, projection='polar')
            # ax6 = fig.add_subplot(4,2,6, projection='polar')
            # ax7 = fig.add_subplot(4,2,7, projection='polar')
            # ax8 = fig.add_subplot(4,2,8, projection='polar')

            ax1.scatter(azs[low:high], zas[low:high], c=np.real(gx[low:high]), s=3, alpha=0.6)
            # ax2.scatter(azs[low:high], zas[low:high], c=np.imag(gx[low:high]), s=3, alpha=0.6)
            ax2.scatter(azs[low:high], zas[low:high], c=np.real(Dx[low:high]), s=3, alpha=0.6)
            # ax4.scatter(azs[low:high], zas[low:high], c=np.imag(Dx[low:high]), s=3, alpha=0.6)
            ax3.scatter(azs[low:high], zas[low:high], c=np.real(Dy[low:high]), s=3, alpha=0.6)
            # ax6.scatter(azs[low:high], zas[low:high], c=np.imag(Dy[low:high]), s=3, alpha=0.6)
            ax4.scatter(azs[low:high], zas[low:high], c=np.real(gy[low:high]), s=3, alpha=0.6)
            # ax8.scatter(azs[low:high], zas[low:high], c=np.imag(gy[low:high]), s=3, alpha=0.6)

            ax1.set_title('Real gx')
            # ax2.set_title('Imag gx')
            ax2.set_title('Real Dx')
            # ax4.set_title('Imag Dx')
            ax3.set_title('Real Dy')
            # ax6.set_title('Imag Dy')
            ax4.set_title('Real gy')
            # ax8.set_title('Imag gy')

            for ax in [ax1, ax2, ax3, ax4]:
                ax.set_yticks([])

            fig.savefig('plots/jones_MWA_analy_gains_azza{:d}_time{:d}_freq{:.3f}.png'.format(num_comps, time, freq/1e+6),bbox_inches='tight')
            plt.close()


            fig = plt.figure(figsize=(12, 12))

            ax1 = fig.add_subplot(2,2,1, projection='polar')
            ax2 = fig.add_subplot(2,2,2, projection='polar')
            ax3 = fig.add_subplot(2,2,3, projection='polar')
            ax4 = fig.add_subplot(2,2,4, projection='polar')

            np.real(gx[low:high])
            np.real(Dx[low:high])
            np.real(Dy[low:high])
            np.real(gy[low:high])

            xx = np.real(gx[low:high])**2 + np.real(Dx[low:high])**2
            xy = np.real(gx[low:high])*np.real(Dy[low:high]) + np.real(Dx[low:high])*np.real(gy[low:high])
            yx = np.real(Dy[low:high])*np.real(gx[low:high]) + np.real(gy[low:high])*np.real(Dx[low:high])

            # XY = (g1xx*g2yx_conj + g1xy*g2yy_conj)*sI
            # YX = (g1yx*g2xx_conj + g1yy*g2xy_conj)*sI

            yy = np.real(gy[low:high])**2 + np.real(Dy[low:high])**2

            ax1.scatter(azs[low:high], zas[low:high], c=xx, s=1, alpha=0.6, vmin=0.0)
            ax2.scatter(azs[low:high], zas[low:high], c=xy, s=1, alpha=0.6)
            ax3.scatter(azs[low:high], zas[low:high], c=yx, s=1, alpha=0.6)
            ax4.scatter(azs[low:high], zas[low:high], c=yy, s=1, alpha=0.6)
            # ax8.scatter(azs[low:high], zas[low:high], c=np.imag(gy[low:high]), s=3, alpha=0.6)

            ax1.set_title('Real XX')
            # ax2.set_title('Imag gx')
            ax2.set_title('Real XY')
            # ax4.set_title('Imag Dx')
            ax3.set_title('Real YX')
            # ax6.set_title('Imag Dy')
            ax4.set_title('Real YY')
            # ax8.set_title('Imag gy')

            for ax in [ax1, ax2, ax3, ax4]:
                ax.set_yticks([])

            fig.savefig('plots/inst_pol_MWA_analy_gains_azza{:d}_time{:d}_freq{:.3f}.png'.format(num_comps, time, freq/1e+6),bbox_inches='tight')
            plt.close()


def plot_jones_square(azs, zas, gx, Dx, Dy, gy, freqs, filename, tag='MWA_analy',
                      num_freqs=1, num_times=1):
    outname = filename.split('/')[-1].split('_')[-1]

    num_comps = int(outname.strip('az.txt'))
    # print(num_comps)

    nside = int(np.sqrt(num_comps))

    freqs = [150e+6, 200e+6]

    for time in range(num_times):
        for freq_ind in range(num_freqs):
            freq = freqs[freq_ind]
            low = time*num_freqs*num_comps + freq_ind*num_comps
            high = low + num_comps

            # print(low,high)

            fig, axs = plt.subplots(2, 2, figsize=(12, 12))

            this_gx = gx[low:high]
            this_Dx = Dx[low:high]
            this_Dy = Dy[low:high]
            this_gy = gy[low:high]

            this_gx.shape = (nside, nside)
            this_Dx.shape = (nside, nside)
            this_Dy.shape = (nside, nside)
            this_gy.shape = (nside, nside)

            im1 = axs[0,0].imshow(this_gx.real, origin='lower')
            im2 = axs[0,1].imshow(this_Dx.real, origin='lower')
            im3 = axs[1,0].imshow(this_Dy.real, origin='lower')
            im4 = axs[1,1].imshow(this_gy.real, origin='lower')

            ims = [im1, im2, im3, im4]

            for im, ax in zip(ims, axs.flatten()):
                add_colourbar(im=im, ax=ax, fig=fig)

            axs[0,0].set_title('Real gx')
            axs[0,1].set_title('Real Dx')
            axs[1,0].set_title('Real Dy')
            axs[1,1].set_title('Real gy')
            #
            for ax in axs.flatten():
                ax.set_yticks([])
                ax.set_xticks([])

            fig.savefig(f'plots/jones_{tag}_gains_nside{nside:d}_t{time:02d}_f{freq/1e+6:.3f}MHz.png',bbox_inches='tight')
            plt.close()


            fig, axs = plt.subplots(2, 2, figsize=(12, 12))

            # xx = this_gx**2 + this_Dx**2
            # xy = this_gx*this_Dy + this_Dx*this_gy
            # yx = this_Dy*this_gx + this_gy*this_Dx
            # yy = this_gy**2 + this_Dy**2


            xx = this_gx*np.conjugate(this_gx) + this_Dx*np.conjugate(this_Dx)
            xy = this_gx*np.conjugate(this_Dy) + this_Dx*np.conjugate(this_gy)
            yx = this_Dy*np.conjugate(this_gx) + this_gy*np.conjugate(this_Dx)
            yy = this_gy*np.conjugate(this_gy) + this_Dy*np.conjugate(this_Dy)

            print(np.max(np.abs(xx)), np.max(np.abs(yy)))

            # im1 = axs[0,0].imshow(np.log10(np.real(xx)), origin='lower') #
            # im2 = axs[0,1].imshow(np.log10(np.real(xy)), origin='lower') #
            # im3 = axs[1,0].imshow(np.log10(np.real(yx)), origin='lower') #
            # im4 = axs[1,1].imshow(np.log10(np.real(yy)), origin='lower') #

            im1 = axs[0,0].imshow(np.real(xx), origin='lower',vmin=0,vmax=0.2)
            im2 = axs[0,1].imshow(np.real(xy), origin='lower')# ,vmin=0,vmax=0.3)
            im3 = axs[1,0].imshow(np.real(yx), origin='lower')# ,vmin=0,vmax=0.3)
            im4 = axs[1,1].imshow(np.real(yy), origin='lower',vmin=0,vmax=0.2)

            ims = [im1, im2, im3, im4]

            for im, ax in zip(ims, axs.flatten()):
                add_colourbar(im=im, ax=ax, fig=fig)

            axs[0,0].set_title('XX')
            axs[0,1].set_title('XY')
            axs[1,0].set_title('YX')
            axs[1,1].set_title('YY')
            #
            for ax in axs.flatten():
                ax.set_yticks([])
                ax.set_xticks([])

            fig.savefig(f'plots/linear_pol_{tag}_gains_nside{nside:d}_t{time:02d}_f{freq/1e+6:.3f}MHz.png',bbox_inches='tight')
            plt.close()

            # fig = plt.figure(figsize=(8,4))
            #
            # ax0 = fig.add_subplot(1,2,1)
            # ax1 = fig.add_subplot(1,2,2)
            # # ax2 = fig.add_subplot(2,1,2)
            #
            # im0 = ax0.imshow(np.log10(np.real(xx)), origin='lower',cmap='gnuplot', vmin=-7, vmax=0)
            # im1 = ax1.imshow(np.log10(np.real(yy)), origin='lower',cmap='gnuplot', vmin=-7, vmax=0)
            # # im2 = ax2.imshow(np.real(xx) - np.real(yy), origin='lower',cmap='gnuplot')
            #
            # add_colourbar(ax=ax0,im=im0, fig=fig, label = '$\log_{10}$(Gain)')
            # add_colourbar(ax=ax1,im=im1, fig=fig, label = '$\log_{10}$(Gain)')
            #
            # # add_colourbar(ax=ax2,im=im2, fig=fig, label = 'Gain difference')
            #
            # ax0.set_title('XX')
            # ax1.set_title('YY')
            # # ax2.set_title('XX - YY')
            # #
            # for ax in [ax0,ax1]:
            #     ax.set_yticks([])
            #     ax.set_xticks([])
            #
            # plt.tight_layout()
            # fig.savefig(filename.split('/')[-1] + '.png', bbox_inches='tight')
            # plt.close()


            # plt.tight_layout()
            # fig.savefig('for_talk_WODEN_GPU_beam_nside{:d}_time{:d}_freq{:.3f}.png'.format(nside, time, freq/1e+6),bbox_inches='tight')
            # plt.close()

if __name__ == '__main__':

    call('mkdir -p plots', shell=True)

    filename = '../../build/cmake_testing/primary_beam_cuda/MWA_analy_gains_azza40401.txt'
    azs, zas, gx, Dx, Dy, gy, freqs = load_data(filename)
    plot_jones_square(azs, zas, gx, Dx, Dy, gy, freqs, filename,
                      num_times=2, num_freqs=2)
    #
    filename = '../../build/cmake_testing/primary_beam_cuda/MWA_FEE_gains_azza40401.txt'
    azs, zas, gx, Dx, Dy, gy, freqs = load_data(filename)
    plot_jones_square(azs, zas, gx, Dx, Dy, gy, freqs, filename, tag='MWA_FEE')

    # filename = 'RTS_analytic_beam_gains_azza40401.txt'
    # azs, zas, gx, Dx, Dy, gy, freqs = load_data(filename)
    # plot_jones_square(azs, zas, gx, Dx, Dy, gy, freqs, filename, tag='RTS_analy')
