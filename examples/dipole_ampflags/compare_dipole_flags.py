from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

def read_uvfits(file_name):
    """Read in the E-W and N-S visibilities from a uvfits file."""
    with fits.open(file_name) as hdulist:
        data = np.squeeze(hdulist[0].data.data)
        
        xx = data[:, :, 0, 0] + 1j * data[:, :, 0, 1]
        yy = data[:, :, 1, 0] + 1j * data[:, :, 1, 1]
        
    return xx.flatten(), yy.flatten()

def read_dipamps(file_name):
    with fits.open(file_name) as hdulist:
        dipamps = hdulist[1].data['Dipamps']
        
    return dipamps



if __name__ == "__main__":
    
    
    
    wod_xx_default, wod_yy_default = read_uvfits("woden_noflags_band01.uvfits")
    wod_xx_deformed, wod_yy_deformed = read_uvfits("woden_dipflags_band01.uvfits")
    
    hyp_xx_default, hyp_yy_default = read_uvfits("hyperdrive_noflags.uvfits")
    hyp_xx_deformed, hyp_yy_deformed = read_uvfits("hyperdrive_dipflags.uvfits")
    
    
    atol=5e-6
    rtol=1e-5
    
    allclose_perf = np.allclose(np.real(wod_xx_default), np.real(hyp_xx_default), atol=atol, rtol=rtol)
    print("Perfect beams real E-W all close? ", allclose_perf)
    allclose_perf = np.allclose(np.imag(wod_xx_default), np.imag(hyp_xx_default), atol=atol, rtol=rtol)
    print("Perfect beams imag E-W all close? ", allclose_perf)
    
    allclose_perf = np.allclose(np.real(wod_yy_default), np.real(hyp_yy_default), atol=atol, rtol=rtol)
    print("Perfect beams real N-S all close? ", allclose_perf)
    allclose_perf = np.allclose(np.imag(wod_yy_default), np.imag(hyp_yy_default), atol=atol, rtol=rtol)
    print("Perfect beams imag N-S all close? ", allclose_perf)
    
    allclose_deform = np.allclose(np.real(wod_xx_deformed), np.real(hyp_xx_deformed), atol=atol, rtol=rtol)
    print("Deformed beams real E-W all close? ", allclose_deform)
    allclose_deform = np.allclose(np.imag(wod_xx_deformed), np.imag(hyp_xx_deformed), atol=atol, rtol=rtol)
    print("Deformed beams imag E-W all close? ", allclose_deform)
    
    allclose_deform = np.allclose(np.real(wod_yy_deformed), np.real(hyp_yy_deformed), atol=atol, rtol=rtol)
    print("Deformed beams real N-S all close? ", allclose_deform)
    allclose_deform = np.allclose(np.imag(wod_yy_deformed), np.imag(hyp_yy_deformed), atol=atol, rtol=rtol)
    print("Deformed beams imag N-S all close? ", allclose_deform)
    
    fig, axs = plt.subplots(2, 2, figsize=(8,8))
    
    # bins = np.linspace(0.9988, 0.99895, 10)
    
    axs[0,0].hist(np.abs(wod_xx_default), histtype='step', 
                  label='E-W WODEN perfect', lw=2, hatch='/')#, bins=bins)
    axs[0,0].hist(np.abs(hyp_xx_default), histtype='step',
                  label='E-W hyperdrive perfect', lw=2)#, bins=bins)
    
    axs[0,1].hist(np.abs(wod_yy_default), histtype='step',
                  label='N-S WODEN perfect', lw=2, hatch='/')#, bins=bins)
    axs[0,1].hist(np.abs(hyp_yy_default), histtype='step',
                  label='N-S hyperdrive perfect', lw=2)#, bins=bins)
    
    bins = np.linspace(0.0, 1.0, 100)
    bins = False
    
    bins = np.linspace(0.65, 1.0, 50)
    
    if type(bins) == np.ndarray:
    
        axs[1,0].hist(np.abs(wod_xx_deformed), histtype='step',
                    label='E-W WODEN deformed', lw=2, bins=bins, hatch='/')
        axs[1,0].hist(np.abs(hyp_xx_deformed), histtype='step',
                    label='E-W hyperdrive deformed', lw=2, bins=bins)
        axs[1,1].hist(np.abs(wod_yy_deformed), histtype='step',
                    label='N-S WODEN deformed', lw=2, bins=bins, hatch='/')
        axs[1,1].hist(np.abs(hyp_yy_deformed), histtype='step',
                    label='N-S hyperdrive deformed', lw=2, bins=bins)
        
    else:
    
        axs[1,0].hist(np.abs(wod_xx_deformed), histtype='step',
                    label='E-W WODEN flagged', lw=2, hatch='/')
        axs[1,0].hist(np.abs(hyp_xx_deformed), histtype='step',
                    label='E-W hyperdrive flagged', lw=2)
        
        axs[1,1].hist(np.abs(wod_yy_deformed), histtype='step',
                    label='N-S WODEN flagged', lw=2, hatch='/')
        axs[1,1].hist(np.abs(hyp_yy_deformed), histtype='step',
                    label='N-S hyperdrive flagged', lw=2)
    
    
    axs[0,0].set_ylabel('Count')
    axs[1,0].set_ylabel('Count')
    
    axs[1,0].set_xlabel('Visi amp (Jy)')
    axs[1,1].set_xlabel('Visi amp (Jy)')
    
    
    for ax in axs.flatten():
        ax.legend()
    
    
    plt.tight_layout()
    # fig.savefig('visi_amp_hists.png', bbox_inches='tight')
    fig.savefig('compare_woden_hyp_dipflags.png', bbox_inches='tight')
    plt.close()

    