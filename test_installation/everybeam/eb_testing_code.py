from subprocess import call
from astropy.io import fits
import numpy as np
from astropy.table import Column, Table
from astropy.time import Time
from astropy.coordinates import SkyCoord
from astropy.coordinates import EarthLocation
from astropy import units as u
import matplotlib.pyplot as plt
import numpy.testing as npt
from astropy.constants import c
from astropy.wcs import WCS
import everybeam as eb
from wodenpy.primary_beam.use_everybeam import load_OSKAR_telescope, load_LOFAR_telescope, get_everybeam_norm, run_everybeam, radec_to_xyz
import erfa
import mwa_hyperbeam

C = c.to('m/s').value

def create_WCS(ra_cent, dec_cent, nside, radec_reso):
    header = fits.Header()

    header['NAXIS'] = 2

    header['NAXIS1'] = nside
    header['NAXIS2'] = nside

    header['CDELT1'] = -radec_reso
    header['CRPIX1'] = int(nside // 2) + 1
    header['CRVAL1'] = ra_cent
    header['CTYPE1'] = 'RA---SIN'
    header['CUNIT1']  = 'deg'

    header['CDELT2'] = radec_reso
    header['CRPIX2'] = int(nside // 2) + 1
    header['CRVAL2'] = dec_cent
    header['CTYPE2'] = 'DEC--SIN'
    header['CUNIT2']  = 'deg'

    header['EQUINOX'] = 2000.0
    
    return header, WCS(header).celestial


def plot_jones_on_sky(all_gx, all_Dx, all_Dy, all_gy, wcs):
    fig, axs = plt.subplots(4, 3, figsize=(12, 12), layout='constrained', subplot_kw={'projection': wcs})

    for row in range(3):
        if row == 2:
            gx = np.abs(all_gx)
            Dx = np.abs(all_Dx)
            Dy = np.abs(all_Dy)
            gy = np.abs(all_gy)
            tag = "Abs"
        
        elif row == 0:
            gx = np.real(all_gx)
            Dx = np.real(all_Dx)
            Dy = np.real(all_Dy)
            gy = np.real(all_gy)
            tag = "Real"
        else:
            gx = np.imag(all_gx)
            Dx = np.imag(all_Dx)
            Dy = np.imag(all_Dy)
            gy = np.imag(all_gy)
            tag = "Imag"

        im0 = axs[0, row].imshow(gx, origin='lower')
        plt.colorbar(im0, ax=axs[0, row])
        axs[0, row].set_title(f'{tag} g_x')
        
        im1 = axs[1, row].imshow(Dx, origin='lower')
        plt.colorbar(im1, ax=axs[1, row])
        axs[1, row].set_title(f'{tag} D_x')
        
        im2 = axs[2, row].imshow(Dy, origin='lower')
        plt.colorbar(im2, ax=axs[2, row])
        axs[2, row].set_title(f'{tag} D_y')
        
        im3 = axs[3, row].imshow(gy, origin='lower')
        plt.colorbar(im3, ax=axs[3, row])
        axs[3, row].set_title(f'{tag} g_y')

    for ax in axs.flatten(): ax.grid()

    plt.show()

def plot_everybeam_on_sky(ra0, dec0, observing_time, freq, station_id, telescope, nside=100, radec_reso=120/100):

    header, wcs = create_WCS(ra0, dec0, nside, radec_reso)
    x_mesh, y_mesh = np.meshgrid(np.arange(nside), np.arange(nside))
    ras, decs = wcs.all_pix2world(x_mesh, y_mesh, 0)

    ras = np.radians(ras.flatten())
    decs = np.radians(decs.flatten())

    all_gx = np.empty(len(ras), dtype=complex)
    all_Dx = np.empty(len(ras), dtype=complex)
    all_Dy = np.empty(len(ras), dtype=complex)
    all_gy = np.empty(len(ras), dtype=complex)

    dir_itrfs = radec_to_xyz(ras, decs, observing_time)
    phase_itrf = radec_to_xyz(np.radians(ra0), np.radians(dec0), observing_time)
    
    norm = get_everybeam_norm(phase_itrf, observing_time, freq, telescope,
                              station_id=station_id)

    ind = 0
    for dir_itrf in dir_itrfs:
        response = run_everybeam(dir_itrf, phase_itrf,
                    observing_time, freq, telescope, station_id=station_id, beam_norms=norm)
        
        all_gx[ind] = response[0,0]
        all_Dx[ind] = response[0,1]
        all_Dy[ind] = response[1,0]
        all_gy[ind] = response[1,1]
        
        ind += 1
        
    all_gx.shape = (nside, nside)
    all_Dx.shape = (nside, nside)
    all_Dy.shape = (nside, nside)
    all_gy.shape = (nside, nside)
    
    plot_jones_on_sky(all_gx, all_Dx, all_Dy, all_gy, wcs)
    
    
def make_sky_models(ra0, dec0):

    for pol in ['I', 'Q', 'U', 'V']:
        
        c_ids = Column(data=np.array([f'{pol}_source']), name='UNQ_SOURCE_ID', dtype='|S20')
        c_names = Column(data=np.array([f'{pol}_source_C00']), name='NAME', dtype='|S20')

        ##Component position
        c_ras = Column(data=np.array([ra0]), name='RA')
        c_decs = Column(data=np.array([dec0]), name='DEC')

        ##This says we have a point source
        c_comp_types = Column(data=np.array(['P'], dtype='|S1'), name="COMP_TYPE", dtype='|S1')
        ##This says we have a Stokes I power-law SED
        c_mod_types = Column(data=np.array(['pl'], dtype='|S3'), name="MOD_TYPE", dtype='|S3')
        
        cols = [c_ids, c_names, c_ras, c_decs, c_comp_types, c_mod_types]
        
        if pol == 'I':
            c_stokes_I_ref = Column(data=np.array([1.0]), name='NORM_COMP_PL')
        if pol == 'Q' or pol == 'U':
            c_stokes_I_ref = Column(data=np.array([0.0]), name='NORM_COMP_PL')
            c_mod_types = Column(data=np.array(['nan'], dtype='|S3'), name="LIN_MOD_TYPE", dtype='|S3')
            
            cols.extend([c_mod_types])
        
            freqs = np.array([150, 250])
            
            if pol == 'Q':
                q_fluxes = np.ones_like(freqs)
                u_fluxes = np.zeros_like(freqs)
            else:
                q_fluxes = np.zeros_like(freqs)
                u_fluxes = np.ones_like(freqs)
                
            q_cols = [c_names]
            u_cols = [c_names]
            
            for freq, q_flux, u_flux in zip(freqs, q_fluxes, u_fluxes):
                q_cols.append(Column(data=np.array([q_flux]), name=f'Q_INT_FLX{freq:.1f}'))
                u_cols.append(Column(data=np.array([u_flux]), name=f'U_INT_FLX{freq:.1f}'))

            
        if pol == 'V':
            c_stokes_I_ref = Column(data=np.array([0.0]), name='NORM_COMP_PL')
            c_stokes_ref = Column(data=np.array([1.0]), name='V_NORM_COMP_PL')
            c_stokes_SI = Column(data=np.array([0.0]), name='V_ALPHA_PL')
            c_mod_types = Column(data=np.array(['pl'], dtype='|S3'), name="V_MOD_TYPE", dtype='|S3')
            cols.extend([c_stokes_ref, c_stokes_SI, c_mod_types])
            
        c_stokes_I_SI = Column(data=np.array([0.0]), name='ALPHA_PL')
        cols.extend([c_stokes_I_ref, c_stokes_I_SI])
        
        main_table = Table(cols)
        
        if pol == 'I' or pol == 'V':
            main_table.write(f'{pol}_source.fits', format='fits', overwrite=True)
        else:
            ##Create an HDU list from the Tables, so we can easily write it to a FITS file
            hdu_list = [fits.PrimaryHDU(),
                        fits.table_to_hdu(main_table),
                        fits.table_to_hdu(Table(q_cols)),
                        fits.table_to_hdu(Table(u_cols)),]

            ##Name the tables so WODEN can understand which is which
            hdu_list[1].name = 'MAIN'
            hdu_list[2].name = 'Q_LIST_FLUXES'
            hdu_list[3].name = 'U_LIST_FLUXES'

            ##convert and write out to a FITS table file
            hdu_list = fits.HDUList(hdu_list)
            hdu_list.writeto(f'{pol}_source.fits', overwrite=True)
            
def read_uvfits(uvfits_name):
    """Real shorthand function to read in visibilities from a WODEN UVFITS file"""
    
    with fits.open(uvfits_name) as hdus:
        data = np.squeeze(hdus[0].data.data)
        
        ##Yup, this is the order of things in a UVFITS file
        XX = data[:, :, 0, 0] + 1j*data[:, :, 0, 1]
        YY = data[:, :, 1, 0] + 1j*data[:, :, 1, 1]
        XY = data[:, :, 2, 0] + 1j*data[:, :, 2, 1]
        YX = data[:, :, 3, 0] + 1j*data[:, :, 3, 1]
            
        
    return XX, XY, YX, YY


def convert_inst_to_stokes(XX, XY, YX, YY):
    """Converts instrumental polarizations to Stokes parameters"""
    
    I = 0.5*(XX + YY)
    Q = 0.5*(XX - YY)
    U = 0.5*(XY + YX)
    V = -0.5j*(XY - YX)
    
    return I, Q, U, V

def test_stokes_recovery(pol, beam, atol = 5e-3):
    
    uvfits_name = f"stokes{pol}_{beam}"
    XX, XY, YX, YY = read_uvfits(f'{uvfits_name}_band01.uvfits')

    ##pick a random baseline to plot, they should all be the same
    baseline = np.random.randint(0, XX.shape[0])

    recover_I, recover_Q, recover_U, recover_V = convert_inst_to_stokes(XX[baseline], XY[baseline], YX[baseline], YY[baseline])
    
    ##There are beam leakages so these shouldn't be perfect, but they
    ##should be close enough as we're near beam centre and zenith
    
    if pol == 'I':
        print('Testing Stokes I')
        npt.assert_allclose(recover_I, 1.0 + 0j, atol=atol)
        npt.assert_allclose(recover_Q, 0.0 + 0j, atol=atol)
        npt.assert_allclose(recover_U, 0.0 + 0j, atol=atol)
        npt.assert_allclose(recover_V, 0.0 + 0j, atol=atol)
        print('Stokes I passed')
        
    elif pol == 'Q':
        print('Testing Stokes Q')
        npt.assert_allclose(recover_I, 0.0 + 0j, atol=atol)
        npt.assert_allclose(recover_Q, 1.0 + 0j, atol=atol)
        npt.assert_allclose(recover_U, 0.0 + 0j, atol=atol)
        npt.assert_allclose(recover_V, 0.0 + 0j, atol=atol)
        print('Stokes Q passed')
        
    elif pol == 'U':
        print('Testing Stokes U')
        npt.assert_allclose(recover_I, 0.0 + 0j, atol=atol)
        npt.assert_allclose(recover_Q, 0.0 + 0j, atol=atol)
        npt.assert_allclose(recover_U, 1.0 + 0j, atol=atol)
        npt.assert_allclose(recover_V, 0.0 + 0j, atol=atol)
        print('Stokes U passed')
        
    elif pol == 'V':
        print('Testing Stokes V')
        npt.assert_allclose(recover_I, 0.0 + 0j, atol=atol)
        npt.assert_allclose(recover_Q, 0.0 + 0j, atol=atol)
        npt.assert_allclose(recover_U, 0.0 + 0j, atol=atol)
        npt.assert_allclose(recover_V, 1.0 + 0j, atol=atol)
        print('Stokes V passed')


def getFDF(dataQ, dataU, freqs, startPhi, stopPhi, dPhi, dType='float32'):
    """
    # Perform RM-synthesis on Stokes Q and U data
    #
    # dataQ, dataU and freqs - contains the Q/U data at each frequency (in Hz) measured.
    # startPhi, dPhi - the starting RM (rad/m^2) and the step size (rad/m^2)
    
    Author: Emil Lenc
    """

    # Calculate the RM sampling
    phiArr = np.arange(startPhi, stopPhi, dPhi)

    # Calculate the frequency and lambda sampling
    lamSqArr = np.power(C / np.array(freqs), 2.0)

    # Calculate the dimensions of the output RM cube
    nPhi = len(phiArr)

    # Initialise the complex Faraday Dispersion Function (FDF)
    FDF = np.ndarray((nPhi), dtype='complex')

    # Assume uniform weighting
    wtArr = np.ones(len(lamSqArr), dtype=dType)

    K = 1.0 / np.nansum(wtArr)

    # Get the weighted mean of the LambdaSq distribution (B&dB Eqn. 32)
    lam0Sq = K * np.nansum(lamSqArr)

    # Mininize the number of inner-loop operations by calculating the
    # argument of the EXP term in B&dB Eqns. (25) and (36) for the FDF
    a = (-2.0 * 1.0j * phiArr)
    b = (lamSqArr - lam0Sq) 
    arg = np.exp( np.outer(a, b) )

    # Create a weighted complex polarised surface-brightness cube
    # i.e., observed polarised surface brightness, B&dB Eqns. (8) and (14)
    Pobs = (np.array(dataQ) + 1.0j * np.array(dataU))

    # Calculate the Faraday Dispersion Function
    # B&dB Eqns. (25) and (36)
    FDF = K * np.nansum(Pobs * arg, 1)
    return FDF, phiArr

def findpeaks(freqs, fdf, phi, rmsf, rmsfphi, nsigma):
    """Find peaks in the FDF
    Author: Emil Lenc"""
    # Create the Gaussian filter for reconstruction
    lam2 = (C / freqs) ** 2.0
    lam02 = np.mean(lam2)
    minl2 = np.min(lam2)
    maxl2 = np.max(lam2)
    width = (2.0 * np.sqrt(3.0)) / (maxl2 - minl2)

    Gauss = np.exp((-rmsfphi ** 2.0) / (2.0 * ((width / 2.355) ** 2.0)))
    components = np.zeros((len(phi)), np.float32)
    peaks = []
    phis = []
    std = 0.0
    rmsflen = int((len(rmsf) - 1) / 2)
    fdflen = len(phi) + rmsflen
    while True:
        std = np.std(np.abs(fdf))
        peak1 = np.max(np.abs(fdf))
        pos1 = np.argmax(np.abs(fdf))
        val1 = phi[pos1]
        if peak1 < nsigma * std :
            break
        fdf -= rmsf[rmsflen - pos1:fdflen - pos1] * fdf[pos1]
        peaks.append(peak1)
        phis.append(val1)
        components[pos1] += peak1
    fdf += np.convolve(components, Gauss, mode='valid')
    return phis, peaks, std


def test_RM_recovery(uvfits_name, expected_RM, expected_pol_frac, freqs, atol=5e-3):


    baseline = 0
    XX, XY, YX, YY = read_uvfits(f'{uvfits_name}_band01.uvfits')
    recover_I, recover_Q, recover_U, recover_V = convert_inst_to_stokes(XX[baseline], XY[baseline], YX[baseline], YY[baseline])

    startPhi = -100.0
    dPhi     = 1.0
    stopPhi = -startPhi+dPhi
    lambda2 = np.power(C / np.array(freqs), 2.0)

    df = []
    for f in range(1, len(freqs)):
        df.append(freqs[f] - freqs[f-1])
    chanBW = np.min(np.array(df))
    fmin = np.min(freqs)
    fmax = np.max(freqs)
    bw = fmax - fmin

    dlambda2 = np.power(C / fmin, 2) - np.power(C / (fmin + chanBW), 2)
    Dlambda2 = np.power(C / fmin, 2) - np.power(C / (fmin + bw), 2)
    phimax = np.sqrt(3) / dlambda2
    dphi = 2.0 * np.sqrt(3) / Dlambda2
    phiR = dphi / 5.0
    Nphi = 2 * phimax / phiR
    # ##These are things that mean something to polarisation peoples
    # print("Input resolutions going into RM-synthesis:")
    # print("\tFrequency range: %7.3f MHz - %7.3f MHz" %(fmin / 1.0e6, (fmin + bw) / 1.0e6))
    # print("\tBandwidth: %7.3f MHz" %(bw / 1.0e6))
    # print("\tChannel width: %.1f KHz" %(chanBW / 1.0e3))
    # print("\tdlambda2: %7.3f" %(dlambda2))
    # print("\tDlambda2: %7.3f" %(Dlambda2))
    # print("\tphimax: %7.3f" %(phimax))
    # print("\tdphi: %7.3f" %(dphi))
    # print("\tphiR: %7.3f" %(phiR))
    # print("\tNphi: %7.3f" %(Nphi))
    fwhm = dphi

    # Determine the FDF using the q, u and ferquency values read from the file.
    dirty, phi = getFDF(recover_Q, recover_U, freqs, startPhi, stopPhi, dPhi)
    FDFqu, phi = getFDF(recover_Q, recover_U, freqs, startPhi, stopPhi, dPhi)
    rstartPhi = startPhi * 2
    rstopPhi = stopPhi * 2 - dPhi
    RMSF, rmsfphi = getFDF(np.ones((len(recover_Q))), np.zeros((len(recover_Q))), freqs, rstartPhi, rstopPhi, dPhi)

    phis, peaks, sigma = findpeaks(np.array(freqs), FDFqu, phi, RMSF, rmsfphi, 6.0)
    snr = peaks / sigma
    phierr = fwhm / (2 * snr)
    
    fig, axs = plt.subplots(1, 1, figsize=(6,4))

    axs.plot(phi, np.real(FDFqu), 'k-', label='FDF')

    axs.set_xlabel("$\phi$ (rad m$^{-2}$)")
    axs.set_ylabel("Polarised Flux (Jy RMSF$^{-1}$)")

    axs.set_xlim(-100, 100)

    axs.axvline(x=expected_RM, color='C1', linestyle='--', label='Expected RM')
    axs.legend()

    plt.show()
    
    polarised_flux = np.abs(recover_Q + 1j*recover_U)
    recovered_pol_frac = (polarised_flux.real / recover_I.real)
    
    print('Recovered RM:', phis[0], 'Expected RM:', expected_RM)
    print('Recovered Pol. Fraction:', np.mean(recovered_pol_frac), 'Expected Pol Fraction', expected_pol_frac)
    
    npt.assert_allclose(recovered_pol_frac, expected_pol_frac, atol=atol)
    npt.assert_allclose(phis[0], expected_RM, atol=atol)
    
    
def make_RM_skymodel(ra0, dec0):    
    ##These are the parameters for the source
    alpha = -0.8
    ##hyperdrive/WODEN/LoBES catalogues are all referenced to 200MHz, so we need
    ##to scale the stokesI flux
    stokesI = 10
    pol_frac = 0.5
    phi_RM = 50
    chi_0 = np.radians(60)


    ##The LoBES/hyperdrive/WODEN catalogues need a unique source ID (UNQ_SOURCE_ID),
    ##a then edge component of that source needs a name (NAME)
    c_ids = Column(data=np.array(['RM_source']), name='UNQ_SOURCE_ID', dtype='|S20')
    c_names = Column(data=np.array(['RM_source_C00']), name='NAME', dtype='|S20')

    ##Component position
    c_ras = Column(data=np.array([ra0]), name='RA')
    c_decs = Column(data=np.array([dec0]), name='DEC')

    ##This says we have a point source
    c_comp_types = Column(data=np.array(['P'], dtype='|S1'), name="COMP_TYPE", dtype='|S1')
    ##This says we have a Stokes I power-law SED
    c_mod_types = Column(data=np.array(['pl'], dtype='|S3'), name="MOD_TYPE", dtype='|S3')
    ##Stokes I power law parameters
    c_stokes_I_ref = Column(data=np.array([stokesI]), name='NORM_COMP_PL')
    c_stokes_SI = Column(data=np.array([alpha]), name='ALPHA_PL')

    ##This says we are using a polarisation fraction linear polarisation model
    c_lin_mod_type = Column(data=np.array(['pf'], dtype='|S3'), name='LIN_MOD_TYPE', dtype='|S3')
    ##linear polarisation parameters
    c_lin_pol_frac = Column(data=np.array([pol_frac]), name='LIN_POL_FRAC')
    c_lin_pol_angle = Column(data=np.array([chi_0]), name='INTR_POL_ANGLE')
    c_rm = Column(data=np.array([phi_RM]), name='RM')

    columns = [c_ids, c_names, c_ras, c_decs, c_comp_types, c_mod_types, c_stokes_I_ref, c_stokes_SI, c_lin_mod_type, c_lin_pol_frac, c_lin_pol_angle, c_rm]

    main_table = Table(columns)
    main_table.write('RM_source.fits', format='fits', overwrite=True)
    
    return phi_RM, pol_frac