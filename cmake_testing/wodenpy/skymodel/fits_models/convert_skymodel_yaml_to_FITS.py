from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, Column
from scipy.optimize import curve_fit
import argparse


D2R = np.pi / 180.0
REF_FREQ = 200.0

class Component_Info():
    """
    Class to contain an RTS source information

    :ivar str comp_type: The component type: either P, G, or S
    :ivar float pa: Position Angle of the component
    :ivar float major: Major angle of the component
    :ivar float minor: Minor angle of the component
    :ivar list shapelet_coeffs: List to contain lists of shapelet coeffs if the source is a S
    """
    def __init__(self):
        self.source_name = None
        self.comp_type = None
        self.pa = np.nan
        self.major = np.nan
        self.minor = np.nan
        self.shapelet_coeffs = []
        self.n1s = []
        self.n2s = []
        self.ra = None
        self.dec = None
        self.flux = None
        self.freq = None

        self.fluxes = []
        self.freqs = []

        self.flux_type = None
        self.curve_q = None
        self.SI = -0.8

        self.norm_comp_pl = np.nan
        self.alpha_pl = np.nan
        self.norm_comp_cpl = np.nan
        self.alpha_cpl = np.nan
        self.curve_cpl = np.nan

    def calc_flux(self, freq):
        """Return the flux of the component at the given frequency
        by scaling via the spectral index"""

        flux = self.flux * (freq/self.freq)**self.SI
        return flux


def read_yaml_srclist(srclist, skip_sources=[]):
    """Reads in sources from a hyperdrive srclist and returns """

    all_freqs = []
    all_names = []

    with open(srclist) as file:

        # current_source = 0

        components = []
        sources = []

        component = False
        source_name = False
        current_source = 0
        comp_count = 0

        source_indexes = []

        freq_count = 0
        freq_indent = 0

        for line in file:
            if line != '---\n' and '#' not in line and line != ''  and line != ' ' and line != '\n':

                if line[0] != ' ' and line[0] != '-':
                    # print(current_source)
                    source_name = line.split('\n')[0].strip(':')
                    all_names.append(source_name)
                    current_source += 1
                    comp_count = -1

                elif 'ra:' in line:

                    ##ra should be the first thing in a component, so we need
                    ##to append all the previously found values and reset the
                    ##counters

                    ##If a previous component exists, append in to the list
                    ##of all components, and then make a new one
                    if component:
                        ##Make some things into np.arrays so we can maths them
                        component.n1s = np.array(component.n1s)
                        component.n2s = np.array(component.n2s)
                        component.shapelet_coeffs = np.array(component.shapelet_coeffs)
                        component.fluxes = np.array(component.fluxes)
                        component.freqs = np.array(component.freqs)
                        component.comp_name = f"{component.source_name}_C{component.comp_count:02d}"

                        if component.source_name in skip_sources:
                            print(f'Skipping component in {component.source_name}')
                        else:
                            components.append(component)


                    component = Component_Info()
                    freq_count = 0

                    component.source_name = source_name
                    comp_count += 1
                    component.comp_count = comp_count
                    
                    ra = float(line.split()[-1])
                    
                    if ra < 0:
                        ra += 360.0
                    
                    component.ra = ra

                    source_indexes.append(current_source)

                elif 'dec:' in line:
                    component.dec = float(line.split()[-1])

                elif 'comp_type: point' in line:
                    component.comp_type = 'P'
                elif 'gaussian:' in line:
                    component.comp_type = 'G'
                elif 'shapelet:' in line:
                    component.comp_type = 'S'

                elif "maj:" in line:
                    component.major = float(line.split()[-1])*(1 / 3600.0)
                elif "min:" in line:
                    component.minor = float(line.split()[-1])*(1 / 3600.0)
                elif "pa:" in line:
                    component.pa = float(line.split()[-1])

                elif 'n1:' in line:
                    component.n1s.append(float(line.split()[-1]))
                elif 'n2:' in line:
                    component.n2s.append(float(line.split()[-1]))
                elif 'value:' in line:
                    component.shapelet_coeffs.append(float(line.split()[-1]))

                elif 'power_law:' in line and 'curved_power_law:' not in line:
                    component.flux_type = 'pl'
                elif 'curved_power_law:' in line:
                    component.flux_type = 'cpl'
                    ##TODO read in curved stuff properly

                elif 'si:' in line:
                    component.SI = float(line.split()[-1])

                elif 'list:' in line:
                    component.flux_type = 'nan'

                elif 'freq:' in line:
                    freq_count += 1
                    component.freqs.append(float(line.split()[-1]))

                    all_freqs.append(float(line.split()[-1]))

                    ##Stick in an empty np.array for Stokes I,Q,U,V
                    component.fluxes.append(np.array([np.nan, np.nan, np.nan,np.nan]))

                    ##See what indent this freq entry starts at - used to
                    ##line up following freq entries, as `q` can either mean
                    ##stokes Q or q curvature param
                    freq_indent = line.index('f')


                elif ' i:' in line:
                    component.fluxes[freq_count - 1][0] = float(line.split()[-1])

                ##Gotta be fancy here to work out if this is a Stokes Q or a
                ##curved power law 'q' param
                elif ' q:' in line:
                    q = float(line.split()[-1])
                    if line.index('q') == freq_indent:
                        component.fluxes[freq_count - 1][1] = q
                    else:
                        if component.flux_type == 'cpl':
                            component.curve_q = q

                elif ' u:' in line:
                    component.fluxes[freq_count - 1][2] = float(line.split()[-1])

                elif ' v:' in line:
                    component.fluxes[freq_count - 1][3] = float(line.split()[-1])


    ##last one doesn't get added to list, so do that
    component.n1s = np.array(component.n1s)
    component.n2s = np.array(component.n2s)
    component.shapelet_coeffs = np.array(component.shapelet_coeffs)
    component.fluxes = np.array(component.fluxes)
    component.freqs = np.array(component.freqs)
    component.comp_name = f"{component.source_name}_C{component.comp_count:02d}"

    if component.source_name in skip_sources:
        print('Skipping component in {component.source_name}')
    else:
        components.append(component)

    all_freqs = np.unique(all_freqs).tolist()

    return components, all_freqs, all_names


def power_law(freqs, SI, ref_flux, ref_freq=REF_FREQ):

    return ref_flux*(freqs / ref_freq)**SI

def fit_power_law(freqs, fluxes):

    # just make sigma 10 percent of the average flux value
    sigma = 0.1*np.mean(fluxes)*np.ones(len(fluxes))


    ##guess everthing has SI = -0.8, is 0.5 Jy at 200Mhz

    param_guess = [-0.8, 0.5]


    popt, pcov = curve_fit(power_law, freqs, fluxes, param_guess, sigma=sigma)

    chi_residish = np.sum(((fluxes - power_law(freqs, *popt)) / sigma)**2)

    return popt, chi_residish


def curved_power_law(freqs, q, SI, ref_flux, ref_freq=REF_FREQ):

    # exp_bit = np.exp(q*np.log(freqs)**2) / np.exp(q*np.log(ref_freq)**2)
    
    si_ratio = (freqs / ref_freq)**SI
    exp_bit = np.exp(q*np.log(freqs / ref_freq)**2)

    return ref_flux*(freqs / ref_freq)**SI*exp_bit

def fit_curved_power_law(freqs, fluxes):

    # just make sigma 10 percent of the average flux value
    sigma = 0.1*np.mean(fluxes)*np.ones(len(fluxes))

    ##guess everthing has SI = -0.8, is 0.5 Jy at 200Mhz
    param_guess = [0, -0.8, 0.5]

    popt, pcov = curve_fit(curved_power_law, freqs, fluxes, param_guess, sigma=sigma)

    chi_residish = np.sum(((fluxes - curved_power_law(freqs, *popt)) / sigma)**2)

    return popt, chi_residish

def fit_component_flux(component, comp_ind=False):

    freqs = component.freqs/1e+6
    ##Only fitting Stokes I here!
    fluxes = component.fluxes
    stokesI = fluxes[:, 0]

    use_indexes = np.where(np.isnan(stokesI) == False)[0]

    ##OK, we can't fit anything if we only have one list entry
    if len(use_indexes) < 2:
        print(f'Only one flux entry for {component.source_name}')
        print("Setting to a power law with SI of -0.8")

        ref_flux = stokesI[use_indexes][0]
        ref_freq = freqs[use_indexes][0]

        component.norm_comp_pl = ref_flux*(REF_FREQ / ref_freq)**-0.8
        component.alpha_pl = -0.8
        component.flux_type = 'pl'

    ##On makes sense to fit a powerlaw with two flux entries
    elif len(use_indexes) == 2:
        # print(f'Only one flux entry for {component.source_name}')
        # print("Setting to a power law with SI of -0.8")

        pl_params, pl_chi = fit_power_law(freqs[use_indexes], stokesI[use_indexes])

        component.norm_comp_pl = power_law(REF_FREQ, *pl_params)
        component.alpha_pl = pl_params[0]
        component.flux_type = 'pl'

    else:

        pl_params, pl_chi = fit_power_law(freqs[use_indexes], stokesI[use_indexes])

        cpl_params, cpl_chi = fit_curved_power_law(freqs[use_indexes], stokesI[use_indexes])

        use_pl = False
        use_cpl = False

        if pl_chi < 2 and pl_chi < cpl_chi:
            use_pl = True

        elif cpl_chi < 2:
            use_cpl = True
        else:
            
            plt.plot(freqs[use_indexes], stokesI[use_indexes], 'o')
            plt.plot(freqs[use_indexes], power_law(freqs[use_indexes], *pl_params), label='PL')
            plt.show()
            ##not fit by models previously, so leave it as a list
            print('Not being fit as reduced chi less than two')
            print('Frequencies (MHz):', freqs[use_indexes])
            print('Stokes I:', stokesI[use_indexes])
            pass

        if use_pl:
            component.flux_type = 'pl'
            component.norm_comp_pl = power_law(REF_FREQ, *pl_params)
            component.alpha_pl = pl_params[0]

        elif use_cpl:
            component.flux_type = 'cpl'
            component.norm_comp_cpl = curved_power_law(REF_FREQ, *cpl_params)
            component.alpha_cpl = cpl_params[1]
            component.curve_cpl = cpl_params[0]

    return component

def calc_pl_norm_at_200MHz(component):

    ref_flux = component.fluxes[0][0]
    ref_freq = component.freqs[0]/1e+6

    component.norm_comp_pl = ref_flux*(REF_FREQ / ref_freq)**component.SI
    component.alpha_pl = component.SI
    component.flux_type = 'pl'

    return component


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='Convert a YAML skymodel to a FITS file')
    parser.add_argument('--input_yaml', type=str, help='Path to the input YAML srclist',
                        required=True)
    parser.add_argument('--output_fits', type=str, help='Name for output FITS srclist',
                        required=True)
    parser.add_argument('--fit_list_entries', default=False, action='store_true',
                        help='Add flag to fit list entries with either a power or curved power-law')
    args = parser.parse_args()

    ##could add an argument to skip sources here
    skip_sources = []
    
    
    components, all_freqs, all_names = read_yaml_srclist(args.input_yaml, skip_sources)

    all_freqs = np.sort(np.unique(all_freqs))
    components = np.array(components)

    num_components = len(components)

    flux_types = np.array([component.flux_type for component in components])

    power_laws = np.where(flux_types == 'pl')[0]
    curve_laws = np.where(flux_types == 'cpl')[0]
    list_laws = np.where(flux_types == 'nan')[0]

    print("Found these types: num power, curved, list", len(power_laws), len(curve_laws), len(list_laws))

    for comp_ind, component in enumerate(components[power_laws]):
        component = calc_pl_norm_at_200MHz(component)
        
    ##If asked for, fit all the list entires with a power-law        
    if args.fit_list_entries:
        for comp_ind, component in enumerate(components[list_laws]):
            component = fit_component_flux(component)

        flux_types = np.array([component.flux_type for component in components])

        power_laws = np.where(flux_types == 'pl')[0]
        curve_laws = np.where(flux_types == 'cpl')[0]
        list_laws = np.where(flux_types == 'nan')[0]

        print("After fitting: num power, curved, list", len(power_laws), len(curve_laws), len(list_laws))

    ##for all components, what SOURCE do they belong to?
    comp_source_names = np.array([component.source_name for component in components])
    comp_names = np.array([component.comp_name for component in components])

    unq_source_ID = Column(data=comp_source_names, name='UNQ_SOURCE_ID')
    name = Column(data=comp_names, name='NAME')


    ras = Column(data=np.array([component.ra for component in components]),
                           name='RA', unit='deg')

    decs = Column(data=np.array([component.dec for component in components]),
                           name='DEC', unit='deg')

    majors = Column(data=np.array([component.major for component in components]),
                           name='MAJOR_DC', unit='deg')

    minors = Column(data=np.array([component.minor for component in components]),
                           name='MINOR_DC', unit='deg')

    pas = Column(data=np.array([component.pa for component in components]),
                           name='PA_DC', unit='deg')

    # main_table.add_columns([unq_source_ID, name])

    mod_type = Column(data=np.array([component.flux_type for component in components]),
                           name='MOD_TYPE', unit='deg')

    norm_comp_pl = Column(data=np.array([component.norm_comp_pl for component in components]),
                          name="NORM_COMP_PL")
    alpha_pl = Column(data=np.array([component.alpha_pl for component in components]),
                          name="ALPHA_PL")
    norm_comp_cpl = Column(data=np.array([component.norm_comp_cpl for component in components]),
                          name="NORM_COMP_CPL")
    alpha_cpl = Column(data=np.array([component.alpha_cpl for component in components]),
                          name="ALPHA_CPL")
    curve_cpl = Column(data=np.array([component.curve_cpl for component in components]),
                          name="CURVE_CPL")
    comp_type = Column(data=np.array([component.comp_type for component in components]),
                          name="COMP_TYPE")


    out_columns = [unq_source_ID, name, ras, decs, majors, minors, pas, mod_type, comp_type,
                   norm_comp_pl, alpha_pl, norm_comp_cpl, alpha_cpl, curve_cpl]
    
    ##gather shapelet type things
    shape_names = []
    shape_n1s = []
    shape_n2s = []
    shape_coeffs = []

    for component in components:
        if component.comp_type == 'S':
            for n1, n2, coeff in zip(component.n1s, component.n2s,
                                     component.shapelet_coeffs):

                shape_names.append(component.comp_name)
                shape_n1s.append(n1)
                shape_n2s.append(n2)
                shape_coeffs.append(coeff)

    s_names = Column(data=np.array(shape_names), name="NAME")
    s_n1s = Column(data=np.array(shape_n1s, dtype=int), name="N1")
    s_n2s = Column(data=np.array(shape_n2s, dtype=int), name="N2")
    s_coeffs = Column(data=np.array(shape_coeffs), name="COEFF")
    
    main_table = Table()
    main_table.add_columns(out_columns)

    shape_table = Table()
    shape_table.add_columns([s_names, s_n1s, s_n2s, s_coeffs])

    # hdu_list = fits.HDUList([
    #     fits.PrimaryHDU(),
    #     fits.table_to_hdu(main_table),
    #     fits.table_to_hdu(shape_table),
    # ])
    
    # hdu_list.writeto('srclist_pumav3_EoR0LoBES_EoR1pietro_CenA-GP_2023-11-07.fits', overwrite=True)
    

    for freq in all_freqs:
        flux_data = np.full(num_components, np.nan)

        for comp_ind, component in enumerate(components):
            # print(component.fluxes.shape)
            fluxes = component.fluxes[np.where(component.freqs == freq)[0]]
            if len(fluxes) == 1:
                stokesI = fluxes[0][0]
                flux_data[comp_ind] = stokesI

        flux_col = Column(data=flux_data, name=f"INT_FLX{freq/1e+6:.3f}", unit='Jy')
        out_columns.append(flux_col)

    main_table = Table()
    main_table.add_columns(out_columns)

    hdu_list = fits.HDUList([
        fits.PrimaryHDU(),
        fits.table_to_hdu(main_table),
        fits.table_to_hdu(shape_table),
    ])

    hdu_list.writeto(args.output_fits, overwrite=True)


