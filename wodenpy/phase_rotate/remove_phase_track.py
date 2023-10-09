import numpy as np

def remove_phase_tracking(frequencies=None, wws_seconds=None,
                          num_time_steps=None, v_container=None,
                          num_baselines=None):
    """
    WARNING - currently does not change the :math:`u,v,w` coordinates, so they
    are still defined via the original phase centre. This function really is
    just to feed uvfits into the RTS (which generates it's own u,v,w using the
    antenna table)

    Undoes phase tracking applied by WODEN - to phase track, a phase was applied
    to counter the delay term caused by :math:`w` term of baseline - so just
    apply the opposite effect of the w term, i.e.

    .. math::
        V^\\prime = V \\exp(2\pi i w)

    where :math:`V` is the phase tracked visibility and :math:`V^\\prime` is
    the visibility after removing phase tracking.

    Parameters
    ----------

    frequencies : float array
        Frequencies of all fine channels (Hz)
    wws_seconds : float array
        The :math:`w` coordinates (seconds)
    num_baselines : int
        Number of baselines
    v_container : float array
        Complex visibility data out of WODEN with phase tracking, with
        `shape=(num_time_steps*num_baselines,1,1,num_freq_channels,4,3))`

    Returns
    -------
    v_container : float array
        Same visibility data as before, with phase tracking returned.

    """

    sign = 1
    PhaseConst = 2j * np.pi * sign

    num_freqs = len(frequencies)

    # print("FREQS", frequencies)

    for time_ind in np.arange(num_time_steps):

        these_wws_secs = wws_seconds[time_ind*num_baselines:(time_ind + 1)*num_baselines]

        for freq_ind, freq in enumerate(frequencies):

            xx_re = v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,0,0]
            xx_im = v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,0,1]

            yy_re = v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,1,0]
            yy_im = v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,1,1]

            xy_re = v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,2,0]
            xy_im = v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,2,1]

            yx_re = v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,3,0]
            yx_im = v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,3,1]

            xx_comp = xx_re + 1j*xx_im
            yy_comp = yy_re + 1j*yy_im
            xy_comp = xy_re + 1j*xy_im
            yx_comp = yx_re + 1j*yx_im

            ##theory - so normal phase delay is caused by path difference across
            ##a base line, which is u*l + v*m + w*n
            ##To phase track, you insert a phase to make sure there is no w contribution at
            ##phase centre; this is when n = 1, so you insert a phase thus:
            ##a base line, which is u*l + v*m + w*(n - 1)
            ##So we just need to remove the effect of the -w term

            wws = these_wws_secs * freq
            phase_rotate = np.exp( PhaseConst * wws)
            xx_comp = xx_comp * phase_rotate
            yy_comp = yy_comp * phase_rotate
            xy_comp = xy_comp * phase_rotate
            yx_comp = yx_comp * phase_rotate

            v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,0,0] = np.real(xx_comp)
            v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,0,1] = np.imag(xx_comp)
            v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,1,0] = np.real(yy_comp)
            v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,1,1] = np.imag(yy_comp)
            v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,2,0] = np.real(xy_comp)
            v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,2,1] = np.imag(xy_comp)
            v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,3,0] = np.real(yx_comp)
            v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,3,1] = np.imag(yx_comp)

    return v_container