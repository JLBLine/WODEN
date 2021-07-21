---
title: '`WODEN`: A CUDA-enabled package to simulate low-frequency radio interferometric data'
tags:
  - Python
  - C
  - CUDA
  - radio astronomy
  - interferometers
authors:
  - name: Jack L. B. Line
    orcid: 0000-0002-9130-5920
    affiliation: "1, 2" # (Multiple affiliations must be quoted)
affiliations:
 - name: International Centre for Radio Astronomy Research, Curtin University, Perth, WA 6845, Australia
   index: 1
 - name: ARC Centre of Excellence for All Sky Astrophysics in 3 Dimensions (ASTRO 3D)
   index: 2
date: XXX
bibliography: paper.bib

---
# Summary

`WODEN` is designed to simulate the response of a class of telescope known as an interferometer, producing output "visibilities" for a given astrophysical sky model. Simulated observations allow us to test other software packages that are designed to calibrate and analyse real interferometric data, including verifying expected behaviour with known inputs, and trialling new sky modelling techniques. The `WODEN` sky model can be specified in dirac-delta like functions on the sky (known in the field as "point sources") elliptical Gaussian models, or built out of "shapelet" basis functions (allowing complicated morphologies to be created). Users are able to input a bespoke layout for the interferometer, vary a number of observational parameters including time of day, length of observation and frequency coverage, and select from a number of predefined primary beams which encode the response of the receiving elements of an interferometer. This allows simulations of a number of telescopes to be undertaken. `WODEN` works with all Stokes $I,Q,U,V$ sky polarisations, simulating telescopes with dual linear polarisations.

The core functionality of `WODEN` is written in CUDA as interferometric simulations are computationally intensive but embarrassingly parallel. The compute performance of CUDA allows for large-scale simulations to be run including emission from all directions in the sky. This is paramount for interferometers with a widefield of view such as the Murchison Widefield Array (MWA, @Tingay2013). A Python wrapper is used to take advantage of community packages such as astropy [@astropy2013, @astropy2018], and to present a user-friendly interface to `WODEN`. Those simulating MWA observations can use the MWA `metafits` file to quickly feed in observational parameters to `WODEN` to match real data.

# Statement of need
An interferometer creates visibilities $V$ by cross-correlating signals detected between pairs of antennas or dishes (baselines), described by coordinates $u,v,w$. Each visibility is sensitive to the entire sky, directions of which we describe by the direction cosines $l,m,n$. The full integral can be discretised as

$$
V(u_i,v_i,w_i) = \\ \sum_j \mathcal{B}(l_j,m_j)I(l_j,m_j) \exp[-2\pi i(u_il_j + v_im_j + w_i(n_j-1))],
$$

where $u_i,v_i,w_i$ are the visibility coordinates of the $i^{\mathrm{th}}$ baseline, $l_j$, $m_j$, $n_j$ is the sky position of the $j^{\mathrm{th}}$ component in the sky model, $I(l_j,m_j)$ is the flux density of that component, and $\mathcal{B}(l_j,m_j)$ the instrument beam pattern.

For a telescope like the MWA, the primary beam $\mathcal{B}(l,m)$ is a complicated pattern on the sky, which is sensitive to emission from directly overhead to all the way down to the horizon. To truly capture the effects of astrophysical foregrounds we therefore have to simulate the entire sky. The MWA Fully Embedded Element (FEE, @Sokolowski2017) model is currently the most accurate representation of the MWA primary beam, and is incorporated into `WODEN`.

Under a formalism like the above, that splits the sky into discrete points, pushing $j>=25\times10^6$ can be required to achieve the angular resolution required. Furthermore, $u,v,w$ are time and frequency dependent, so to sample in frequency of order 500 times and 100 samples in time, there are of order $10^{12}$ visibility calculations to make. This makes CUDA acceleration paramount.

Alternative approaches to interferometric simulations exist, such as [pyuvsim](https://github.com/RadioAstronomySoftwareGroup/pyuvsim) (which sacrifices speed for excellent precision), and [RIMEz](https://github.com/upenneor/rimez) (which decomposes the sky into spherical harmonics rather than discrete points). `WODEN` was designed with the Australian MWA Epoch of Reionisation (EoR) processing pipeline in mind, which uses a calibration and foreground removal software called the `RTS` [@Mitchell2008] in search of signals from the very first stars (see @Yoshiura2021 for a recent use of this pipeline). The `RTS` creates a sky model using the same formalism above, however the code is not optimised enough to handle the volume of sources to simulate the entire sky. To test the `RTS` method of sky generation, we therefore needed a fast and discretised method. Another excellent CUDA accelerated simulation package, [OSKAR](https://github.com/OxfordSKA/OSKAR), addresses these two points. However, the `RTS` also generates parts of the sky model via shapelets (see @Line2020 for an overview), which `OSKAR` cannot. Furthermore, in real data, the precession of the Earth's rotational axis causes sources to move from the sky coordinates as specified in the RA, DEC J2000 coordinate system. The `RTS` is designed to undo this precession, and so a simulation fed into the `RTS` should *contain* precession. `WODEN` adds in this precession using the same method as the `RTS` to be consistent. This unique combination of CUDA, shapelet foregrounds, the MWA FEE primary beam, along with source precession, created the need for `WODEN`. These effects should not preclude other calibration packages from using `WODEN` outputs however, meaning `WODEN` is not limited to feeding data into the `RTS` alone.

# Example application

In @Line2020, we compared two methods to model Fornax A: a combination of point and elliptical Gaussians, compared to shapelets (see \autoref{fig:ForA}). We were able to quickly compare the computational efficiency of the methods using a desktop, and comment on their respective strengths and weaknesses in regard to foreground removal for EoR purposes. Furthermore, as we could control the simulations, we could compare the methods in the absence of other processing systematics that are present in the real data from the MWA, which dominated the comparison when using the `RTS` alone.

![Two methods to simulate Fornax A visibilities are compared here (both imaged using `WSClean` [@Offringa2014; @Offringa2017]), with point and elliptical Gaussians on the left, and shapelets on the right.\label{fig:ForA}](FornaxA_model_comparison.png)

# Documentation

The documentation for `WODEN` can be found [here on readthedocs](https://woden.readthedocs.io/en/latest/), including a detailed installation guide, ways to test a local installation, details of the calculations `WODEN` makes under the hood, and worked examples, which are also included in the `github` repo.

# Acknowledgements

I acknowledge direct contributions from Tony Farlie (who taught me how pointer arithmetic works in C) and indirect contributions from Bart Pindor and Daniel Mitchell (through their work in the `RTS` and through advising me on CUDA). I'd like to thank Chris Jordan who acted as a sounding board as I learnt C and CUDA.

This research was supported by the Australian Research Council Centre of Excellence for All Sky Astrophysics in 3 Dimensions (ASTRO 3D), through project number CE170100013.

# References
