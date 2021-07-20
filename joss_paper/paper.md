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

`WODEN` is designed to simulate the response of a class of telescope known as an interferometer, producing output "visibilities" for a given astrophysical sky model. The sky model can be specified in dirac-delta like functions on the sky (known in the field as "point sources") elliptical Gaussian models, or built out of "shapelet" basis functions (allowing complicated morphologies to be created). Users are able to input a bespoke layout for the interferometer, vary a number of observational parameters including time of day, length of observation and frequency coverage, and select from a number of predefined primary beams which encode the response of the receiving elements of an interferometer.

Simulated observations allow us to test other software packages that are designed to calibrate and analyse real interferoteric data, including verifying expected behaviour with known inputs, and trialling new sky modelling techniques.

The core functionality of `WODEN` is written in CUDA as interferometric simulations are computationally intensive but embarrassingly parallel. The compute performance of CUDA allows for large-scale simulations to be run including emission from all directions in the sky. This is paramount for interferometers with a widefield of view such as the Murchison Widefield Array (MWA, @Tingay2013).

A Python wrapper is used to take advantage of community packages such as astropy(get citation) and to present a user-friendly interface to `WODEN`.

# Statement of need
Stick in the measurement equation, and say how it scales with number of baselines / frequencies / time steps / number of sources on sky.

The MWA has receiving elements made up of dipole antennas, which are sensitive to substantial areas of the sky at any time. The MWA combines 16 dipoles in a four by four grid to make a single receiving element, creating a distinct pattern on the sky or "primary beam".

Explain all sky at MWA resolution means millions up millions of point sources.

Mention:

* [OSKAR](https://github.com/OxfordSKA/OSKAR)
* [pyuvsim](https://github.com/RadioAstronomySoftwareGroup/pyuvsim)
* [RIMEz](https://github.com/upenneor/rimez)

State the unique part of `WODEN` is the combination of CUDA, shapelets, and the MWA FEE beam, along with being able to plug straight into the `RTS` [@Mitchell2008].

<!-- An interferometer is made from up to hundreds of receivers, which range from large dishes to small dipole antennas (each called a receiving element in the following). For astronomy at low radio-wavelengths (radiation of the range ~50 to 500$\,$MHz), an interferometer is typically necessary to resolve interesting astrophysical sources on the sky.

Incoming radio waves from astrophysical sources induce a voltage in every receiving element, which is then cross-correlated with the signal from every other receiving element. A pair of receiving elements is called as 'baseline', and the cross-correlated signals 'visibilities'; both of these quantities are frequency dependent. If we have $N$ receiving elements, at any time and frequency step, we produce $N(N-1)) / 2$ -->


<!-- Many astrophysical sources can be observed at low radio-wavelengths, including: black holes undergoing accretion; spinning neutron stars; cosmic rays spiralling around magnetic fields lines; red-shifted signals from neutral hydrogen early in the Universe, to name a few. Visibilities can be used to reconstruct an image of the radio-loud night sky (amongst many other things),  allowing us to study these astrophysical sources.

These  -->



<!-- As such, radio emission that comes areas of the sky far from




An example is

We need it to feed into .

That citation didn't work, does the example one work? [@Tingay2013] -->

# Example application

`WODEN` paper @Line2020

# Documentation

The documentation for `WODEN` can be found [on readthedocs here](https://woden.readthedocs.io/en/latest/), including a installation guide, ways to test a local installation, details of the calculations `WODEN` makes under the hood, and worked examples which are also included in the `github` repo.

# Acknowledgements

I acknowledge direct contributions from Tony Farlie and indirect contributions from Bart Pindor and Daniel Mitchell.

# References
