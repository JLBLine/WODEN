``skymodel``
=========================
Tests for the functions in ``wodenpy.skymodel``. These functions handle reading in sky models from various formats, counting components, and splitting sky models into chunks. As ``WODEN`` is intended to simulate all-sky sky models with a large number of components, these models are lazy-loaded. The code below therefore goes to great pains to map how this loading will occur. 

test_crop_below_horizon.py
*********************************************
Given input RA/Dec coords and an LST,Latitude, ``crop_below_horizon`` crops everything below the horizon. There are also two modes of cropping: one that crops by component, and one by source. Test by running using three different LST, each cropping either by component or by source. Test against stored expected outcomes.

test_create_skymodel_chunk_map.py
*********************************************
``create_skymodel_chunk_map`` makes a map of how to split a given number of sky model components into bit size chunks that fit onto the GPU, and also fit in RAM. This splitting is done based on the number of visibilities being simulated and user settings. Test by running with seven different combinations of inputs, varying the numbers of points, gaussians, shapelets, maximum visibilities, number time steps and frequencies. Test against a separate set of test functions that produce the desired outputs.

test_map_chunk_pointgauss.py
*********************************************
This function is used by ``create_skymodel_chunk_map`` to correctly map the point and gaussian components into chunks, given a maximum number of components per chunk. Test by inputting multiple combinations of inputs, varying the numbers of points, gaussians, maximum visibilities, number time steps and frequencies, different numbers of list flux-type values, and different types and numbers of polarised components (which effects how much ctype "malloc"ing needs to be done). Test against a separate set of test functions written in C that produce the desired outputs.

test_map_chunk_shapelets.py
*********************************************
This function is used by ``create_skymodel_chunk_map`` to correctly map the shapelets components into chunks, given a maximum number of shapelet coefficients per chunk (shapelets are chunked by the number of shapelet basis functions to use, rather than number of components). Test by inputting 4 different combinations of inputs, varying the numbers of shapelet components and coeffs, number time steps and frequencies, different numbers of list flux-type values, and different types and numbers of polarised components (which effects how much ctype "malloc"ing needs to be done). Test against a separate set of test functions written in C that produce the desired outputs.

test_read_FITS_skymodel_chunk.py
*********************************************
Tests reading in a LoBES-style FITS sky model.

This actually tests both ``read_fits_radec_count_components`` and ``read_fits_skymodel_chunks``, by writing a FITS sky model, and then reading that model in. We can then check the resultant chunked sky model matches the values we wrote to the FITS sky model, and is chunked how we expect. This test uses ``crop_below_horizon`` and 
``create_skymodel_chunk_map`` to do the sky model cropping and mapping for the sky model chunking. Functions in the test code ``common_skymodel_test.py`` and ``read_skymodel_common.py`` have been setup to create expected outcomes for the test to compare against. Test by running multiple combinations of inputs, each which create a different number of RA/Dec coords in the catalogue, have different LSTs, different chunking sizes, and different numbers of components per source. For each RA/Dec location, three point, gaussian, and shapelet components are included, each with a flux models of either a power law, curved power law, or list type flux. Furthermore, all types of polarisation models are included, at varying cadences within the sky model. This means all combinations possible in the sky model are tested. Each parameter is some fraction/multiple of the index of the component in the sky model, allowing for simpler predictions of what the resultant chunked component values should be.

test_read_text_skymodel_chunk.py
*********************************************
This does exactly the same as the ``test_read_FITS_skymodel_chunk.py``, but creates and reads a WODEN-style text sky model. This is only really for backwards compatibility, as the WODEN sky model is deprecated in favour of the FITS sky model.


test_read_yaml_skymodel_chunk.py
*********************************************
This does exactly the same as the ``test_read_FITS_skymodel_chunk.py``, but creates and reads a hyperdrive-style yaml sky model. `WODEN` only reads Stokes I from the ``yaml`` format, so this again is for backwards compatibility.

test_read_skymodel_chunk.py
*********************************************
This tests ``read_skymodel_chunk``, which is a wrapper function that detects the type of sky model being read in, and calls the appropriate function to read it in. This is tested by running the same tests as in ``test_read_FITS_skymodel_chunk.py``, ``test_read_text_skymodel_chunk.py``, and ``test_read_yaml_skymodel_chunk.py``.