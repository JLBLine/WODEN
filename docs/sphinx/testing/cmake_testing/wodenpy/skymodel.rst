``skymodel``
=========================
Tests for the functions in ``wodenpy.skymodel``. These functions handle reading in sky models from various formats, counting components, and splitting sky models into chunks. As ``WODEN`` is intended to simulate all-sky sky models with a large number of components, these models are lazy-loaded. The code below therefore goes to great pains to map how this loading will occur. 

.. warning:: From version 2.0 onward, the ``LoBES``-style FITS sky model format is preferred over others. At the time of writing, the FITS format is limited to Stokes I only. Internally, the code reads in all formats and converts them to the FITS format. The same mapping and chunking code is performed upon them for simplicity. The skymodel machinery to read in/store Stokes Q,U,V still exists however, so future versions of ``WODEN`` may support these formats.


test_create_skymodel_chunk_map.py
*********************************************


test_crop_below_horizon.py
*********************************************


test_map_chunk_pointgauss.py
*********************************************


test_map_chunk_shapelets.py
*********************************************


test_read_FITS_skymodel_chunk.py
*********************************************


test_read_skymodel_chunk.py
*********************************************


test_read_text_skymodel_chunk.py
*********************************************


test_read_yaml_radec_count_components.py
*********************************************


test_read_yaml_skymodel_chunk.py

