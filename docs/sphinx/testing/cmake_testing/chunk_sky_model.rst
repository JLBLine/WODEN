``chunk_sky_model``
=========================
Tests for the functions in ``WODEN/src/chunk_sky_model.c``. These functions
handle splitting the sky model up into chunks that can fit into GPU memory
during simulation. There is a variable (which can be controlled by the user
via ``--chunking_size``) that controls the maximum number of measurement
equations that can be calculated simultaneously, so the some of the functions
below need to know the number of baselines, frequencies etc to calculate
the maximum number of COMPONENTs to put into each chunk.


``test_null_comps.c``
****************************
Tests the function ``chunk_sky_model::null_components``.
This function sets either POINT, GAUSSIAN, or SHAPELET COMPONENT attributes of
the sky model to ``NULL``. Tested here by setting up a dummy populated sky model,
and making sure each function sets the correct attributes to ``NULL``.
Also tests that other COMPONENT attributes are left as they are.


``test_fill_chunk_src_with_pointgauss.c``
***********************************************
Tests ``chunk_sky_model::fill_chunk_src_with_pointgauss``, which calculates
how many POINT and GAUSSIAN COMPONENTs should go into each chunk. Tests by
trying 12 different combinations of number of POINTS, number of GAUSSIANs,
chunk size, number time steps, number baselines, and number of frequencies.
Iterates over the function and checks that each chunked sky model has the
correct number of output COMPONENTs of each type. Checks by setting attributes
of the sky model to the array index from the original full sky model, and checks
that every attribute is being split correctly into the resultant chunked sky
models. Furthermore, as azimuth and zenith angles are stored for all time steps
in the sky model, the dummy sky model is setup with repeating/tiled arrays
to make sure the correct az/za coords are being copied from the full sky
model into the cropped sky models.

``test_fill_chunk_src_with_shapelets.c``
***********************************************
Tests ``chunk_sky_model::fill_chunk_src_with_shapelets``, which handles
chunking the SHAPELET COMPONENTs of the sky model. Works the same as pointgauss
in that it sets sky model attributes to array indexes, to check whether the
correct values are being chunked into the smaller chunked sky models. On top
of that, due to the way SHAPELETs can have multiple basis functions per
COMPONENT, has an extra layer of sky model trickery to have repeating indexes
inside the basis function attributes, which can then be traced and checked
once the chunking has been performed. Tests here vary the number of SHAPELETs,
the number of basis functions per SHAPELET, and the number of time steps.

``test_create_chunked_sky_models.c``
***********************************************
``chunk_sky_model::create_chunked_sky_models`` uses all functions above to take
in a full sky model and created an array of chunked sky models. This test
runs the same testing for both functions above with a sky model containing
various combinations of POINT, GAUSSIAN, and SHAPELET COMPONENTs, as well
as different chunking sizes and time settings.


``test_remap_source_for_gpu.c``
***********************************************
At some point, we have to take all the information in the chunked sky models
and copy it across to the GPU, where we do all the heavy lifting. As there
can be millions of COMPONENTs (if simulating a full sky healpix projection for
example), we likely won't be able to copy across all information to the GPU.
``chunk_sky_model::remap_source_for_gpu`` takes a ``source_t`` struct as spat
out by ``chunk_sky_model::create_chunked_sky_models``, and does any necessary
memory allocation and mapping of information, such that the new remapped
``source_t`` can just be straight copied across into GPU memory. This way, we
can loop over many chunked sky models, grab how much memory we need at the
time, and then free it once we're done. If we did this remapping before
iterating over chunks, we double the necessary host memory (or have to creatively
free the full-sky model memory, which I cba to do.)

This test loops over many chunked sources by first running
``chunk_sky_model::create_chunked_sky_models`` on the same set of sky models
as tested in ``test_create_chunked_sky_models.c``. It the runs
``chunk_sky_model::remap_source_for_gpu`` for each chunk, checks that the
outputs are as expected (I've basically set all the ra,dec etc etc in the
original full sky models to indexes or arrays, to make testing expected outcomes
straight forward). After each remapping, ``chunk_sky_model::free_remapped_source_for_gpu``
is called to free up memory which is assigned in each loop, so we also test this
freeing function at the same time.
