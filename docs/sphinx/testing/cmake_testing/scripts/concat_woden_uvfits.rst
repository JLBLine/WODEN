``concat_woden_uvfits.py``
===========================
``concat_woden_uvfits.py`` is supposed to concatenate multiple coarse band ``uvfits`` into one (aka 24 coarse bands of 32 channels in one ``uvfits`` of length 384 channel). This is useful to simulate across 24 GPUs but combine outputs into one ``uvfits`` that can be input to calibration software.


test_concat_woden_uvfits.py
****************************
This test writes out some example ``uvfits`` files with known contents, creates a mock of a command line argument object, and runs ``concat_woden_uvfits.main`` to create concatenated outputs. Simply we just make 4 uvfits files and stich em together. There are two options of reverse the X and Y polarisations, and halving the power of the visibilities (for different Stokes conventions). So run the test 4 times to test all options.