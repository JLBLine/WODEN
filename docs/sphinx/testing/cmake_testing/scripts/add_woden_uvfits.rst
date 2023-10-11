``add_woden_uvfits.py``
===========================
``add_woden_uvfits.py`` is supposed to sum two different uvfits files, which is useful for example to add foregrounds onto a EoR signal simulation.


test_add_woden_uvfits.py
***************************
This test writes out some example ``uvfits`` files with known contents, creates a mock of a command line argument object, and runs ``add_woden_uvfits.main`` to create combined outputs. These outputs are checked to be the sum of the inputs. A number of input arguments that should cause failures, including attempting to add to uvfits with different dimensions, are also tested.