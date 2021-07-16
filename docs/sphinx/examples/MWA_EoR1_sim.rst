.. _`this googledoc link to woden-srclist_pumav3.txt`: https://drive.google.com/file/d/1GFnQPXVGsS_7eE5EKTuRp6naNO6IHNFI/view?usp=sharing

MWA EoR1 simulation
====================

.. note:: Running the simulation and making all the images will take up around 1.8 GB storage.

You'll need to download the skymodel from `this googledoc link to woden-srclist_pumav3.txt`_ and put it in the correct directory. If you're comfortable with ``wget`` you can do::

  $ cd WODEN/examples/MWA_EoR1
  $ wget 'https://docs.google.com/uc?export=download&id=1GFnQPXVGsS_7eE5EKTuRp6naNO6IHNFI' -O woden-srclist_pumav3.txt

This skymodel contains over 300,000 sources, and is based on the GLEAM catalogue, with some embellishment. The simulation we're going to run is of the MWA 'EoR1' field, centred at RA, Dec = :math:`60^\circ, -30^\circ`. This field contains Fornax A, and is a good way to demonstrate a sky with point, Gaussian, and shapelet models, as well as the effect of the MWA FEE primary beam. To run the simulation, simply run::

  $ ./MWA_EoR1_simulation.sh

Which contains the command::

  time run_woden.py \
    --ra0=60.0 --dec0=-27.0 \
    --num_freq_channels=16 --num_time_steps=14 \
    --freq_res=80e+3 --time_res=8.0 \
    --cat_filename=woden-srclist_pumav3.txt \
    --metafits_filename=../metafits/1136380296_metafits_ppds.fits \
    --band_nums=1,2,3,4,5 \
    --output_uvfits_prepend=./data/MWA_EoR1 \
    --primary_beam=MWA_FEE \
    --sky_crop_components

Running this took 5 mins 30 seconds on my GPU. I've reduced the time and frequency resolution as specified in the ``metafits`` file to keep the size of the outputs smaller on your machine. If you wanted to run the full resolution data of this observation, (2s, 40kHz), you can just remove the ``--num_freq_channels, --num_time_steps, --freq_res, --time_res`` arguments.

I've included two imaging commands::

  $ ./MWA_EoR1_smaller_imaging.sh
  $ ./MWA_EoR1_larger_imaging.sh

The 'larger' image is 6000 by 6000 pixels and has large *w*-terms, so it took about an hour on my desktop. If you don't want to hang around, just run the smaller imaging. The smaller image looks like this:

.. image:: MWA_EoR1_plot_smaller.svg
   :width: 600px

Here I've blown up the colour scale just to highlight how many sources there are (especially FornaxA booming in the corner). The CLEAN isn't great here as I haven't made the image big enough (you can see some aliasing of Fornax A around the image).

If you have the patience to make the bigger image, it looks like this:

.. image:: MWA_EoR1_plot_larger.svg
   :width: 600px

On the left we see the full image, which clearly shows the main lobe of the MWA primary beam. On the right I have zoomed into north of the main lobe, and you can see sources sat in the northern beam sidelobe of the primary at the top, and the edge of the primary lobe at the bottom. This is as expected due to the MWA primary beam shape (see the :ref:`MWA Fully Embedded Element` section if you are unfamiliar with the beam).
