Call Graphs
-------------------

Ok, here is an attempt to map out the structure of ``WODEN``. These call graphs split into two sets: one for the Python side, the other for the C/GPU side. They are generated using the ``WODEN/docs/sphinx/code_graphs/run_make_graph.sh``, via the ``pyan`` Python module, and ``Doxygen``.


Python Call Graphs
~~~~~~~~~~~~~~~~~~~~

Firstly, the entire ``wodenpy`` call graph looks like this. Obviously this is a mess, so scroll on for individual submodules call graphs.

.. image:: wodenpy_all.svg
   :width: 100%
   :align: center

``wodenpy/array_layout``:

.. image:: wodenpy_array_layout.svg
    :width: 100%
    :align: center


``wodenpy/observational``:

.. image:: wodenpy_observational.svg
    :width: 100%
    :align: center


``wodenpy/skymodel``:

.. image:: wodenpy_skymodel.svg
    :width: 100%
    :align: center


``wodenpy/use_libwoden``:

.. image:: wodenpy_use_libwoden.svg
    :width: 100%
    :align: center


``wodenpy/uvfits``:

.. image:: wodenpy_uvfits.svg
    :width: 100%
    :align: center


``wodenpy/wodenpy_setup``:

.. image:: wodenpy_wodenpy_setup.svg
    :width: 100%
    :align: center



C/GPU Call Graphs
~~~~~~~~~~~~~~~~~~~~
Eventually, ``run_woden.py`` calls the C/GPU function ``calculate_visibilities``. This is the call graph for ``calculate_visibilities``; note that ``dot`` has truncated some boxes (which are rendered red), as it has a maximum width. Scroll further for a second graph that starts the function ``source_component_common``, which includes the missing calls.

.. image:: calculate__visibilities_8cpp_a40d3a09b3be87b97684afc84864248fd_cgraph.png
   :width: 100%
   :align: center


Call graph starting at ``source_component_common``:

.. image:: source__components_8cpp_a577483371b5f50c345ac3a9d141b09aa_cgraph.png
   :width: 100%
   :align: center
