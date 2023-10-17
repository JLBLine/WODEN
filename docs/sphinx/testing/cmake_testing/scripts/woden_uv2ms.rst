``woden_uv2ms.py``
===========================
``woden_uv2ms.py`` converts a uvfits file to a measurement set via ``pyuvdata``.


test_add_woden_uvfits.py
***************************
This test writes out some example ``uvfits`` files, and attempts to convert them into measurement sets using ``woden_uv2ms.py``. It then simply checks the file exists, as we assume error handling is done by ``pyuvdata``.