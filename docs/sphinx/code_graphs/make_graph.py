import pyan



def make_dot(files, output):
    callgraph = pyan.create_callgraph(files, format='dot')

    with open(output, 'w') as f:
        f.write(callgraph)


if __name__ == "__main__":
    woden_dir = "/home/jack-line/software/WODEN"
    
    files = [f'{woden_dir}/scripts/run_woden.py', f'{woden_dir}/wodenpy/*/*.py']
    make_dot(files, "wodenpy_all.dot")

    files = [f'{woden_dir}/scripts/run_woden.py',
            f'{woden_dir}/wodenpy/array_layout/create_array_layout.py',
            f'{woden_dir}/wodenpy/array_layout/precession.py']
    make_dot(files, "wodenpy_array_layout.dot")

    files = [f'{woden_dir}/scripts/run_woden.py',
            f'{woden_dir}/wodenpy/observational/calc_obs.py']
    make_dot(files, "wodenpy_observational.dot")
    
    files = [f'{woden_dir}/scripts/run_woden.py',
            f'{woden_dir}/wodenpy/skymodel/read_skymodel.py',
            f'{woden_dir}/wodenpy/skymodel/read_fits_skymodel.py'
            f'{woden_dir}/wodenpy/skymodel/woden_skymodel.py',
            f'{woden_dir}/wodenpy/skymodel/chunk_sky_model.py']
    make_dot(files, "wodenpy_skymodel.dot")
    
    
    files = [f'{woden_dir}/scripts/run_woden.py',
            f'{woden_dir}/wodenpy/use_libwoden/array_layout_struct.py',
            f'{woden_dir}/wodenpy/use_libwoden/beam_settings.py'
            f'{woden_dir}/wodenpy/use_libwoden/create_woden_struct_classes.py',
            f'{woden_dir}/wodenpy/use_libwoden/shapelets.py',
            f'{woden_dir}/wodenpy/use_libwoden/skymodel_structs.py',
            f'{woden_dir}/wodenpy/use_libwoden/use_libwoden.py',
            f'{woden_dir}/wodenpy/use_libwoden/visibility_set.py',
            f'{woden_dir}/wodenpy/use_libwoden/woden_settings.py']
    make_dot(files, "wodenpy_use_libwoden.dot")
    
    files = [f'{woden_dir}/scripts/run_woden.py',
            f'{woden_dir}/wodenpy/uvfits/wodenpy_uvfits.py']
    make_dot(files, "wodenpy_uvfits.dot")
    
    files = [f'{woden_dir}/scripts/run_woden.py',
            f'{woden_dir}/wodenpy/wodenpy_setup/git_helper.py',
            f'{woden_dir}/wodenpy/wodenpy_setup/run_setup.py']
    make_dot(files, "wodenpy_wodenpy_setup.dot")
    
    