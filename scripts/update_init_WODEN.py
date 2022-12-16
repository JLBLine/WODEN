"""Script that reads in the library locations of ERFA, HDF5, JSONC and PAL
from the CMakeCache.txt file, then appends those locations to the
LD_LIBRARY_PATH inside of init_WODEN.sh so `woden` can easily find them at
run time"""


def grab_path(line, paths):
    """Takes a line from a CMakeCache.txt file, and splits and strips
    to find the path to the library so we can append to LD_LIBRARY_PATH.
    Appends the path found into the list `paths`"""

    path = line.split('=')[1]
    path = '/'.join(path.split('/')[:-1])

    paths.append(path)

    return paths


if __name__ == '__main__':
    ##Read in the CMakeCache.txt file, which contains all the paths to the
    ##libraries `woden` needs
    with open('CMakeCache.txt','r') as cmakecache:
        lines = [line for line in cmakecache.read().split('\n') if line != '']

    ##Something to hold all the paths in
    paths = []

    ##Following libraries are the main dependencies that might live in
    ##unusual locations. Other libraries are system libraries that should
    ##be handled by system owner
    for line in lines:
        if 'ERFA_LIB' in line:
            paths = grab_path(line, paths)

        elif 'HDF5_LIB' in line:
            paths = grab_path(line, paths)

        elif 'JSONC_LIB' in line:
            paths = grab_path(line, paths)

        elif 'PAL_LIB' in line:
            paths = grab_path(line, paths)

    with open('init_WODEN.sh', 'w') as outfile:
        outfile.write('##This line finds the current directory at sets the env variable WODEN_DIR\n')
        outfile.write('export WODEN_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"\n')
        outfile.write('##This adds the line to PATH\n')
        outfile.write('export PATH=$WODEN_DIR:$PATH\n')

        ##Taking a set of paths avoids duplication. For each unique path,
        ##add the path to LD_LIBRARY_PATH so they can be found at runtime
        outfile.write('##Add library paths to LD_LIBRARY_PATH so the can be found at runtime\n')
        for path in set(paths):
            outfile.write(f"export LD_LIBRARY_PATH={path}/:$LD_LIBRARY_PATH\n")
