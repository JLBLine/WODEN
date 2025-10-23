from casacore.tables import table
from wodenpy.array_layout.create_array_layout import convert_ecef_to_enh
import numpy as np

lofar_lat = 52.905329712
lofar_long = 6.867996528

MWA_LAT=-26.703319405555554
MWA_LONG=116.67081523611111


def get_station_layout_in_enh(station_id : int,  ms_path : str):
    """
    Get the station layout in east,north,height coordinates for a given station ID from a Measurement Set.
    Args:
        station_id (int): The ID of the station to retrieve the layout for.
        ms_path (str): Path to the Measurement Set.
    Returns:
        east_station (np.ndarray): East coordinates of antennas in the station.
        north_station (np.ndarray): North coordinates of antennas in the station.
        height_station (np.ndarray): Height coordinates of antennas in the station.

    """
    
    with table(ms_path + '/ANTENNA', readonly=True) as t: 
        num_stations = len(t)
        ant_locations = np.array([t.getcell('POSITION', ant) for ant in range(num_stations)])
    
    with table(ms_path+'::LOFAR_ANTENNA_FIELD', readonly=True) as antenna:
        xyz_station = antenna.getcell('ELEMENT_OFFSET', station_id) + ant_locations[station_id, :]
        
        east_station, north_station, height_station = convert_ecef_to_enh(xyz_station[:,0],
                                            xyz_station[:,1], xyz_station[:,2],
                                            np.radians(lofar_long), np.radians(lofar_lat))
        
    return east_station, north_station, height_station

if __name__ == "__main__":
    east_001, north_001, height_001 = get_station_layout_in_enh(0, 'pointed_LBA.ms')
    east_302, north_302, height_302 = get_station_layout_in_enh(20, 'pointed_LBA.ms')
    
    np.savez('LOFAR_LBA_station_layouts.npz',
             east_001=east_001, north_001=north_001, height_001=height_001,
             east_302=east_302, north_302=north_302, height_302=height_302)
    
    
    with table('create_OSKAR-SKA_ms/OSKAR-SKA-layout.ms' + '/ANTENNA') as t: 
        
        num_ants = len(t)
        ant_locations = np.array([t.getcell('POSITION', ant) for ant in range(num_ants)])
        ##convert from ECEF to ENH, as WODEN starts with enh coords
        east, north, height = convert_ecef_to_enh(ant_locations[:,0],
                                        ant_locations[:,1], ant_locations[:,2],
                                        np.radians(MWA_LONG),
                                        np.radians(MWA_LAT))
        
    np.savez('OSKAR-SKA-layout.npz', east=east, north=north, height=height)