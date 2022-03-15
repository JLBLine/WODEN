import mwa_hyperbeam
beam = mwa_hyperbeam.FEEBeam()


delays = [6, 4, 2, 0, 8, 6, 4, 2, 10, 8, 6, 4, 12, 10, 8, 6]


az = 3.7997899
za = 0.3080986

jones = beam.calc_jones(az, za, 150e6, delays, [1]*16, True, True)


print(jones)
