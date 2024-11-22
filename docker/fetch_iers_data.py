from astropy.utils.iers import IERS_Auto

# Pre-fetch the IERS-A data to ensure it's available
iers_a = IERS_Auto.open()