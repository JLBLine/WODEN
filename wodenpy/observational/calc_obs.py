"""Functions for calculating observational parameters based on date and location."""
from astropy.time import Time, TimeDelta
import numpy as np
from astropy.coordinates import EarthLocation
from astropy import units as u

def calc_jdcal(date):
    """Takes a string calendar date-time and returns julian date by using
    `astropy.time.Time`_, so date can be formatted anyway `astropy.time.Time`_
    accepts. Returns the date in two parts, an integer day, and a fractional day.
    The header of a uvfits file takes the integer part, and the fractional
    part goes into the DATE array

    .. _astropy.time.Time: https://docs.astropy.org/en/stable/time/

    Parameters
    ----------
    date : string
        UTC date/time (e.g. in format YYYY-MM-DDThh:mm:ss or similar)

    Returns
    -------
    jd_midnight : float
        midnight value of the julian date, which is the UVFITS standard
        to go into the header
    jd_fraction : float
        remaining fraction of the julian date

    """

    t = Time(date)
    jd = t.jd
    
    ##Get the midnight value of the julian date
    midnight = Time(t.datetime.date().isoformat())
    jd_midnight = midnight.jd
    
    jd_fraction = (jd - jd_midnight)

    ##The header of the uvdata file takes the integer, and
    ##then the fraction goes into the data array for PTYPE4
    return jd_midnight, jd_fraction

def get_uvfits_date_and_position_constants(latitude=None,longitude=None,
                                           date=None,height=None):
    """
    Returns a number of date and time based values that are needed for uvfits
    headers. For the given Earth location and UTC date return the local sidereal
    time (deg), the Greenwich sidereal time at 0 hours on the given date (deg),
    the rotational speed of Earth on the given date (in degrees per day), and
    the difference between UT1 and UTC.
    Uses `astropy.time.Time`_ and `astropy.coordinates.EarthLocation`_ to make
    the calculations.

    .. _astropy.time.Time: https://docs.astropy.org/en/stable/time/
    .. _astropy.coordinates.EarthLocation: https://docs.astropy.org/en/stable/api/astropy.coordinates.EarthLocation.html?highlight=EarthLocation


    Parameters
    ----------
    latitude : float
        Latitude of location on Earth (deg)
    longitude : float
        Longitude of location on Earth (deg)
    height : float
        Height above sea level of location on Earth (m)
    date : string
            UTC date/time (in format YYYY-MM-DDThh:mm:ss)

    Returns
    -------
    LST_deg : float
        Local sidereal time (degrees)
    GST0_deg : float
        Greenwich sidereal time at 0 hours on the given date (degrees)
    DEGPDY : float
        Rotational speed of Earth on the given date (degrees per day)
    ut1utc : float
        Difference between UT1 and UTC (secs)
    """

    ##Setup location
    observing_location = EarthLocation(lat=latitude*u.deg, lon=longitude*u.deg, height=height)
    ##Setup time at that location

    lst_type = 'mean'

    observing_time = Time(date, scale='utc', location=observing_location)
    ##Grab the LST
    LST = observing_time.sidereal_time(lst_type)
    LST_deg = LST.value*15.0

    ##uvfits file needs to know the greenwich sidereal time at 0 hours
    ##on the date in question
    zero_date = date.split('T')[0] + "T00:00:00"
    zero_time = Time(zero_date, scale='utc', location=observing_location)
    GST0 = zero_time.sidereal_time(lst_type, 'greenwich')
    GST0_deg = GST0.value*15.0

    ##It also needs to know the rotational rate of the Earth on that day, in
    ##units of degrees per day
    ##Do this by measuring the LST exactly a day later

    date_plus_one_day =  observing_time + TimeDelta(1*u.day)
    LST_plusone = date_plus_one_day.sidereal_time(lst_type)

    LST_plusone_deg = LST_plusone.value*15.0

    DEGPDY = 360.0 + (LST_plusone_deg - LST_deg)

    ut1utc = float(observing_time.delta_ut1_utc)

    return LST_deg, GST0_deg, DEGPDY, ut1utc