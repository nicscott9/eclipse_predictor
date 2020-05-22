import matplotlib.pyplot as plt
import numpy as np
import math 
import datetime
import julian
from astropy.time import Time
import astropy.units as u
import astropy.constants as const
import astropy.coordinates as coord
from astropy.coordinates import ICRS, Galactic, FK4, FK5  # Low-level frames
from astropy.coordinates import Angle, Latitude, Longitude  # Angles

#https://docs.astropy.org/en/stable/time/ 

#put in target name
#target="ip_peg"
#or coords
ra = "15:02:41.00"
dec = "+33:34:24.00"

#warning: this code is not yet using the proper motions in the time calculation
pmra = -77.759 
pmdec = -36.881
pmtime = Time(2000, format='decimalyear') #epoch

parallax = 5.3552

#enter observed eclipse time and period
t_0 = 53799.1406073
p = 0.058908615
t = Time(t_0, format='mjd')  

#enter telescope location
telescope = "gemini_north"
#telescope = "gemini_south"
#telescope = "kitt peak"
#l = coord.EarthLocation.get_site_names()

#enter day of interest in mjd
mjd = 58990
interesttime = Time(mjd, format='mjd') 

#enter number of eclipses you want to display from day of interest 
n_eclipses = 10.

#calc number of eclipses from t_0 to day of interest
n_since = (mjd-t_0)/p
num_eclipses = math.floor(n_since+n_eclipses)

#generate array of eclipses from t_0 to day of interest+n_eclipses
p_array = np.arange(0, num_eclipses+1)
eclipses = t_0+(p_array*p) #BJD eclipses from t_0 to time of interest + n_eclipses

# nt = Time.now()
# ut = Time(datetime.utcnow(), scale='utc')

#calculate obj coords
#obj=coord.SkyCoord.from_name(target)  
obj = coord.SkyCoord(ra, dec, 1000/parallax, unit=(u.hourangle, u.deg, u.pc), pm_ra_cosdec=pmra * u.mas/u.yr, pm_dec=pmdec * u.mas/u.yr, obstime=pmtime,frame='fk5')
newcoords = obj.apply_space_motion(interesttime) 
#print(newcoords)

#calc light travel time
times = Time(mjd, format='mjd', location=coord.EarthLocation.of_site(telescope))  
ltt_bary = times.light_travel_time(newcoords, ephemeris='jpl')  

#eclipses are in bjd so should subtract ltt to get them in mjd
mjdeclipses = eclipses - ltt_bary
#eclipse_times= Time(mjdeclipses, format='jd', scale='utc')  
#utc_eclipses = eclipse_times.to_datetime()

print('obj coords = '+ra+' '+dec)
print('pm corrected coords = '+newcoords.to_string('hmsdms'))
print('location = '+telescope)
print('prev obs eclipse = '+str(t_0)+' UTC:'+str(t.to_datetime()))
print('obj period = '+str(p)+'d'+' = '+str(p*24*60)+' mins')
print('entered time = '+str(times.to_datetime()))
print('light travel time = '+str(ltt_bary))
print('next eclipses (mjd) from location = ')
print(mjdeclipses[-10:])

#t = mjdeclipses[-1]
#times=t-mjd
#timesseries = mjdeclipses - mjd
#print(Time(times, format='datetime',scale='tdb'))

print('for '+Time(interesttime, out_subfmt='date_hm').iso+' the next '+str(int(n_eclipses))+' eclipses will occur at UTC times: ')
for i in range(0, int(n_eclipses)):
    t = mjdeclipses[-int(n_eclipses)+i]
    times = t - math.floor(mjd)
    print(Time(times, format='datetime'))
