from skyfield.api import load, wgs84, Topos
from skyfield.api import Distance
from skyfield.toposlib import ITRSPosition

# Most of the calculations of this code come from the skyfield documentation: https://rhodesmill.org/skyfield/
# This is the personal website of Brandon Rhodes, skyfield library's creator
# We have mainly used the following subsections: https://rhodesmill.org/skyfield/earth-satellites.html (Earth Satellites)
# and https://rhodesmill.org/skyfield/positions.html (Positions)

# Some useful theory used in the calculations:
# ECEF = Earth-centered Earth-fixed reference system. It is inertial relative to the Earth, assimilated to ITRS
# GCRS = Geocentric Celestial Reference System
# Both are different coordinate reference frames used to define positions of variables
# Conversion between them will be necessary for our calculations

stations_url = 'http://celestrak.org/NORAD/elements/stations.txt'
satellites = load.tle_file(stations_url)

by_name = {sat.name: sat for sat in satellites}
satellite = by_name['ISS (ZARYA)']
ts = load.timescale()


#SPECIFIC LOCATION
# We have (latX,lonX) = (37.833959, -75.487833)
# This is our launch point, we are going to find the moment during the day where ISS is closer to the point


def shortest_distance(lat_given, lon_given, date_list): #date list = year,month,day
    dist_min = 500000 #QU: voir quelle valeur mettre
    latmin=0
    lonmin=0
    altkmin=0
    for s in range(0,86400): #s determine la seconde
        hour = s//3600
        p = s%3600
        minute = p//60
        seconde = p%60
        t1 = ts.utc(date_list[0],date_list[1],date_list[2],hour,minute,seconde)
        #position iss:
        geocentric = satellite.at(t1)
        lat, lon = wgs84.latlon_of(geocentric)
        alt_iss = wgs84.height_of(geocentric).m #returns a Distance

        #To find the distance bewteen ISS and our given point, we first need to convert the coordinates of the point into gcrs coordinates
        #This is done to be able to use Skyfield's attribute of distance, available for positions in the GCRS reference frame
        #point en grcs:
        topos = wgs84.latlon(lat_given, lon_given, alt_iss) #We create a position in lat, lon, alt #we use the altitude of the ISS (distance ortho)
        carte = topos.itrs_xyz #We convert the position to a distance in the ECEF reference frame (ITRS)
        posi = ITRSPosition(carte) #we convert the distance to a position in ECEF(ITRS)
        carte_gcrs = posi.at(t1) #We convert this position ECEF into GCRS

        #distance:
        v = satellite.at(t1) - posi.at(t1)
        dist = v.distance().km
        a=float(dist)
        
        if a<int(dist_min):
            dist_min=a
            latmin = lat
            lonmin = lon
            altkmin = alt_iss/1000 #to have km
            h_min = hour
            min_min = minute
            sec_min = seconde
           
    return dist_min, latmin.degrees, lonmin.degrees, altkmin, h_min, min_min, sec_min


specific_day = [2023,12,22]
d, lati, loni, alti, hi, mi, si = shortest_distance(37.833959, -75.487833, specific_day)


#We redirect the outputs to a textfile, to take the answer to a c file after that
with open("output_coords.txt", "w") as f: #with closes the file alone
    print(d, lati, loni, alti, hi, mi, si, file=f)
