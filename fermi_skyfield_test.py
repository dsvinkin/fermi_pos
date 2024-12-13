import sys
import datetime

from numpy import pi, sin, cos, arcsin, arccos, arctan2, deg2rad, rad2deg, size, zeros, sum

from astropy.io import fits
from astropy.table import Table

from skyfield.api import load, utc

import clock

import skyfield_tle

def DEC_to_GES(r_DEC):
    """
    from DEC to GES
    """
    r = (r_DEC[0]*r_DEC[0] + r_DEC[1]*r_DEC[1] + r_DEC[2]*r_DEC[2])**0.5

    alpha = arctan2(r_DEC[1],r_DEC[0])
    if (alpha < 0):
        alpha += 2*pi
    
    delta = arcsin(r_DEC[2]/r)

    r_GES = zeros(3)
    r_GES[0]=r
    r_GES[1]=alpha
    r_GES[2]=delta

    return r_GES

def GES_to_DEC(r_GES):
    """
    from GES to DEC
    """

    r = r_GES[0]
    alpha = r_GES[1]
    delta = r_GES[2]

    r_DEC = zeros(3)
    r_DEC[0] = r*cos(alpha)*cos(delta)
    r_DEC[1] = r*sin(alpha)*cos(delta)
    r_DEC[2] = r*sin(delta)

    return r_DEC

def GES_deg_to_DEC(r_GES):
    """
    from GES to DEC
    """

    r = r_GES[0]
    alpha = deg2rad(r_GES[1])
    delta = deg2rad(r_GES[2])

    r_DEC = zeros(3)
    r_DEC[0] = r*cos(alpha)*cos(delta)
    r_DEC[1] = r*sin(alpha)*cos(delta)
    r_DEC[2] = r*sin(delta)

    return r_DEC

def read_pos(poshist_file):

    data = fits.getdata(poshist_file, ext=1)

    sc_pos = zeros((size(data),3), float)

    sc_time = data.SCLK_UTC

    sc_pos[:,0] = data.POS_X / 1e3 # km
    sc_pos[:,1] = data.POS_Y / 1e3 # km
    sc_pos[:,2] = data.POS_Z / 1e3 # km

    return sc_time, sc_pos

def main():

    poshist_file = './glg_poshist_all_240703_v00.fit'

    sc_time, sc_pos = read_pos(poshist_file)
    n_records = len(sc_time)

    tle_file = './20240703_Fermi_tle.txt'
    sc_name = 'Fermi'

    ts = load.timescale()

    tle_data = skyfield_tle.TLE(tle_file, sc_name)
    
    str_out = 'Time SOD  RA_poshist Dec_poshist R_poshist RA_sf Dec_sf R_sf  Delta_epoch Delta_R'

    for i in range(0, n_records, 200):

        time = clock.fermi2utc(sc_time[i])
        time = time.replace(tzinfo=datetime.timezone.utc)
        str_t = time.strftime('%Y-%m-%d %H:%M:%S')
      
   
        t = ts.from_datetime(time)
        ra, dec, distance, r_DEC_sf, days_from_epoch = tle_data.get_equat_pos(t)

        #print(ra, dec, distance, r_DEC_sf)
        #sys.exit()

        dr = (sum((r_DEC_sf - sc_pos[i,:])**2))**0.5

        r_GES = DEC_to_GES(sc_pos[i,:])
    
        sod = (time - time.replace(hour=0, minute=0, second=0, microsecond=0)).total_seconds()

        str_out += f'{str_t:20s} {sod:10.3f} {rad2deg(r_GES[1]):9.4f} {rad2deg(r_GES[2]):8.4f} {r_GES[0]:12.4f}' +\
                   f' {ra:9.4f} {dec:8.4f} {distance:12.4f} {days_from_epoch:5.2f} {dr:10.3f}\n'


    with open('pos_comp.txt', 'w') as f:
        f.write(str_out)

if __name__ == '__main__':

    main()