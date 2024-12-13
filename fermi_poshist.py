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

def read_eph_sgp8():

    file_name = '20240703_Fermi_sgp8.txt'
    with open(file_name) as f:
        lines = f.read().split('\n')

    lst_data = []
    for l in lines:

       lst_ = l.split()

       if len(lst_) != 7:
           continue

       lst_data.append((lst_[0], float(lst_[1]), float(lst_[3]), float(lst_[4]), float(lst_[5]), float(lst_[6])))
    

    tab = Table(rows=lst_data, names=('Date', 'Time_SOD', 'RA', 'Dec', 'R', 'dT'))
    #print(len(tab))
    #sys.exit()

    return tab

def test():

    r_GES_sgp4 = [6886.2084, 317.5682, -16.2330]
    r_GES_sgp8 = [6886.2000, 317.9091, -16.1338]
    r_GES_sc =   [6886.2575, 317.5577, -16.2300] 

    # SGP4
    # 317.5682066275038 -16.232959331342563 6886.208392270856 
    # [ 4879.95294521 -4460.97760221 -1924.99457319] - calc in skyfield
    # [ 4879.95142619 -4460.97724979 -1924.99926834] - calc


    r_DEC_sgp4 = GES_deg_to_DEC(r_GES_sgp4)
    r_DEC_sgp8 = GES_deg_to_DEC(r_GES_sgp8)
    r_DEC_sc = GES_deg_to_DEC(r_GES_sc)

    print(r_DEC_sgp4)

    dr = (sum((r_DEC_sgp8 - r_DEC_sc)**2))**0.5
    dr2 = (sum((r_DEC_sgp4 - r_DEC_sc)**2))**0.5
    dr3 = (sum((r_DEC_sgp4 - r_DEC_sgp8)**2))**0.5

    print(f'{dr:.3f} {dr2:.3f} {dr3:.3f}')

def main():

    poshist_file = './glg_poshist_all_240703_v00.fit'

    sc_time, sc_pos = read_pos(poshist_file)
    n_records = len(sc_time)

    tle_file = './20240703_Fermi_tle.txt'
    sc_name = 'Fermi'

    ts = load.timescale()

    tle_data = skyfield_tle.TLE(tle_file, sc_name)

    tab = read_eph_sgp8()
    
    str_out = ''

    ii = 0
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

        r_GES_sgp8 = zeros(3)
        r_GES_sgp8[0] = tab['R'][ii] # km
        r_GES_sgp8[1] = deg2rad(tab['RA'][ii])
        r_GES_sgp8[2] = deg2rad(tab['Dec'][ii])

        r_DEC_sgp8 = GES_to_DEC(r_GES_sgp8)

        dr_sgp8 = (sum((r_DEC_sgp8 - sc_pos[i,:])**2))**0.5
    
        sod = (time - time.replace(hour=0, minute=0, second=0, microsecond=0)).total_seconds()

        str_out += f'{str_t:20s} {sod:10.3f} {rad2deg(r_GES[1]):9.4f} {rad2deg(r_GES[2]):8.4f} {r_GES[0]:12.4f}' +\
                   f' {ra:9.4f} {dec:8.4f} {distance:12.4f} {days_from_epoch:5.2f} {dr:10.3f} '+\
                   f' {tab["Date"][ii]} {tab["Time_SOD"][ii]:10.3f} {tab["RA"][ii]:9.4f} {tab["Dec"][ii]:8.4f} {tab["R"][ii]:12.4f}'+\
                   f' {dr_sgp8:10.3f}\n'

        ii += 1

    print(str_out)

if __name__ == '__main__':

    main()
    #test()