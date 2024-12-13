import sys

from skyfield.api import load, utc, EarthSatellite

class TLE:
    """
    TLE manipulation using EarthSatellite class from skyfield
    https://rhodesmill.org/skyfield/earth-satellites.html
    """

    def __init__(self, tle_file, sc_name):
        
        self.tle_file = tle_file
        self.lst_sat, self.lst_epoch = self.read_tle(tle_file, sc_name)
        

    def read_text_file(self, file_name):

        with open(file_name) as f:
            return f.read()

    def read_tle(self, file_name, sc_name):

        str_tle = self.read_text_file(file_name)
        lst_lines = [s.strip() for s in str_tle.split('\n') if len(s.strip()) > 0]
    
        lst_tle = [(l1, l2) for l1, l2 in zip(lst_lines[:-1:2], lst_lines[1::2])]
    
        lst_sat = []
        lst_epoch = []
        for tle in lst_tle:
            sat = EarthSatellite(tle[0], tle[1], sc_name)
            lst_sat.append(sat)
            lst_epoch.append(sat.epoch)
    
        return lst_sat, lst_epoch

    def get_nearest_time(self, time_ts, lst_ts):
        """
        https://stackoverflow.com/questions/9706041/
            finding-index-of-an-item-closest-to-the-value-in-a-list-thats-not-entirely-sort
        """

        idx, val = min(enumerate(lst_ts), key=lambda x: abs(x[1]-time_ts))
        return idx, val

    def get_pos(self, time_ts):
        """
        time_ts = ts.utc(2014, 1, 23, 11, 18, 7)
        """
        
        idx, t_epoch = self.get_nearest_time(time_ts, self.lst_epoch)

        days_from_epoch = time_ts - t_epoch
        #print('{:.3f} days away from epoch'.format(days_from_epoch))

        geoposition = self.lst_sat[idx].at(time_ts)

        return geoposition, days_from_epoch

    def get_geo_pos(self, time_ts):

        geoposition, days_from_epoch = self.get_pos(time_ts)

        cur_long = geoposition.subpoint().longitude.radians * 180 / np.pi
        cur_lat = geoposition.subpoint().latitude.radians * 180 / np.pi

        return cur_long, cur_lat, days_from_epoch

    def get_equat_pos(self, time_ts):

        geoposition, days_from_epoch = self.get_pos(time_ts)

        ra, dec, distance= geoposition.radec()

        r_DEC = geoposition.position.km

        #print(f'{type(ra)} {type(dec)} {type(distance)}')
        #print(f'{ra} {dec} {distance}')
        #print(f'{ra._degrees} {dec._degrees} {distance.km}')
        #sys.exit()

        return ra._degrees, dec._degrees, distance.km, r_DEC, days_from_epoch

def main():

    tle_file = './20240703_Fermi_tle.txt'
    sc_name = 'Fermi'

    ts = load.timescale()

    t = ts.utc(2024, 7, 3, 5, 24, 26)

    tle_data = TLE(tle_file, sc_name)
    ra, dec, distance, r_DEC, days_from_epoch = tle_data.get_equat_pos(t)

    str_t = t.utc_strftime('%Y-%m-%d %H:%M:%S.%f')
    print(f'{str_t} {ra:.4f} {dec:.4f} {distance:.4f}')

if __name__ == '__main__':
    main()