'''
######################################################################
# reactor_list_selector.py
# author: Jeff Lidgard <jeffrey.lidgard@physics.ox.ac.uk>
# Grabs information from the REACTORS and REACTORS_STATUS ratdb files
# for a specified reactor list and formats it into an output text file
#
# Revision History:
#  - 2018/10/19: first implementation
######################################################################
'''

import os
import os.path
import sys
import argparse
import numpy as np

def get_reactor_list(reactor_list_name, filename):
    '''
    Returns the list of reactors for the specified reactor list
    '''
    reactor_list_entries = None
    with open(filename, 'r') as ratdb_file:
        for line in ratdb_file:
            if 'index: "'+reactor_list_name+'"' in line:
                next_lines = [next(ratdb_file) for x in xrange(12)]
                next_lines = "".join(next_lines).rstrip()
                next_lines = next_lines.split("}")[0].rstrip()

                reactor_list_entries = (next_lines.split('active_reactors: ['))[1].split('],')[0]
                reactor_list_entries = (reactor_list_entries.replace('"', '')).split(',')
                reactor_list_entries = map(str.strip, reactor_list_entries)
                break
    return reactor_list_entries

def get_reactor_info(reactor_list_name, filename, filename_status, filename_output):
    '''
    Returns the list of reactors for the specified reactor list
    '''
    # init variables
    reactor_info = {}

    # get the reactors in the specified list
    reactor_list = get_reactor_list(reactor_list_name, filename)

    # open both the REACTORS and the REACTORS_STATUS files
    with open(filename, 'r') as ratdb_file, open(filename_status, 'r') as ratdb_status_file:

        # for each reactor in the list, lookup values
        for reactor_name in reactor_list:

            # go back to start of files, re-init vars
            ratdb_file.seek(0)
            ratdb_status_file.seek(0)
            cores = None
            latitudes = None
            longitudes = None
            core_powers = None
            core_types = None
            next_lines = None
            next_lines2 = None

            # search REACTORS file
            for line in ratdb_file:
                if 'index: "'+reactor_name+'"' in line:
                    next_lines = [next(ratdb_file) for x in xrange(12)]
                    next_lines = "".join(next_lines).rstrip()
                    next_lines = next_lines.split("}")[0].rstrip()
                    break

            # search REACTORS_STATUS file
            for line in ratdb_status_file:
                if 'index: "'+reactor_name+'"' in line:
                    next_lines2 = [next(ratdb_status_file) for x in xrange(12)]
                    next_lines2 = "".join(next_lines2).rstrip()
                    next_lines2 = next_lines2.split("}")[0].rstrip()
                    break

            #ensure the reactor info is found, otherwise raise an exception
            if next_lines is None:
                raise Exception("Didn't find reactor information in REACTORS file")
            if next_lines2 is None:
                raise Exception("Didn't find reactor information in REACTORS_STATUS file")

            # Now go through data and pull out figures
            # get number of cores
            cores = ((next_lines.split('no_cores:'))[1].split(',\n')[0]).strip()
            cores = int(cores)

            # get core power information
            core_powers = (next_lines2.split('core_power: ['))[1].split('],')[0]
            core_powers = core_powers.split(',')
            core_powers = map(str.strip, core_powers)
            core_powers = map(float, core_powers)
            # rounded to nearest (and set as int)
            core_power = int(round(sum(core_powers)/len(core_powers), 0))

            # get core type information
            core_types = (next_lines2.split('core_spectrum: ['))[1].split('],')[0]
            core_types = (core_types.replace('"', '')).split(',')
            core_types = map(str.strip, core_types)
            core_type = core_types[0] # assumption that all cores are the same type

            # get latitude information
            latitudes = (next_lines.split('latitude: ['))[1].split('],')[0]
            latitudes = (latitudes.strip()).split(',')
            latitudes = map(str.strip, latitudes)
            latitudes = map(float, latitudes)
            try:
                #average weighted for core power
                latitude = np.average(latitudes, weights=core_powers)
            except ZeroDivisionError:
                #and if it fails for some reason then don't weight
                latitude = np.average(latitudes)

            # get longitude information
            longitudes = (next_lines.split('longitude: ['))[1].split('],')[0]
            longitudes = (longitudes.strip()).split(',')
            longitudes = map(str.strip, longitudes)
            longitudes = map(float, longitudes)
            try:
                #average weighted for core power
                longitude = np.average(longitudes, weights=core_powers)
            except ZeroDivisionError:
                #and if it fails for some reason then don't weight
                longitude = np.average(longitudes)

            # convert lat and long to distance
            distance = lat_long_to_distance(latitude, longitude)

            # write values to dictionary
            reactor_info[reactor_name] = [distance, core_type, cores, core_power]

    # write data to output file
    if not os.path.isfile(filename_output):
        write_header = True
    else:
        write_header = False
    with open(filename_output, 'ab+') as file_out:
        if write_header:
            file_out.write\
                ('reactor_name,distance_km,spectrum_type,number_cores,average_core_power\n')
        for key, values in reactor_info.items():
            file_out.write(key+',')
            values = map(str, values)
            values = ','.join(values)
            file_out.write(values+'\n')

    return None

def lat_long_to_ecef(latitude, longitude, altitude):
    '''
    Returns the distance (km) from a lat,long to SNOLAB
    '''
    # reference http://www.mathworks.co.uk/help/aeroblks/llatoecefposition.html
    #print latitude, longitude, altitude
    to_rad = np.pi/180
    earth_radius = 6378137.0 #Radius of the Earth (in meters)
    flattening = 1./298.257223563 #Flattening factor WGS84 Model
    latitude = np.radians(latitude)
    longitude = np.radians(longitude)
    lfactor = np.arctan(pow((1. - flattening), 2)*np.tan(latitude))/to_rad
    rs_factor = np.sqrt(pow(earth_radius, 2)/(1. + (1./pow((1. - flattening), 2) - 1.) \
        *pow(np.sin(np.radians(lfactor)), 2)))
    x_coord = (rs_factor*np.cos(np.radians(lfactor))*np.cos(longitude) \
        + altitude*np.cos(latitude)*np.cos(longitude))/1000 #in km
    y_coord = (rs_factor*np.cos(np.radians(lfactor))*np.sin(longitude) \
        + altitude*np.cos(latitude)*np.sin(longitude))/1000 #in km
    z_coord = (rs_factor*np.sin(np.radians(lfactor)) + altitude*np.sin(latitude))/1000 #in km
    #print x,y,z
    return np.array([x_coord, y_coord, z_coord])

def lat_long_to_distance(latitude, longitude, altitude=0):
    '''
    Returns the distance (km) from a lat,long to SNOLAB
    '''
    #SNOLLA  = np.array([-81.2014, 46.4753, -1766.0])
    sno_ecef = np.array([672.87, -4347.18, 4600.51]) # converted numbers
    ecef = lat_long_to_ecef(latitude, longitude, altitude)
    displacement = np.subtract(sno_ecef, ecef)
    distance = np.linalg.norm(displacement)
    return round(distance, 2)

def main(args):
    '''
    main - pass args
    '''
    parser = argparse.ArgumentParser("Pulls reactor info from ratdb files, " \
        +"output txt file contains selected reactor info.")
    parser.add_argument("-n", dest="reactor_list_name", required=True,
                        help="reactor list (from REACTORS.ratdb) to use")
    parser.add_argument("-i", dest='REACTORS_filename', type=str, nargs='?',
                        help='filename & path to REACTORS.ratdb',
                        default="/home/lidgard/rat3/data/REACTORS.ratdb")
    parser.add_argument("-s", dest='REACTORS_STATUS_filename', type=str, nargs='?',
                        help='filename & path to REACTORS_STATUS.ratdb',
                        default="/home/lidgard/rat3/data/REACTORS_STATUS.ratdb")
    parser.add_argument("-o", dest='output_filename', type=str, nargs='?',
                        help='filename & path to output file.ratdb',
                        default="reactor_list_selection.csv")
    args = parser.parse_args(args)

    # check if the specified files exist
    if os.path.isfile(args.REACTORS_filename) and os.path.isfile(args.REACTORS_STATUS_filename):
        if os.path.isfile(args.output_filename):
            exists_query = str(raw_input("File exists, OK to append? 'y' to append," \
                " any other key to exit..")).lower().strip()
            if exists_query != "y":
                print "Exiting..."
                sys.exit()
        print "Getting reactor info..."
        get_reactor_info(args.reactor_list_name, args.REACTORS_filename, \
            args.REACTORS_STATUS_filename, args.output_filename)
    else:
        print "One of the specified ratdb files cannot be found, check paths. Exiting..."

if __name__ == "__main__":
    main(sys.argv[1:])
