'''
######################################################################
# reactor_list_selector.py
# author: Jeff Lidgard <jeffrey.lidgard@physics.ox.ac.uk>
# Grabs information from the REACTORS and REACTORS_STATUS ratdb files
# for a specified reactor list and formats it into an output text file
#
# currently using kamland position
#
# Revision History:
#  - 2018/10/19: first implementation
#  - 2018/12/12: modification, reshaping
######################################################################
'''

import os
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
                next_lines = []
                next_line = next(ratdb_file)
                # read lines until a line containing a "}" is found
                while next_line.split("}")[0] is not "":
                    next_line = next(ratdb_file)
                    next_lines.append(next_line)
                next_lines = "".join(next_lines).rstrip()
                next_lines = next_lines.split("}")[0].rstrip()

                reactor_list_entries = (next_lines.split('active_reactors: ['))[1].split('],')[0]
                reactor_list_entries = (reactor_list_entries.replace('"', '')).split(',')
                reactor_list_entries = map(str.strip, reactor_list_entries)
                break
    return reactor_list_entries

def get_reactor_ratdb_info(reactor_ratdb_filename):
    '''
    Returns all the info from the reactors.ratdb file
    '''
    # init variables
    reactor_info = {}

    # open the REACTORS file
    with open(reactor_ratdb_filename, 'r') as ratdb_file:

        # search REACTORS file
        for line in ratdb_file:

            # re-init vars
            reactor_name = None
            cores = None
            latitudes = None
            longitudes = None
            next_lines = None

            # look for REACTOR entries
            if 'type: "REACTOR",' in line:
                next_lines = []
                next_line = next(ratdb_file)
                # read lines until a line containing a "}" is found
                while next_line.split("}")[0] is not "":
                    next_line = next(ratdb_file)
                    next_lines.append(next_line)
                next_lines = "".join(next_lines).rstrip()
                next_lines = next_lines.split("}")[0].rstrip()

                #ensure the reactor info is found, otherwise raise an exception
                if next_lines is None:
                    raise Exception("Didn't find reactor information in REACTORS file")

                # Now go through data and pull out figures
                # get reactor name (index)
                reactor_name = ((next_lines.split('index:'))[1].split(',\n')[0]).strip()
                reactor_name = reactor_name.replace('"', '')

                # get number of cores
                cores = ((next_lines.split('no_cores:'))[1].split(',\n')[0]).strip()
                cores = int(cores)

                # get latitude information
                latitudes = (next_lines.split('latitude: ['))[1].split('],')[0]
                latitudes = (latitudes.strip()).split(',')
                latitudes = map(str.strip, latitudes)
                latitudes = map(float, latitudes)

                # get longitude information
                longitudes = (next_lines.split('longitude: ['))[1].split('],')[0]
                longitudes = (longitudes.strip()).split(',')
                longitudes = map(str.strip, longitudes)
                longitudes = map(float, longitudes)

                # write values to dictionary
                reactor_info[reactor_name] = {}
                reactor_info[reactor_name]["latitudes"] = latitudes
                reactor_info[reactor_name]["longitudes"] = longitudes

    return reactor_info

def get_reactor_status_ratdb_info(reactor_ratdb_status_filename):
    '''
    Returns all the info from the reactors_status.ratdb file
    '''
    # init variables
    reactor_status_info = {}

    # open the REACTORS_STATUS file
    with open(reactor_ratdb_status_filename, 'r') as ratdb_status_file:

        # search REACTORS_STATUS file
        for line in ratdb_status_file:

            # re-init vars
            reactor_name = None
            cores = None
            core_powers = None
            core_powers_err = None
            core_types = None
            next_lines = None

            # look for REACTOR_STATUS entries
            if 'type: "REACTOR_STATUS",' in line:
                next_lines = []
                next_line = next(ratdb_status_file)
                # read lines until a line containing a "}" is found
                while next_line.split("}")[0] is not "":
                    next_line = next(ratdb_status_file)
                    next_lines.append(next_line)
                next_lines = "".join(next_lines).rstrip()
                next_lines = next_lines.split("}")[0].rstrip()

                # ensure the reactor info is found, otherwise raise an exception
                if next_lines is None:
                    raise Exception("Didn't find reactor information in REACTORS_STATUS file")

                # Now go through data and pull out figures
                # get reactor name (index)
                reactor_name = ((next_lines.split('index:'))[1].split(',\n')[0]).strip()
                reactor_name = reactor_name.replace('"', '')

                # get number of cores
                cores = ((next_lines.split('no_cores:'))[1].split(',\n')[0]).strip()
                cores = int(cores)

                # get core power information
                core_powers = (next_lines.split('core_power: ['))[1].split('],')[0]
                core_powers = core_powers.split(',')
                core_powers = map(str.strip, core_powers)
                core_powers = map(float, core_powers)
                
                # get core power error information
                core_powers_err = (next_lines.split('core_power_err: ['))[1].split('],')[0]
                core_powers_err = core_powers_err.split(',')
                core_powers_err = map(str.strip, core_powers_err)
                core_powers_err = map(float, core_powers_err)

                # get core type information
                core_types = (next_lines.split('core_spectrum: ['))[1].split('],')[0]
                core_types = (core_types.replace('"', '')).split(',')
                core_types = map(str.strip, core_types)

                # write values to dictionary
                reactor_status_info[reactor_name] = {}
                reactor_status_info[reactor_name]["cores"] = cores
                reactor_status_info[reactor_name]["core_powers"] = core_powers
                reactor_status_info[reactor_name]["core_powers_err"] = core_powers_err
                reactor_status_info[reactor_name]["core_types"] = core_types

    return reactor_status_info

def get_reactor_info(reactor_ratdb_info, reactor_status_ratdb_info, reactor_list="All", position = "sno"):
    '''
    Returns the combined, averaged info from the reactors.ratdb and reactors_status.ratdb files
    '''
    # init variables
    reactor_info = {}

    for reactor_name in reactor_ratdb_info:

        #ensure the reactor info is found, otherwise raise an exception
        if reactor_name is None:
            raise Exception("Reactor is None, something is wrong")
        try:
            reactor_status_ratdb_info[reactor_name]
        except KeyError:
            raise Exception("Didn't find corresponding reactor status information")

        # combine dictionaries
        reactor_info[reactor_name] = reactor_ratdb_info[reactor_name]
        reactor_info[reactor_name]["cores"] = reactor_status_ratdb_info[reactor_name]["cores"]
        reactor_info[reactor_name]["core_powers"] = \
            reactor_status_ratdb_info[reactor_name]["core_powers"]
        reactor_info[reactor_name]["core_powers_err"] = \
            reactor_status_ratdb_info[reactor_name]["core_powers_err"]
        reactor_info[reactor_name]["core_types"] = \
            reactor_status_ratdb_info[reactor_name]["core_types"]

        # Now go through data and add some average values
        # rounded core power (and set as int)
        reactor_info[reactor_name]["core_power"] = \
            int(round(
                sum(reactor_info[reactor_name]["core_powers"])\
                            /len(reactor_info[reactor_name]["core_powers"])
                , 0))
        reactor_info[reactor_name]["core_power_err"] = \
            int(round(
                sum(reactor_info[reactor_name]["core_powers_err"])\
                            /len(reactor_info[reactor_name]["core_powers_err"])
                , 0))

        # get latitude information
        try:
            #average weighted for core power
            latitude = np.average(reactor_info[reactor_name]["latitudes"],\
                weights=reactor_info[reactor_name]["core_powers"])
        except ZeroDivisionError:
            #and if it fails for some reason then don't weight
            latitude = np.average(reactor_info[reactor_name]["latitudes"])
        reactor_info[reactor_name]["latitude"] = latitude

        # get longitude information
        try:
            #average weighted for core power
            longitude = np.average(reactor_info[reactor_name]["longitudes"],\
                weights=reactor_info[reactor_name]["core_powers"])
        except ZeroDivisionError:
            #and if it fails for some reason then don't weight
            longitude = np.average(reactor_info[reactor_name]["longitudes"])
        reactor_info[reactor_name]["longitude"] = longitude

        # convert lat and long to distance
        reactor_info[reactor_name]["distance"] = lat_long_to_distance(position, latitude, longitude)

        # get core type information
        # assumption that all cores are the same type
        core_type = max(set(reactor_info[reactor_name]["core_types"]),\
            key=reactor_info[reactor_name]["core_types"].count)
        # if assumption not true, throw error
        if any(i != core_type for i in reactor_info[reactor_name]["core_types"]):
            print "Warning: mix of core types: "\
                +reactor_name+", core_type:"+core_type+", All cores:"\
                +str(reactor_info[reactor_name]["core_types"])\
                +", distance:"+str(reactor_info[reactor_name]["distance"])
        reactor_info[reactor_name]["core_type"] = core_type

    if reactor_list == "All":
        return reactor_info
    else:
        #Return the reactor info for the specified reactor list
        #is this broken? (i always seem to be using 'all')
        reactor_info_select = {}
        for reactor_name in reactor_list: # select out reactor from the list
            reactor_info_select[reactor_name] = reactor_info[reactor_name]
        return reactor_info_select

def write_output_file(filename_output, reactor_info):
    '''
    write output file
    '''
    # write data to output file
    if not os.path.isfile(filename_output):
        write_header = True
    else:
        write_header = False
    with open(filename_output, 'ab+') as file_out:
        if write_header:
            file_out.write(
                'reactor_name,distance_km,spectrum_type,number_cores,average_core_power,core_power_stderr\n')
        for key, values in reactor_info.items():
            file_out.write(key\
                +','+str(values["distance"])\
                +','+str(values["core_type"])\
                +','+str(values["cores"])\
                +','+str(values["core_power"])\
                +','+str(values["core_power_err"])\
                +'\n')

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

def lat_long_to_distance(position, latitude, longitude, altitude=0):
    '''
    Returns the distance (km) from a lat,long to SNOLAB
    '''
    #SNOLLA  = np.array([-81.2014, 46.4753, -1766.0])
    kamland_ecef = np.array([-3777.14425893, 3483.58137383, 3766.0181443]) #using kamland (lat, long, alt) = (36.4225, 137.3153, 0.358)
    sno_ecef = np.array([672.87, -4347.18, 4600.51]) # converted numbers
    ecef = lat_long_to_ecef(latitude, longitude, altitude)
    if position == "sno":
        displacement = np.subtract(sno_ecef, ecef)
    if position == "kamland":
        displacement = np.subtract(kamland_ecef, ecef)
    distance = np.linalg.norm(displacement)
    return round(distance, 2)

def main(args):
    '''
    main - pass args
    '''
    #print args
    parser = argparse.ArgumentParser("Pulls reactor info from ratdb files, " \
        +"output txt file contains selected reactor info.")
    parser.add_argument("-n", dest="reactor_list_name",
                        help="reactor list (from REACTORS.ratdb) to use, default is 'All'",\
                        default="All")
    parser.add_argument("-i", dest='REACTORS_filename', type=str, nargs='?',
                        help='filename & path to REACTORS.ratdb',
                        default="/home/lidgard/rat-7/data/REACTORS.ratdb")
    parser.add_argument("-s", dest='REACTORS_STATUS_filename', type=str, nargs='?',
                        help='filename & path to REACTORS_STATUS.ratdb',
                        default="/home/lidgard/rat-7/data/REACTORS_STATUS.ratdb")
    parser.add_argument("-o", dest='output_filename', type=str, nargs='?',
                        help='filename & path to output file.ratdb',
                        default="reactor_list_selection.csv")
    parser.add_argument("-p", dest='position', type=str, nargs='?',
                        help='location to take distance to, options: \'sno\' or \'kamland\'',
                        default="sno")
    parser.add_argument('--fearless',action='store_true',\
                        help='Do not prompt the user for input - run without fear!')
    args = parser.parse_args(args)

    # check if the specified files exist
    if os.path.isfile(args.REACTORS_filename) and os.path.isfile(args.REACTORS_STATUS_filename):
        if os.path.isfile(args.output_filename) and not args.fearless:
            exists_query = str(raw_input("File exists, OK to append? 'y' to append," \
                " any other key to exit..")).lower().strip()
            if exists_query != "y":
                print "Exiting..."
                sys.exit()
        print "Getting reactor info..."

        # get the reactors in the specified list
        if args.reactor_list_name is not "All":
            reactor_list_name = get_reactor_list(args.reactor_list_name, args.REACTORS_filename)
        else:
            reactor_list_name = args.reactor_list_name

        # then get info
        reactor_ratdb_info = get_reactor_ratdb_info(args.REACTORS_filename)
        reactor_status_ratdb_info = get_reactor_status_ratdb_info(args.REACTORS_STATUS_filename)

        reactor_info = get_reactor_info(reactor_ratdb_info, \
            reactor_status_ratdb_info,\
            reactor_list_name, \
            position)
        write_output_file(args.output_filename, reactor_info)
    else:
        print "One of the specified ratdb files cannot be found, check paths. Exiting..."

if __name__ == "__main__":
    main(sys.argv[1:])
