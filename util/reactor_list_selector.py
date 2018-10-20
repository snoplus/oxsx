######################################################################
# reactor_list_selector.py
# author: Jeff Lidgard <jeffrey.lidgard@physics.ox.ac.uk>
# Grabs information from the REACTORS and REACTORS_STATUS ratdb files
# for a specified reactor list and formats it into an output text file
#
# Revision History:
#  - 2018/10/19: first implementation
######################################################################

import os, sys, string, csv, argparse
from math import pow
import numpy as np

def getReactorList(reactorListName, filename):
    """Returns the list of reactors for the specified reactor list"""
    reactorListEntries = None
    with open(filename, 'r') as ratDBFile:
        for line in ratDBFile:
            if 'index: "'+reactorListName+'"' in line:
                nextlines = [next(ratDBFile) for x in xrange(12)]
                nextlines="".join(nextlines).rstrip()
                nextlines = nextlines.split("}")[0].rstrip()

                reactorListEntries = (nextlines.split('active_reactors: ['))[1].split('],')[0]
                reactorListEntries = (reactorListEntries.replace('"','')).split(',')
                reactorListEntries = map(str.strip, reactorListEntries)

    return reactorListEntries

def getReactorInfo(reactorListName, filename, filenameStatus, filenameOutput):
    """Returns the list of reactors for the specified reactor list"""
    # init variables
    reactorInfo = {}
    reactorStatusInfo = {}

    # get the reactors in the specified list
    reactorList = getReactorList(reactorListName, filename)

    # open both the REACTORS and the REACTORS_STATUS files
    with open(filename, 'r') as ratDBFile, open(filenameStatus, 'r') as ratDBStatusFile:

        # for each reactor in the list, lookup values
        for reactorName in reactorList:

            # go back to start of files, re-init vars
            ratDBFile.seek(0)
            ratDBStatusFile.seek(0)
            cores = None
            latitudes = None
            longitudes = None
            corePowers = None
            coreTypes = None
            nextlines = None
            nextlines2 = None

            # search REACTORS file
            for line in ratDBFile:
                if 'index: "'+reactorName+'"' in line:
                    nextlines = [next(ratDBFile) for x in xrange(12)]
                    nextlines="".join(nextlines).rstrip()
                    nextlines = nextlines.split("}")[0].rstrip()

            # search REACTORS_STATUS file
            for line in ratDBStatusFile:
                if 'index: "'+reactorName+'"' in line:
                    nextlines2 = [next(ratDBStatusFile) for x in xrange(12)]
                    nextlines2="".join(nextlines2).rstrip()
                    nextlines2 = nextlines2.split("}")[0].rstrip()

            #ensure the reactor info is found, otherwise raise an exception
            if nextlines is None:
                raise Exception("Didn't find reactor information in REACTORS file")
            if nextlines2 is None:
                raise Exception("Didn't find reactor information in REACTORS_STATUS file")

            # Now go through data and pull out figures
            # get number of cores
            cores = ((nextlines.split('no_cores:'))[1].split(',\n')[0]).strip()
            cores = int(cores)

            # get core power information
            corePowers = (nextlines2.split('core_power: ['))[1].split('],')[0]
            corePowers = corePowers.split(',')
            corePowers = map(str.strip, corePowers)
            corePowers = map(float, corePowers)
            corePower = int(round(sum(corePowers)/len(corePowers),0)) # rounded to nearest (and set as int)

            # get core type information
            coreTypes = (nextlines2.split('core_spectrum: ['))[1].split('],')[0]
            coreTypes = (coreTypes.replace('"','')).split(',')
            coreTypes = map(str.strip, coreTypes)
            coreType = coreTypes[0] # here is the assumption that all cores are the same type (they are in the current list)

            # get latitude information
            latitudes = (nextlines.split('latitude: ['))[1].split('],')[0]
            latitudes = latitudes.split(',')
            latitudes = map(str.strip, latitudes)
            latitudes = map(float, latitudes)
            latitude = np.average(latitudes, weights=corePowers) #average weighted for core power

            # get longitude information
            longitudes = (nextlines.split('longitude: ['))[1].split('],')[0]
            longitudes = longitudes.strip().split(',')
            longitudes = map(str.strip, longitudes)
            longitudes = map(float, longitudes)
            longitude = np.average(longitudes, weights=corePowers) #average weighted for core power

            # convert lat and long to distance
            distance = latLongToDistance(latitude,longitude)

            # write values to dictionary
            reactorInfo[reactorName] = [distance,coreType,cores,corePower]

    # write data to output file
    with open(filenameOutput,'wb') as fileOut:
        fileOut.write('reactorName,distance_km,spectrum_type,number_cores,average_core_power\n')
        for key, values in reactorInfo.items():
            fileOut.write( key+',' )
            values = map(str, values)
            values = ','.join(values)
            fileOut.write( values+'\n' )

    return None

def latLongToECEF(latitude, longitude, altitude):
    """Returns the distance (km) from a lat,long to SNOLAB"""
    # reference http://www.mathworks.co.uk/help/aeroblks/llatoecefposition.html
    #print latitude, longitude, altitude
    toRad = np.pi/180
    earthRadius = np.double(6378137.0) #Radius of the Earth (in meters)
    flattening = np.double(1./298.257223563) #Flattening factor WGS84 Model
    latitude = np.radians(latitude)
    longitude = np.radians(longitude)
    lfactor = np.arctan( pow((1. - flattening),2)*np.tan(latitude))*180./np.pi
    rs = np.sqrt( pow(earthRadius,2)/(1. + (1./pow((1. - flattening),2) - 1.)*pow(np.sin(np.radians(lfactor)),2)))
    x = (rs*np.cos(np.radians(lfactor))*np.cos(longitude) + altitude*np.cos(latitude)*np.cos(longitude))/1000 #in km
    y = (rs*np.cos(np.radians(lfactor))*np.sin(longitude) + altitude*np.cos(latitude)*np.sin(longitude))/1000 #in km
    z = (rs*np.sin(np.radians(lfactor)) + altitude*np.sin(latitude))/1000 #in km
    #print x,y,z
    return np.array([x,y,z])

def latLongToDistance(latitude, longitude, altitude=0):
    """Returns the distance (km) from a lat,long to SNOLAB"""
    SNOLLA  = np.array([-81.2014, 46.4753, -1766.0])
    SNOECEF = np.array([ 672.87, -4347.18, 4600.51]) # converted numbers
    ECEF = latLongToECEF(latitude, longitude, altitude)
    displacement = np.subtract(SNOECEF,ECEF)
    distance = np.linalg.norm(displacement)
    return round(distance,2)

if __name__=="__main__":
    parser = argparse.ArgumentParser("Pulls reactor info from ratdb files, output txt file contains selected reactor info.")
    parser.add_argument("-n", dest="reactorListName", required=True,
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
    args = parser.parse_args()

    getReactorInfo(args.reactorListName, args.REACTORS_filename, args.REACTORS_STATUS_filename, args.output_filename)
