'''
reactor_list_selector.py
author: Jeff Lidgard <jeffrey.lidgard@physics.ox.ac.uk>
Grabs information from the REACTORS and REACTORS_STATUS ratdb files
for a specified reactor list and formats it into an output text file

Revision History:
- 2018/10/19: first implementation

'''

import os
import sys
import argparse
import numpy as np
import pandas as pd
import ROOT
ROOT.gROOT.SetBatch(True) # ROOT in batch mode

def process(filename_input, filename_output):
    '''
    Returns the list of reactors for the specified reactor list
    '''
    file_input = ROOT.TFile.Open(filename_input);

    # check to see if file contains the tree, otherwise exit
    try:
        file_input.nt
    except AttributeError:
        raise Exception("Tree not found in input file")

    # get tree from root file
    tree = ROOT.gDirectory.Get('nt')

    # get vectors from tree
    myDict = {}
    for i_ev in xrange(tree.GetEntries()):

        #load entry
        tree.GetEntry(i_ev)

        # load branches, convert to numpy arrays
        mc_quench = np.asarray(tree.mc_quench)
        mc_neutrino_energy = np.asarray(tree.mc_neutrino_energy)
        mc_positron_energy = np.asarray(tree.mc_positron_energy)
        mc_neutron_energy = np.asarray(tree.mc_neutron_energy)
        mc_neutrino_position_r = np.asarray(tree.mc_neutrino_position_r)
        mc_neutrino_position_x = np.asarray(tree.mc_neutrino_position_x)
        mc_neutrino_position_y = np.asarray(tree.mc_neutrino_position_y)
        mc_neutrino_position_z = np.asarray(tree.mc_neutrino_position_z)
        mc_positron_position_r = np.asarray(tree.mc_positron_position_r)
        mc_positron_position_x = np.asarray(tree.mc_positron_position_x)
        mc_positron_position_y = np.asarray(tree.mc_positron_position_y)
        mc_positron_position_z = np.asarray(tree.mc_positron_position_z)
        mc_neutron_position_r = np.asarray(tree.mc_neutron_position_r)
        mc_neutron_position_x = np.asarray(tree.mc_neutron_position_x)
        mc_neutron_position_y = np.asarray(tree.mc_neutron_position_y)
        mc_neutron_position_z = np.asarray(tree.mc_neutron_position_z)
        ev_fit_validity = np.asarray(tree.ev_fit_validity)
        ev_fit_energy = np.asarray(tree.ev_fit_energy)
        ev_fit_position_r = np.asarray(tree.ev_fit_position_r)
        ev_fit_position_x = np.asarray(tree.ev_fit_position_x)
        ev_fit_position_y = np.asarray(tree.ev_fit_position_y)
        ev_fit_position_z = np.asarray(tree.ev_fit_position_z)
        ev_nhit = np.asarray(tree.ev_nhit)
        ev_time_days = np.asarray(tree.ev_time_days)
        ev_time_seconds = np.asarray(tree.ev_time_seconds)
        ev_time_nanoseconds = np.asarray(tree.ev_time_nanoseconds)
        reactor_info_latitude = np.asarray(tree.reactor_info_latitude)
        reactor_info_longitude = np.asarray(tree.reactor_info_longitude)
        reactor_info_altitude = np.asarray(tree.reactor_info_altitude)
        reactor_info_distance = np.asarray(tree.reactor_info_distance)

        # mulitindex dataframe labels (indicies)
        particles = ['neutrino','positron','neutron']
        #mc
        fields = ['energy','position_r','position_x','position_y','position_z']
        fields_l = fields*len(particles)
        particles_l=[]
        for i in xrange(len(particles)):
            particles_l += [particles[i]]*len(fields)
        df_mc_columns=[np.array(particles_l),np.array(fields_l)]
        #ev
        fields = ['ev_fit_validity', 'ev_fit_energy', 'ev_fit_position_r',\
            'ev_fit_position_x', 'ev_fit_position_y', 'ev_fit_position_z',\
            'ev_nhit', 'ev_time_days', 'ev_time_seconds', 'ev_time_nanoseconds']
        df_ev_columns=[np.array(fields)]
        #info
        fields = ['reactor_info_latitude', 'reactor_info_longitude',\
            'reactor_info_altitude', 'reactor_info_distance']
        df_info_columns=[np.array(fields)]

        # put list together into arrays
        mc_data_i = np.array([mc_neutrino_energy, mc_neutrino_position_r,\
            mc_neutrino_position_x, mc_neutrino_position_y, mc_neutrino_position_z,\
            mc_positron_energy, mc_positron_position_r,\
            mc_positron_position_x, mc_positron_position_y, mc_positron_position_z,\
            mc_neutron_energy, mc_neutron_position_r,\
            mc_neutron_position_x, mc_neutron_position_y, mc_neutron_position_z])

        ev_data_i = np.array([ev_fit_validity, ev_fit_energy, ev_fit_position_r,\
            ev_fit_position_x, ev_fit_position_y, ev_fit_position_z,\
            ev_nhit, ev_time_days, ev_time_seconds, ev_time_nanoseconds])

        info_data_i = np.array([reactor_info_latitude, reactor_info_longitude,\
            reactor_info_altitude, reactor_info_distance])

        # add new row to array for each event (using dataframe rows as events, columns as particle properties)
        if (i_ev==0):
            mc_data = mc_data_i
            ev_data = ev_data_i
            info_data = info_data_i
        else:
            mc_data = np.column_stack((mc_data, mc_data_i))
            ev_data = np.column_stack((ev_data, ev_data_i))
            info_data = np.column_stack((info_data, info_data_i))

        # create pandas dataframe for this event, append this to the total dataframe
        df_mc_i = pd.DataFrame(np.transpose(mc_data), columns=df_mc_columns)
        df_ev_i = pd.DataFrame(np.transpose(ev_data), columns=df_ev_columns)
        df_info_i = pd.DataFrame(np.transpose(info_data), columns=df_info_columns)
        
        for i_ev in [0,2]: 
            print df_ev_i
        
        





    # # save to msgpack
    # filename_output = os.path.splitext(filename_input)[0]+".msg"
    # df.to_msgpack(filename_output)

    # # save to pickle
    # filename_output = os.path.splitext(filename_input)[0]+".pkl"
    # df.to_pickle(filename_output)

    # # save to json
    # filename_output = os.path.splitext(filename_input)[0]+".json"
    # df.to_json(filename_output)

def main(args):
    '''
    main - pass args
    '''
    parser = argparse.ArgumentParser("Pulls reactor info from ratdb files, " \
        +"output txt file contains selected reactor info.")
    parser.add_argument("-i", dest='filename_input', type=str, nargs='?',
                        help='filename & path to input ROOT file',
                        default="/data/snoplus/lidgard/OXSX/flux1/bruce_flux1_day360.root")
    parser.add_argument("-o", dest='filename_output', type=str, nargs='?',
                        help='filename & path to output parquet file',
                        default="./output.parquet")
                        #default=os.path.splitext(args.filename_input)[0]+".parquet")
    args = parser.parse_args(args)

    # check if the specified files exist
    if os.path.isfile(args.filename_input):
        process(args.filename_input, args.filename_output)
    else:
        print "The specified input file cannot be found."

if __name__ == "__main__":
    main(sys.argv[1:])

        # myDict['event'+str(i_ev)] = {}

        # myDict['event'+str(i_ev)]

        # myDict['event'+str(i_ev)]['mc'] = {'mc_quench':mc_quench}
        # myDict['event'+str(i_ev)]['mc']['neutrino'] = {"energy":mc_neutrino_energy,
                                                # "position_r":mc_neutrino_position_r,
                                                # "position_x":mc_neutrino_position_x,
                                                # "position_y":mc_neutrino_position_y,
                                                # "position_z":mc_neutrino_position_z}
        # myDict['event'+str(i_ev)]['mc']['positron'] = {"energy":mc_positron_energy,
                                                # "position_r":mc_positron_position_r,
                                                # "position_x":mc_positron_position_x,
                                                # "position_y":mc_positron_position_y,
                                                # "position_z":mc_positron_position_z}
        # myDict['event'+str(i_ev)]['mc']['neutron'] = {"energy":mc_neutron_energy,
                                                # "position_r":mc_neutron_position_r,
                                                # "position_x":mc_neutron_position_x,
                                                # "position_y":mc_neutron_position_y,
                                                # "position_z":mc_neutron_position_z}
        # myDict['event'+str(i_ev)]['mc']['reactor_info'] = {"latitude":reactor_info_latitude,
                                                # "longitude":reactor_info_longitude,
                                                # "altitude":reactor_info_altitude,
                                                # "distance":reactor_info_distance}
        # myDict['event'+str(i_ev)]['ev'] = {"validity":ev_fit_validity,
                                        # "energy":ev_fit_energy,
                                        # "position_r":ev_fit_position_r,
                                        # "position_x":ev_fit_position_x,
                                        # "position_y":ev_fit_position_y,
                                        # "position_z":ev_fit_position_z,
                                        # "nhit":ev_nhit,
                                        # "time_days":ev_time_days,
                                        # "time_seconds":ev_time_seconds,
                                        # "time_nanoseconds":ev_time_nanoseconds}
