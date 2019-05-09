# Create a class that represent a parcellation scheme
import csv
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from builtins import any
import os
from collections import OrderedDict


class net_file(object):

    def __init__(self, filename):
        self.filename = filename
        self.num_of_header_lines = 2
        self.get_subject()
        self.get_parcellation()

    def get_subject(self):
        file = os.path.basename(self.filename)
        subject = file.split('sub-')[1].split('_')[0]
        self.subject = subject

    def get_parcellation(self):
        file = os.path.basename(self.filename)
        parcellation = file.split('.')[0].split('_')[-1]
        self.parcellation = parcellation


class parcellation(object):
    def __init__(self, parcellation, network_dict, key_file, num_of_parcels, header_size=2):
        self.parcellation = parcellation
        self.network_dict = network_dict
        self.key_file = key_file
        self.header_size = header_size
        self.num_of_parcels = num_of_parcels
        self.get_line_num_dict()
        self.get_new_order()

    # Search through the Key to look for the search terms given in the network_dict
    # So this will return the a dictionary with networdks as the Keys and a list of parcels numbers that belong to that network
    def get_line_num_dict(self):

        line_num_dict = {}
        for key, value in network_dict.items():
            search_term = value
            line_list = []

            with open(self.key_file) as file:
                reader = csv.reader(file, delimiter=',')
                line_number = 0
                for line in reader:
                    line_number += 1
                    if any(search_term in x for x in line):
                        line_list.append(line_number - self.header_size)

            line_num_dict[key] = line_list

        self.line_num_dict = line_num_dict

    # Get the new ordering of the parcels so that we can make a cohesive graph.
    # This will return 2 values
    # a list of the new order of the parcels based on the order of the network dictionary
    # a list of the number of parcels in each network, this will be used to draw boundaries on the graph
    def get_new_order(self):

        new_order = []
        length_list = [0]
        for key, value in self.line_num_dict.items():
            length_list.append(len(self.line_num_dict[key]) + length_list[-1])
            new_order.extend(self.line_num_dict[key])
        length_list.pop(0)

        self.new_order = new_order
        self.length_list = length_list


#This will take the .net files and place the values of the correlations at the
# correct location according to the new ordering
def build_Matrix(parcellation, net_file, fisherz=False, unordered=False):

    #open the input file
    with open(net_file.filename) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=' ')
        skip_line = net_file.num_of_header_lines
        line_count = 0
        #Create a Zero matix the size of num of parcels x num of parcels
        a = np.zeros((parcellation.num_of_parcels, parcellation.num_of_parcels))


        for row in csv_reader:

            #Skip the headers
            if line_count < skip_line:
                line_count += 1
            else:
                if not unordered:
                    a[parcellation.new_order.index(int(row[0]) - 1), parcellation.new_order.index(
                        int(row[1]) - 1 )] = float(row[2])
                else:
                    a[int(row[0]) - 1, int(row[1]) - 1] = float(row[2])
                line_count += 1


    a = (a + a.T) / 2

    if fisherz:
        max = np.amax(a)
        min = np.amin(a)
        a = np.arctanh(a)

        np.fill_diagonal(a, np.inf)
    else:
        max = np.amax(a)
        min = np.amin(a)
        np.fill_diagonal(a, 1.)
    return a, max, min


def write_new_order(new_order,outfile):
    print("Writing the new ordering to file: " + outfile)
    with open(outfile, 'w') as f:
        for item in new_order:
            f.write("%s\n" % item)


#Plot the matrix to a corr matrix
def plot_correlation_matrix(parcellation, subject, output_file, matrix, max, min, fisherz=False):

    #Place in a data frame
    df = pd.DataFrame(matrix)
    plt.matshow(df)

    #Create the Title
    if fisherz:
        plt.title(str(subject) + ' ' + parcellation.parcellation + ' fisher-z')
    else:
        plt.title(str(subject) + ' ' + parcellation.parcellation)
    if abs(min) > abs(max):
        max = min * -1
    else:
        min = max * -1

    plt.imshow(df, vmin=min, vmax=max)

    #Setup the boundaries for the network
    for i in range(len(parcellation.length_list) - 1):
        plt.axhline(y=parcellation.length_list[i], color='black', linestyle='-')
        plt.axvline(x=parcellation.length_list[i], color='black', linestyle='-')
    names = list(parcellation.network_dict.keys())

    #set the tick location to be in the middle of the network boundaries
    tick_location = []
    tick_location.append((parcellation.length_list[0]) / 2)
    for i in range(1, len(parcellation.length_list)):
        tick_location.append(parcellation.length_list[i - 1] + ((parcellation.length_list[i] -
                                                                parcellation.length_list[i - 1]) / 2))

    #make the xticks slant
    plt.xticks(tick_location, names, rotation=45, horizontalalignment='right')
    ax = plt.gca()
    ax.tick_params(axis="x", direction="out", top=0, bottom=1, labelbottom=1, labeltop=0)
    plt.yticks(tick_location, names)

    plt.colorbar()
    plt.set_cmap('seismic')
    plt.savefig(output_file,bbox_inches='tight',pad_inches=.25, dpi=400)

#write the matrix to a csv
def write_out_correlation(matrix, out_file):
    print("Writing out matix: " + out_file)
    np.savetxt(out_file, matrix, delimiter=",")

#Check to see if we have setup that parcellation scheme
def check_input(net_file):
    if net_file.parcellation == 'aal116':
        print("Sorry " + net_file.parcellation + " has not yet been implemented")
        exit()
    elif net_file.parcellation == 'glasser360':
        print("Sorry " + net_file.parcellation + " has not yet been implemented")
        exit()
    elif net_file.parcellation == 'gordon333':
        print("Sorry " + net_file.parcellation + " has not yet been implemented")
        exit()
    elif net_file.parcellation == 'power264':
        print("Using " + net_file.parcellation)
        parcellation = Power264
    elif net_file.parcellation == 'schaefer100':
        print("Sorry " + net_file.parcellation + " has not yet been implemented")
        exit()
    elif net_file.parcellation == 'schaefer200':
        print("Sorry " + net_file.parcellation + " has not yet been implemented")
        exit()
    elif net_file.parcellation == 'schaefer400':
        print("Using " + net_file.parcellation)
        parcellation = Schaefer400
    else:
        print("Could not recognize parcellation: " + net_file.parcellation)
        exit()

    return parcellation


#Encapsulates all the processing
def run(net_file, parcellation, outputDir, fisherz, scaleBounds, outputUnordered, outputNewOrder):
    matrix, max, min = build_Matrix(parcellation=parcellation, net_file=net_file,
                                                  fisherz=False, unordered=outputUnordered)
    if (scaleBounds):
        max = float(scaleBounds)
        min = -1 * max

    filename = net_file.subject + '_' + net_file.parcellation

    if fisherz:
        filename = filename + '_fisherz'
    if outputUnordered:
        filename = filename + '_unordered'

    output_file = os.path.join(outputDir, filename + '.csv')
    write_out_correlation(matrix=matrix, out_file=output_file)

    if not outputUnordered:
        output_file = os.path.join(outputDir, filename + '.png')
        plot_correlation_matrix(parcellation=parcellation, output_file=output_file, subject=net_file.subject,
                                          matrix=matrix, max=max, min=min, fisherz=fisherz)
    if outputNewOrder:
        write_new_order(parcellation.new_order, os.path.join(outputDir, net_file.parcellation + '_reordering'))


#########################3Set up the different parcellation schemes
#The set up for the network dict is:
#"Label on plot":"what to look for in csv"
network_dict = {'Default': 'Default mode', 'Uncertain': 'Uncertain', 'Somatomotor': 'somatomotor', 'Cingulo': 'Cingulo-opercular Task Control',
                'Auditory': 'Auditory', 'Memory': 'Memory', 'Vent Attn': 'Ventral attention',
                'Fronto-parietal': 'Fronto-parietal', 'Salience': 'Salience', 'Subcortical': 'Subcortical',
                'Cerebellar': 'Cerebellar', 'DorsalAttn': 'Dorsal attention', 'Visual': 'Visual'}
key_name = '/usr/ParcellationKeys/Consensus264.csv'
Power264 = parcellation(parcellation='Power264', network_dict=network_dict, key_file=key_name,
                        num_of_parcels=264)


# The Network Dict is a network network name with the key value of what to look for in the key
network_dict = {'Default': 'Default', 'Visual': 'Vis', 'Cont': 'Cont', 'DorsalAttn': 'DorsAttn', 'VentAttn':
    'SalVentAttn', 'SomatoMotor': 'SomMot', 'Limbic': 'Limbic'}
#The file name of the key
key_name = '/usr/ParcellationKeys/schaeferKey400.csv'
Schaefer400 = parcellation(parcellation='Schaefer400', network_dict=network_dict, key_file=key_name,
                            num_of_parcels=400)
