import argparse
from collections import OrderedDict
import Parcellations
import os.path

#Find out if the fie exist
def is_valid_file(parser, arg):
    if not os.path.exists(arg):
        parser.error("The file %s does not exist!" % arg)
    else:
        return open(arg, 'r')  # return an open file handle


parser = argparse.ArgumentParser()
#Add valid arguments to take in
parser.add_argument('--inputFile', '-i', help='The .net which contains all the correlation Values', required=True)
parser.add_argument('--outputDir', '-o', help='The input directory where the dicoms live', required=True)
parser.add_argument('--fisherz', '-f', help='Ouput correlations though fisher-z transform', action="store_true")
parser.add_argument('--autoscale', help='Automatically find the scaling for the color used on graph', action="store_true")
parser.add_argument('--scaleBounds', help='Float value that will be used for both the max and the negative will be used as the minimum on the colorbar')
parser.add_argument('--outputUnordered', help='Output an unordered csv. the parcels will not be rearranged', action='store_true')
parser.add_argument('--outputNewOrder', help='Output a file with the ordering that was used to create the csv and the graph', action='store_true')
args = parser.parse_args()

input_file = args.inputFile
outputDir = args.outputDir
fisherz = args.fisherz
autoscale = args.autoscale
scaleBounds = args.scaleBounds
outputUnordered = args.outputUnordered
outputNewOrder = args.outputNewOrder


#input_file = '/home/mitchell/Desktop/sub-155938_power264.net'
print(input_file)
net_file = Parcellations.net_file(input_file)

parcellation = Parcellations.check_input(net_file)


network_dict = OrderedDict()
Parcellations.run(net_file=net_file, parcellation=parcellation, outputDir=outputDir, fisherz=fisherz
                 ,scaleBounds=scaleBounds, outputUnordered=False, outputNewOrder=outputNewOrder)
if outputUnordered:
    Parcellations.run(net_file=net_file, parcellation=parcellation, outputDir=outputDir, fisherz=fisherz
                      , scaleBounds=scaleBounds, outputUnordered=outputUnordered,outputNewOrder=outputNewOrder)
