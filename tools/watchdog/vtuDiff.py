#!/usr/bin/env python3

# This script compares the DataArrays of 2 VTU files.
# Integer fields are directly compared. Floating point numbers
# are compared up to a provided tolerance.
#
# Please note:
#     The VTU files have to be in ASCII.
#     The DataArrays are related by their Name, thus these should be unique.
#
# Requirements:
#   - Adapted to Python 3.6.0, original was for Python 2.6.4

# packages
import os, sys, getopt
from xml.etree import ElementTree

# help message 
def help_message():
    print ('''vtuDiff.py compares the DataArrays of two VTU files. Integer fields are directly
compared. Floating point numbers are compared up to a provided tolerance.

Options: -h,--help          displays this help message
         --ref=<vtufile>    reference VTU file
         --new=<vtufile>    newly computed VTU file
         --comp=(abs|rel)   compare with absolute (default) or relative tolerance
         --tol=<double>     tolerance for comparison''')
    sys.exit(0)

# initialise constants
xmlFileRef = ""
xmlFileNew = ""
tolerance = 0.0
comparison = 'abs'

# read options
try:
    options, xarguments = getopt.getopt(sys.argv[1:],
    'h', ['help', 'ref=', 'new=', 'comp=', 'tol='])
except getopt.error:
    print ('''Error: You tried to use an unknown option or the
    argument for an option that requires it was missing. Try
    `vtuDiff.py -h\' for more information.''')
    sys.exit(0)
for a in options[:]:
    if (a[0] == '-h') or (a[0] == '--help'):
        help_message()
for a in options[:]:
    if (a[0] == '--ref') and (a[1] != ''):
        xmlFileRef = a[1]
        options.remove(a)
        break
    elif (a[0] == '--ref') and (a[1] == ''):
        print ('--ref expects an argument')
        sys.exit(0)
for a in options[:]:
    if (a[0] == '--new') and (a[1] != ''):
        xmlFileNew = a[1]
        options.remove(a)
        break
    elif (a[0] == '--new') and (a[1] == ''):
        print('--new expects an argument')
        sys.exit(0)
for a in options[:]:
    if (a[0] == '--comp') and (a[1] != ''):
        comparison = a[1]
        options.remove(a)
        break
    elif (a[0] == '--comp') and (a[1] == ''):
        print ('--comp expects an argument')
        sys.exit(0)
for a in options[:]:
    if (a[0] == '--tol') and (a[1] != ''):
        tolerance = float(a[1])
        options.remove(a)
        break
    elif (a[0] == '--tol') and (a[1] == ''):
        print ('--tol expects an argument')
        sys.exit(0)

# check existence of files etc
if (not os.path.isfile(xmlFileRef)):
    print ("Error: Reference file %s cannot be found" % xmlFileRef)
    sys.exit(-1)

if (not os.path.isfile(xmlFileNew)):
    print ("Error: New file %s cannot be found" % xmlFileNew)
    sys.exit(-1)

if (tolerance<=0.0):
    print ("Error: Tolerance %g must be positive" % tolerance)
    sys.exit(-1)


# parse reference VTU file
try:
    treeRef = ElementTree.parse(xmlFileRef)
except Exception(inst):
    print ("Unexpected error opening %s: %s" % (xmlFileRef, inst))
    sys.exit(-1)

# parse newly computed VTU file
try:
    treeNew = ElementTree.parse(xmlFileNew)
except Exception(inst):
    print ("Unexpected error opening %s: %s" % (xmlFileNew, inst))
    sys.exit(-1)

# collect all DataArrays of reference file
allArrays = {}
for node in treeRef.findall('.//DataArray'):
    name = node.attrib.get('Name')
    allArrays[name] = { 'type' : node.attrib.get('type'),
                        'data' : node.text.replace('\n',' ') }
#print allArrays

# loop all DataArrays in computed file and
# compare them with reference data
for nodeNew in treeNew.findall('.//DataArray'):
    nameNew = nodeNew.attrib.get('Name')
    try:
        array = allArrays[nameNew]
    except Exception:
        print ("Error: DataArray Name=\"%s\": Not existent in reference file" % (nameNew))
        sys.exit(-1)
    if (nodeNew.attrib.get('type') == array['type']):
        # integer data
        if (array['type'] == 'Int32'):
            intRef = array['data'].split()
            intNew = nodeNew.text.replace('\n',' ').split()
            if (len(intRef) != len(intNew)):
                print ("Error: DataArray Name=\"%s\": Length does not match" % (nameNew))
                sys.exit(-1)
            else:
                for i in range(len(intRef)):
                    if (int(intNew[i]) != int(intRef[i])):
                        print ("Error: DataArray Name=\"%s\": Entry %i does not match" % (nameNew, i))
                        sys.exit(-1)
        # floating point data
        elif (array['type'] == 'Float64'):
            floatRef = array['data'].split()
            floatNew = nodeNew.text.replace('\n',' ').split()
            if (len(floatRef) != len(floatNew)):
                print ("Error: DataArray Name=\"%s\": Length does not match" % (nameNew))
                sys.exit(-1)
            else:
                for i in range(len(floatRef)):
                    if (comparison == 'abs'):
                        if (abs(float(floatNew[i])-float(floatRef[i])) > tolerance):
                            print ("Error: DataArray Name=\"%s\": Entry %i does not match up to tolerance %g: New=%s, Ref=%s" % (nameNew, i, tolerance, floatNew[i], floatRef[i]))
                            sys.exit(-1)
                    elif (comparison == 'rel'):
                        # compare only larger entries
                        mx = max(abs(float(floatNew[i])),abs(float(floatRef[i])))
                        if ( mx > 1.0e-8 and abs(float(floatNew[i])-float(floatRef[i]))/mx > tolerance ):
                                print ("Error: DataArray Name=\"%s\": Entry %i does not match up to tolerance %g: New=%s, Ref=%s" % (nameNew, i, tolerance, floatNew[i], floatRef[i]))
                                sys.exit(-1)
                    else:
                        print ("Error: Cannot handle comparison=\"%s\"" % (comparison))
                        sys.exit(-1)
        # unknown number format
        else:
            print ("Error: Unknown number format")
            sys.exit(-1)

# successful exit
sys.exit(0)
        

