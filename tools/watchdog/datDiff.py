#!/usr/bin/env python3

# This script compares the data tables of 2 ASCII files.
#
# Note:
#   - Comment lines are skipped. The default comment symbol is "#",
#     a different one can be optionally set.
#   - The entries in the table are separated by white spaces
#     and line breaks. All entries are converted to floating point
#     values. Floating point numbers are compared up to a provided tolerance.
#   - There can be several data tables of different lay-out.
#   - Empty lines are skipped.
#
# Example data file:
#  +-------------------------------+
#  |# Gnuplot data file            |
#  |# time force_x force_y         |
#  |  0.   0.      0.              |
#  |  0.5  0.1     0.              |
#  |  1.0  0.1234  3.141           |
#  |...............................|
#
# Requirements:
#   - Adapted to Python 3.6.0, original was for Python 2.6.4

# packages
import os, sys, getopt
import string

# help message 
def help_message():
    print ('''datDiff.py compares the data tables of two ASCII files.
Comment lines are skipped. All other entries are compared
as floating point numbers up to a provided tolerance.

Options: -h,--help          displays this help message
         --ref=<vtufile>    reference VTU file
         --new=<vtufile>    newly computed VTU file
         --comp=(abs|rel)   compare with absolute (default) or relative tolerance
         --tol=<double>     tolerance for comparison
         --csym=<char>      skip lines starting with comment symbol <char> (default \'#\')''')
    sys.exit(0)

# initialise constants
datFileRef = ""
datFileNew = ""
commSymbol = "#"
tolerance = 0.0
comparison = 'abs'

# read options
try:
    options, xarguments = getopt.getopt(sys.argv[1:],
    'h', ['help', 'ref=', 'new=', 'comp=', 'tol=', 'csym='])
except getopt.error:
    print ('''Error: You tried to use an unknown option or the
    argument for an option that requires it was missing. Try
    `datDiff.py -h\' for more information.''')
    sys.exit(0)
for a in options[:]:
    if (a[0] == '-h') or (a[0] == '--help'):
        help_message()
for a in options[:]:
    if (a[0] == '--ref') and (a[1] != ''):
        datFileRef = a[1]
        options.remove(a)
        break
    elif (a[0] == '--ref') and (a[1] == ''):
        print ('--ref expects an argument')
        sys.exit(0)
for a in options[:]:
    if (a[0] == '--new') and (a[1] != ''):
        datFileNew = a[1]
        options.remove(a)
        break
    elif (a[0] == '--new') and (a[1] == ''):
        print ('--new expects an argument')
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
for a in options[:]:
    if (a[0] == '--csym') and (a[1] != ''):
        commSymbol = a[1]
        options.remove(a)
        break
    elif (a[0] == '--csym') and (a[1] == ''):
        print ('--csym expects an argument')
        sys.exit(0)

# check existence of files etc
if (not os.path.isfile(datFileRef)):
    print ("Error: Reference file %s cannot be found" % datFileRef)
    sys.exit(-1)

if (not os.path.isfile(datFileNew)):
    print ("Error: New file %s cannot be found" % datFileNew)
    sys.exit(-1)

if (tolerance <= 0.0):
    print ("Error: Tolerance %g must be positive" % tolerance)
    sys.exit(-1)

# open files
refFile = open(datFileRef)
newFile = open(datFileNew)

# parse data files
lineNumber = 0
while (True):
    # read lines
    refLine = refFile.readline()
    newLine = newFile.readline()

    lineNumber += 1
    # check end-of-file
    if (not refLine and not newLine):
        break
    if (not refLine and newLine):
        print ("Error: Line %i: Reference file does not have this line" % lineNumber)
        sys.exit(-1)
    if (refLine and not newLine):
        print ("Error: Line %i: New file does not have this line" % lineNumber)
        sys.exit(-1)
    # strip leading and trailing white spaces of lines
    refLine = refLine.strip()
    newLine = newLine.strip()
    # skip lines made up only of white space
    if (refLine == '' and newLine == ''):
        continue
    if (refLine != '' and newLine == ''):
        print ("Error: Line %i: New file does not have white line" % lineNumber)
        sys.exit(-1)
    if (refLine == '' and newLine != ''):
        print ("Error: Line %i: Reference file does not have white line" % lineNumber)
        sys.exit(-1)
    # check for comment lines
    if (commSymbol != ""):
        refIsComment = (refLine.find(commSymbol) == 0)
        newIsComment = (newLine.find(commSymbol) == 0)
        if (refIsComment and newIsComment):
            continue
        if (not refIsComment and newIsComment):
            print ("Error: Line %i: New file is not commented there" % lineNumber)
            sys.exit(-1)
        if (refIsComment and not newIsComment):
            print ("Error: Line %i: Reference file is not commented there" % lineNumber)
            sys.exit(-1)
    # extract entries on lines and convert them to float
    floatRef = [float(e) for e in refLine.split()]
    floatNew = [float(e) for e in newLine.split()]
    if (len(floatRef) != len(floatNew)):
        print ("Error: Line %i: Different number of entries" % lineNumber)
        sys.exit(-1)
    # loop (column) entries
    for e in range(len(floatRef)):
        if (comparison == 'abs'):
            if (abs(float(floatNew[e])-float(floatRef[e])) > tolerance):
                print ("Error: Line %i: Entry %i does not match up to tolerance %g: New=%s, Ref=%s" \
                      % (lineNumber, e, tolerance, floatNew[e], floatRef[e]))
                sys.exit(-1)
        elif (comparison == 'rel'):
            # compare only larger entries
            mx = max(abs(float(floatNew[e])),abs(float(floatRef[e])))
            if ( mx > 1.0e-8 and abs(float(floatNew[e])-float(floatRef[e]))/mx > tolerance ):
                print ("Error: Line %i: Entry %i does not match up to tolerance %g: New=%s, Ref=%s" \
                      % (lineNumber, e, tolerance, floatNew[e], floatRef[e]))
                sys.exit(-1)
        else:
            print ("Error: Cannot handle comparison=\"%s\"" % comparison)
            sys.exit(-1)
    # passed everything so far -- increase line counter
    pass

# successful exit
sys.exit(0)
        

