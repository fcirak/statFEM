#!/usr/bin/env python

# script to extract the partial openFTL necessary to compile an application
# Jan 2012

# packages
import os, sys, string, shutil, getopt

# help message 
def helpMessage():
    print '''This script extract the partial openFTL necessary to compile an application

Options:
    -h,--help                      Displays this help message
    -q,--quiet                     Do not display verbose messages
    --srcopenftlroot=<directory>   Location of existing full openFTL tree
                                   from which extraction is done.
                                   If the option is not provided the enironment
                                   variable OPENFTLROOT is used.
    --dstopenftlroot=<directory>   New destination at which the partial
                                   openFTL tree is based
    --appcpp=<cppfile>             Location of application for which
                                   extraction is carried out
    --copyconfig=(0|1)             Copy openFTL make configuration files
                                   (default 0)'''
    sys.exit(0)

# set parameters
scriptname = 'bundle.py'
try:
    openftlroot = os.environ['OPENFTLROOT']
except KeyError:
    openftlroot = ''
dstopenftlroot = ''
appcpp = ''
copyconfig = 0
beQuiet = False

# read options
try:
    options, xarguments = getopt.getopt(sys.argv[1:],
    'hq', ['help', 'quiet', 'srcopenftlroot=', 'dstopenftlroot=', 'appcpp=', 'copyconfig='])
except getopt.error:
    print '''Error: You tried to use an unknown option or the
    argument for an option that requires it was missing. Try
    `'''+scriptname+''' -h\' for more information.'''
    sys.exit(0)
for a in options[:]:
    if (a[0] == '-h') or (a[0] == '--help'):
        helpMessage()
for a in options[:]:
    if (a[0] == '-q') or (a[0] == '--quiet'):
        beQuiet = True
for a in options[:]:
    if (a[0] == '--srcopenftlroot') and (a[1] != ''):
        openftlroot = a[1]
        options.remove(a)
        break
    elif (a[0] == '--srcopenftlroot') and (a[1] == ''):
        print '--srcopenftlroot expects an argument'
        sys.exit(0)
for a in options[:]:
    if (a[0] == '--dstopenftlroot') and (a[1] != ''):
        dstopenftlroot = a[1]
        options.remove(a)
        break
    elif (a[0] == '--dstopenftlroot') and (a[1] == ''):
        print '--dstopenftlroot expects an argument'
        sys.exit(0)
for a in options[:]:
    if (a[0] == '--appcpp') and (a[1] != ''):
        appcpp = a[1]
        options.remove(a)
        break
    elif (a[0] == '--appcpp') and (a[1] == ''):
        print '--appcpp expects an argument'
        sys.exit(0)
for a in options[:]:
    if (a[0] == '--copyconfig') and (a[1] != ''):
        copyconfig = int(a[1])
        options.remove(a)
        break
    elif (a[0] == '--copyconfig') and (a[1] == ''):
        print '--copyconfig expects an argument'
        sys.exit(0)

# sanity checks
openftlroot = os.path.expanduser(openftlroot)
if not os.path.exists(openftlroot):
    print 'Cannot find openftlroot \''+openftlroot+'\''
    sys.exit(-1)
appcpp = os.path.expanduser(appcpp)
if not os.path.isfile(appcpp):
    print 'Cannot find appcpp \''+appcpp+'\''
    sys.exit(-1)
dstopenftlroot = os.path.expanduser(dstopenftlroot)
approot = os.path.dirname(appcpp)
cppname = os.path.basename(appcpp)

# construct full file name
def completeFilename(name):
    if os.path.isfile(os.path.join(openftlroot,name)):
        return os.path.join(openftlroot,name)
    elif os.path.isfile(os.path.join(approot,name)):
        return os.path.join(approot,name)
    else:
        print 'Cannot locate file '+name
        sys.exit(-1)

# collect header files from 
def collectHeaders(dstHeaders, relFileName):
    if not beQuiet:
        print 'Browsing file '+relFileName
    absFileName = completeFilename(relFileName)
    fileHandle = open(absFileName, 'r')
    fileHeaders = set([])
    for line in fileHandle:
        if string.find(line, '#include') >= 0:
            angleOpen = string.find(line, '<')
            if angleOpen == -1:
                angleOpen = string.find(line, '"')
            angleClose = string.find(line, '>')
            if angleClose == -1:
                angleClose = string.rfind(line, '"')
            relName = line[angleOpen+1:angleClose]
            isHpp = (string.find(relName, '.hpp') >= 0)
            isIpp = (string.find(relName, '.ipp') >= 0)
            isBoost = (string.find(relName, 'boost/') >= 0)
            if isHpp and not isBoost:
                fileHeaders.add(relName)
            if isIpp:
                moddir = os.path.dirname(relFileName)
                fileHeaders.add(os.path.join(moddir,relName))
    fileHandle.close()
    fileHeaders -= dstHeaders  # only unread headers remain
    dstHeaders |= fileHeaders  # all headers so far
    for h in fileHeaders:  # read previously unread headers
        collectHeaders(dstHeaders, h)  # recursive call
    return

# collect all included OpenFTL headers
dstHeaders = set([])
collectHeaders(dstHeaders, cppname)
#print dstHeaders

# copy hpp and ipp headers to destination
for h in dstHeaders:
    hpp = os.path.join(openftlroot, h)
    dsthpp = os.path.join(dstopenftlroot, h)
    if not os.path.exists(hpp):
        if not beQuiet:
            print 'Skip '+h
        continue
    elif not beQuiet:
        print 'Copy '+hpp+'\n ==> '+dsthpp
    dstdir = os.path.dirname(dsthpp)
    if not os.path.exists(dstdir):
        os.makedirs(dstdir)
    shutil.copyfile(hpp, dsthpp)

# copy config items
if copyconfig:
    openftldef = 'openFTL_definitions.mk'
    defin = os.path.join(openftlroot, openftldef)
    dstdefin = os.path.join(dstopenftlroot, openftldef)
    if not beQuiet:
        print 'Copy '+defin+'\n ==> '+dstdefin
    shutil.copyfile(defin, dstdefin)
    openftlconfig = 'config'
    config = os.path.join(openftlroot, openftlconfig)
    dstconfig = os.path.join(dstopenftlroot, openftlconfig)
    if not beQuiet:
        print 'Copy '+config+'\n ==> '+dstconfig
    if os.path.exists(dstconfig):
        shutil.rmtree(dstconfig)
    shutil.copytree(config, dstconfig, ignore=shutil.ignore_patterns('.*'))

# quit
sys.exit(0)
