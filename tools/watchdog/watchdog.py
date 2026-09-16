# -*- using: utf-8 -*-

import os
import sys
import subprocess
import datetime


def watchdog(statFEMRoot, ignoreError = False):
    '''
    statFEMRoot: The path to statFEM directory set using $STATFEMROOT. 
                 When using Circle continuous integration, $STATFEMROOT
                 is '/root/project' (see Circle's file `config.yml`).
    ignoreError: False - default: to let watchdog stop when an error occurs.
                 True  - to let watchdog continue until all checks are done. 
                         Failed application outputs red-coloured [Fail].
    '''

    print('\033[0m' + " -- TEST START --")
    
    # check on the file list "appList"
    watchdogDir = statFEMRoot + "tools/watchdog/"
    listDir = watchdogDir + "appList.txt"
    
    # define name of the temporary file
    suffix = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    tmpLogFileDir = watchdogDir + suffix + ".tmp"
    
    # clean dormant temporary files
    cleanTmpCommand = "rm -f " + watchdogDir + "*.tmp"
    subprocess.call(cleanTmpCommand, shell = True)

    # loop over the application list    
    with open(listDir, "r") as appList:
        for app in appList:
            # locate the app
            app = app.replace("\n", "") # remove the "\n" at the end of the line
            dir = statFEMRoot + app

            # run watchdog
            sys.stdout.write('\033[0m' + " - " + app)
            sys.stdout.flush()

            # save the terminal message to temporary log file
            testCommand = "cd " + os.path.dirname(dir) + "; STATFEMROOT=" + statFEMRoot + \
                          " make -f ./" + os.path.basename(dir) + " > " + tmpLogFileDir + " 2>&1"
            subprocess.call(testCommand, shell = True)

            # read log file to see if there is any problem
            cleanPass = compilePass = runPass = comparisonPass = False
            with open(tmpLogFileDir, "r") as logFile:
                for info in logFile:
                    if "WatchDog clean successful" in info:
                        cleanPass = True
                    if "WatchDog build successful" in info:
                        compilePass = True
                    if "WatchDog run successful"   in info:
                        runPass = True
                    if "WatchDog test successful"  in info:
                        comparisonPass = True

            # if one of the test fail
            if not cleanPass or not compilePass or not runPass or not comparisonPass:
                tmpStream  = open(tmpLogFileDir, 'r')
                tmpContent = tmpStream.read()
                sys.stdout.write(tmpContent)
            
            # test passes
            if cleanPass and compilePass and runPass and comparisonPass: 
                # remove the tempary file       
                subprocess.call(cleanTmpCommand, shell = True)
                sys.stdout.write('\033[92m' + " [Pass]\n")
            # test fails    
            else:
                sys.stdout.write('\033[91m' + " [Fail]\n")
                if not cleanPass:   
                    sys.stdout.write("   ** ERROR: Cleaning fails **\n")
                elif not compilePass:
                    sys.stdout.write("   ** ERROR: Compiling fails **\n")
                elif not runPass:
                    sys.stdout.write("   ** ERROR: Running fails **\n")
                elif not comparisonPass:
                    sys.stdout.write("   ** ERROR: Comparison with the benchmark fails **\n")
                
                # remove the tempary file
                subprocess.call(cleanTmpCommand, shell = True)
                if not ignoreError:
                    print('\033[0m' + " -- TEST FAILS --") 
                    sys.exit(1)
                    
    print('\033[0m' + " -- TEST FINISH --")        
    sys.exit(0)
  
  
if __name__ == "__main__":
    if len(sys.argv) < 2:
        watchdog(os.environ.get("STATFEMROOT")+"/")
    else:
        watchdog(os.environ.get("STATFEMROOT")+"/", True)
