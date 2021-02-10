# Scons build file
######################################################################
# User section - please adjust the lines here
######################################################################

# CLHEP setup - where CLHEP is installed

CLHEP_VERSION = '2.4.4.1'
CLHEPSYS='/home/sebouh/CLHEP_install'

# ROOT setup - where your ROOT is installed
ROOTSYS='/home/miguel/root'

######################################################################
# end of user input - do not change below this line
######################################################################

# Python imports 
import os;

# Setup build environment
env = Environment(CCFLAGS = '-O2', # -pg for gprof
                  LINKFLAGS = '-Wl' # -pg for gprof
                  )

# retrieve standard libraries and options from ROOT
env.MergeFlags('!' + ROOTSYS + '/bin/root-config --libs --glibs --cflags')
# add extra libraries for 3D Visualization
#env.MergeFlags('-lX3d -lRGL')
env.MergeFlags('-lRGL')

# setup CLHEP includes and libs
env.MergeFlags('-I'+CLHEPSYS+'/include/'
               + ' -L/usr/lib64'
               + ' -lCLHEP-Matrix-'+CLHEP_VERSION
               + ' -lCLHEP-Vector-'+CLHEP_VERSION)
# setup local needs
env.MergeFlags('-I#event -L#event -I#util -L#util -lAlignEvent -lUtilities')

# list other scons files (hierarchical setup)
#SConscript(['event/SConscript'], exports='env')
#SConscript(['util/SConscript'], exports='env')
SConscript(['simulation/SConscript'], exports='env')
SConscript(['kfa/SConscript'], exports='env')
SConscript(['validation/SConscript'], exports='env')
