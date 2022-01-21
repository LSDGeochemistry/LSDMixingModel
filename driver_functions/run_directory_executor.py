#Loop through and run all the folders in a file containing multiple runs
import subprocess
import os

#Set the top folder for the model to run through
root = '/exports/csce/datastore/geos/users/s0933963/github/LSDMixingModel/Runs'
#Assumes that the directory is in the runs folder
# DataDirectory = root + '/changing_base_level/0_1mm_linear'
DataDirectory = root + '/Model_testing'

#For each folder run the mixing model, simple stuff that saves some time
for subdirs, dirs, files in os.walk(DataDirectory):
    for dirs in dirs:
        temp_name = './mixing_column_variable_erosion_new_elev_bc.out  '+DataDirectory+'/'+str(dirs)
        print(temp_name)
        subprocess.call(temp_name,shell=True)


