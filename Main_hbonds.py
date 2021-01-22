#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@Author: Andres Vodopivec
@Date: 2020-12-08
@Filename: Main_hbonds.py
@Last modified by: Andres Vodopivec
@Last modified time: 2020-12-08

This program can be used to calculate quantify and plot the location of hydrogen bonds in any
Molecular Dynamic's system compatible with Gromacs5.x or higher, independent of the the
force-field used, as long as it is an Atomistic force-field.

The purpose of this program is to quantify and plot the hydrogen bond interaction between two
different species of molecules. The two different species of molecules are referred as:
molec0 and molec1.

PROGRAM REQUIREMENTS:
--------------------

a) This program requires: FileManager.py, Topology.py, Additionals.py, Molecule.py, H_bond.py,
    PlotManger.py, and INPUT_USER_H_BOND_CALCS.json.
    All entries must be made in the INPUT_USER_H_BOND_CALCS.json file.

b) Entries can be made to the command line. Please use python h_bonds_calc.py --help for usage
    information.

c) The GRO files to calculate the hydrogen bonds. This grofiles should be frames from a
    simulation.

The following command line can be used to extract your simulations frames:

gmx_mpi trjconv -sep -o frame.gro -skip 1 -f sim.xtc -s sim.tpr -pbc whole.
GMX will generate a series of frameX.gro file, where X is a number from 0 to how many frames
are present in your xtc file.

If a list is supplied using the -g flag in the python commnand line, it is important that your GRO
files are named using alpha-numeric (i.e. mygro0.gro, mygro1.gro, mygro2.gro, mygro3.gro, ....).
Otherwise, the code will fail.

d) The two itp files (one for each of the two molecule species).
    We refer in them in the program as molec0_itp and molec1_itp.

e) The xtc and tpr files from which the GRO files were extracted

f) All the files mentioned in this section MUST be located in the same directory where the
    Main_hbonds.py is located.



UNDERSTANDING PROGRAM NOTATION AND CONSIDERATIONS:
-------------------------------------------------

To consider a hydrogen bond, the following considerations must be taken into account:

a) The hydrogen from a molecule species has to be chemically bonded to another atom in the same
    molecule. This is represented by X-H. Where H is the hydrogen atom, X is the atom chemically
    bonded to the hydrogen atom, and "-" represents the chemical bond.

b) The third atom has to come the other molecule species, and it is represented by Y. Therefore,
    the notation X-H---Y is the Y atom interacting with the hydrogen from the other molecule species,
    and "---" represent the dipole-dipole interaction.

c) The Y atom must be within a maximum distance threshold from the Hydrogen atom only
    (designated as dist_thres).

d) The configuration X-H---Y should comply with the cosine rule. In this cosine rule we need all
    there distances (XtoH, XtoY, and HtoY), and the angle that Y forms with X-H should have a minimum
    threshold (designated as angle_thres) in order to be considered a hydrogen bond.

e) Hydrogen bonds cannot take place if Y atom is an Alkali (i.e. sodium, potassium, or calcium) or
    another hydrogen atom.



IMPORTANT:
----------

a) All files must be compatible with Gromacs in order for this program to work.

b) Please be sure to read carefully all the help info and comments to clear any doubts.

c) The hydrogen bond criteria in this program are supported in the publication:

Experimental and Computational Studies of Choline Chloride-Based Deep Eutectic Solvents
Sasha L. Perkins, Paul Painter, and Coray M. Colina
Journal of Chemical & Engineering Data 2014 59 (11), 3652-3662
DOI: 10.1021/je500520h
PROGRAM ADDITIONALS:
This program also calculates and plot a density profile plot between molec0, molec1, and the locaiton
of the hydrogen bonds. Gromacs must be installed in order to determine the density profile.

If there is no interest in doing this extra analysis, just uncomment lines from the sections:
    'Getting the density profile using gmx density' and 'Doing the density profile plot'.


USEFUL TIP:
-----------

If you do not want to carry so many files, you can create a python executable using pyinstaller.
Make sure all the python files required are in the same directory and run:
    pyinstaller --onefile Main_hbonds.py

The executable will be located inside the dist folder. This excutable is all you need and you
can run it in any directory, just make sure you have the gromacs files and Json file in the
same directory.


'''
# ==================================================================================================
# Libraries
# ==================================================================================================

# Importing modules
# -----------------
import os
import traceback
import argparse
import concurrent.futures
from FileManager import FileUtil
from Topology import Topology
from Additionals import JsonConfig as conf
from Additionals import AdditionalFunctions as add
from PlotManger import HbondsPlotter
from PlotManger import XvgPlotter

# ==================================================================================================
# Multiprocessing function
# ==================================================================================================

def run_multiprocessing(a_grofileName):

    # --------------------------------------------------------------------
    # Getting topology for molec0 and molec1
    # --------------------------------------------------------------------
    itpmolec0 = Topology.fromItpFile(workdir, molec0_ipt, Hydros0)
    itpmolec0.fromGroFile(workdir, a_grofileName)

    itpmolec1 = Topology.fromItpFile(workdir, molec1_itp, Hydros1)
    itpmolec1.fromGroFile(workdir, a_grofileName)

    # --------------------------------------------------------------------
    # Getting the molecules0 and molecules1 corresponding to the itpmolec0 and itpmolec1, respectvely
    # --------------------------------------------------------------------
    grofile = FileUtil.FileReader(workdir, a_grofileName).file

    print('\nGetting all molec0 molecules present in grofile: {}\n'.format(a_grofileName))
    molecs0_dict = add.get_molecules(itpmolec0, grofile)

    print('\nGetting all molec1 molecules present in grofile: {}\n'.format(a_grofileName))
    molecs1_dict = add.get_molecules(itpmolec1, grofile)

    # --------------------------------------------------------------------
    # Getting h-bonds interactions
    # --------------------------------------------------------------------

    # h-bonds interactions between molec0 and molec1. Where in molec0 we get X-H bond and in molec 1 we get Y
    # ------------------------------
    print('\nThe hydrogen bonds found between X-H from molec0 and Y from molec1 present in grofile {} are:\n'.format(a_grofileName))
    hbonds0_list = add.get_bonded_interactions(molecs0_dict, itpmolec0, molecs1_dict, Hydros1, Alkali1, dist_thres, angle_thres)

    if not hbonds0_list:
        print('\nNo hydrogen bonds found between X-H from molec0 and Y from molec1 present in grofile {}\n'.format(a_grofileName))
    else:
        print('\nA total of {} hydrogen bonds found between X-H from molec0 and Y from molec1 present in grofile {}\n'.format(str(len(hbonds0_list)), a_grofileName))

    # h-bonds interactions between molec0 and molec1. Where in molec1 we get X-H bond and in molec0 we get Y
    # ------------------------------
    print('\nThe hydrogen bonds found between X-H from molec1 and Y from molec0 present in grofile {} are:\n'.format(a_grofileName))
    hbonds1_list = add.get_bonded_interactions(molecs1_dict, itpmolec1, molecs0_dict, Hydros0, Alkali0, dist_thres, angle_thres)

    if not hbonds1_list:
        print('\nNo hydrogen bonds found between X-H from molec1 and Y from molec0 present in grofile {}\n'.format(a_grofileName))
    else:
        print('\nA total of {} hydrogen bonds found between X-H from molec1 and Y from molec0 present in grofile {}\n'.format(str(len(hbonds1_list)), a_grofileName))

    # adding both list
    # ------------------------------
    hbond_list = hbonds0_list + hbonds1_list # we just add them

    del(itpmolec0, itpmolec1, grofile, molecs0_dict, molecs1_dict, hbonds0_list, hbonds1_list) # delete things to free space in case we run out of memory during the multiprocessing

    return hbond_list

# ==================================================================================================
# MAIN (PROGRAM IMPLEMENTATION)
# ==================================================================================================

# -------------------------------------------------------------------------------------------------
# Input variables from INPUT_USER_H_BOND_CALCS.json
# -------------------------------------------------------------------------------------------------

# Hydrongen bonds considerations
# ------------------------------
dist_thres = conf.dist_thres.value # maximum distance required to be considered a hydrogen bond interaction
angle_thres = conf.angle_thres.value # minimum angle required to be considered a hydrogen bond interaction

# itp files
# ------------------------------
molec0_ipt = conf.molec0_itp.value # molec0 itp file name
molec1_itp = conf.molec1_itp.value # molec1 itp file name

# Hydrongen atomtypes lists
# ------------------------------
Hydros0 = conf.molec0_H_atomTypes.value.split(', ') # List of hydrogen atomtypes present in the molec0 itp file (i.e. potassium, calcium)
Hydros1 = conf.molec1_H_atomTypes.value.split(', ') # List of hydrogen atomtypes present in the molec1 itp file (i.e. potassium, calcium)

# Alkali atomtypes lists.
# ------------------------------
# This atomtypes are not be considered in case they show up as a hydrogen bond.
Alkali0 = conf.molec0_AlkaliAtomtypes.value.split(', ') # List of alkali atomtypes present in the molec0 itp file (i.e. potassium, calcium)
Alkali1 = conf.molec1_AlkaliAtomtypes.value.split(', ') # List of alkali atomtypes present in the molec1 itp file (i.e. potassium, calcium)

# Gromacs files
# ------------------------------
gmx_call = conf.gmx_call.value # to execute gmx
xtc_filename = conf.xtc_filename.value # name of the xtc file (see Additionals.py for more info)
tpr_filename = conf.tpr_filename.value # name of the tpr file (see Additionals.py for more info)
outdens_xvg = conf.outdens_xvg.value # name of the output file when doing gmx density commandline


# -------------------------------------------------------------------------------------------------
# Using argparse to get the names of the grofiles and folder to save results
# -------------------------------------------------------------------------------------------------
parser = argparse.ArgumentParser()

parser.add_argument('-g', '--gro', nargs='+',
                    default=[file for file in os.listdir(os.getcwd()) if file.find('.gro') >= 0],
                    help='grofile list input (example: system0.gro, system0.gro, ...)',
                    type=str)

parser.add_argument('-t', '--time', default=1,
                    help='Total simulation time in ns (default = number of grofiles)',
                    type=float)

args = parser.parse_args()
grofiles = sorted(args.gro, key=lambda x: int("".join([i for i in x if i.isdigit()])))

if args.time == 1: # if --time has no entry
    sim_time = len(grofiles) # sim_time will be just the grofiles supplied. Meaning it will just be simulation frame number
else:
    sim_time = args.time

# sim_time = 1
# grofiles = ['test.gro']#, 'test1.gro']
frame_time = sim_time / len(grofiles)

# -------------------------------------------------------------------------------------------------
# Running the program in multiprocessing
# -------------------------------------------------------------------------------------------------
workdir = os.getcwd()
all_hbonds = [] # all hydrogen bonds present in the grofiles will be saved
print('\nRunning multiprocessing function for all grofiles\n')

try:

    print('\n\nThe gro files to use are: \n')

    for grofile in grofiles: # printing the list of grofiles that will be analyzed
        print(grofile)

    with concurrent.futures.ProcessPoolExecutor() as executor: # Doing multiprocessing
        results = executor.map(run_multiprocessing, grofiles) # Results is a list of H_Bond objects (one for every h-bond found)

        for result in results: # saving the results. H_bond objects in this case.
            all_hbonds.append(result) # a list of results

except Exception:
    print('\n\nCaught an error:\n----------------')
    traceback.print_exc()

# -------------------------------------------------------------------------------------------------
# Plotting the number of hydrogen bonds for each one of the different grofiles
# -------------------------------------------------------------------------------------------------
plotObj = HbondsPlotter(all_hbonds) # creating the HbondsPlotter object from the all_bonds list
plotObj.plot(workdir, frame_time) # method to print all the h-bonds plots.
print('\nFinished all plots!\n')

# -------------------------------------------------------------------------------------------------
# Getting the density profile using gmx density
# -------------------------------------------------------------------------------------------------

# Getting the residue name for molec0 and molec1 (see Topology.py for more info)
# ------------------------------

molec0_residue= Topology.fromItpFile(workdir, molec0_ipt, Hydros0).residue
molec1_residue = Topology.fromItpFile(workdir, molec1_itp, Hydros1).residue

# Getting the z-axis length of the grofile (any grofile can be used)
# ------------------------------

zBoxLen = float(FileUtil.FileReader(workdir, grofiles[-1]).file[-1].split()[2]) # Any grofile can be used. The last was selected.

# Getting the density file file
# ------------------------------

add.get_gmx_index_density_files(workdir, all_hbonds, molec0_residue, molec1_residue,
                            zBoxLen, gmx_call, grofiles[-1], xtc_filename,
                            tpr_filename, outdens_xvg)

print('\nDensity profile was calculated sucessfully!\n')

# -------------------------------------------------------------------------------------------------
# Doing the density profile plot
# -------------------------------------------------------------------------------------------------

print('\nPlotting the density profile!\n')

xbgPlotter_Obj = XvgPlotter.from_XVGfile(workdir, outdens_xvg)
xbgPlotter_Obj.featureset_scaling()
xbgPlotter_Obj.plot()

print('\nDensity profile was plotted sucessfully!\n')
print('\nProgram ended sucessfully!\n')




