#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@Author: Andres Vodopivec
@Date: 2020-12-08
@Filename: Additionals.py
@Last modified by: Andres Vodopivec
@Last modified time: 2020-12-08

Usage:
------
Thia file constist only of static functions for different activities required.
Also the class to read the JSON file for the input variables.

'''

from enum import Enum
import json
from Molecule import Molecule
from H_bond import H_bond
import math
import subprocess
import os
from FileManager import FileUtil

# ==================================================================================================
# AdditionalFunctions abstract class
# ==================================================================================================

class AdditionalFunctions:
    '''AdditionalFunctions is a class composed of only staticmethods'''

# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

    @staticmethod
    def get_molecules(a_itpmolecObj, a_grofile):
        '''
        AdditionalFunctions.get_molecules is a method that returns a dictionary
        that has all the molecules of the same species present in the a_grofile
        provided.
        Example of dictionary:
            molecs_in_system = {[Mol0]: {Molecule Object}, [Mol1]: {Molecule Object},
                                [Mol2]: {Molecule Object}, .....}
        where Mol0 to MolN belong to the same molecule species (same itp file description).
        '''

        molecs_in_system = {}
        number_molecs = a_itpmolecObj.number_molecs
        atoms_per_molec = a_itpmolecObj.atoms_per_molec
        count = 0

        for i in range(0, number_molecs):
            molec_start = i * atoms_per_molec
            molec_end = molec_start + atoms_per_molec
            indexes_slice = a_itpmolecObj.lines_in_grofile[molec_start : molec_end]

            molec_id = a_itpmolecObj.residue + '_' + str(count)
            molecOb = Molecule.fromList(indexes_slice, a_grofile[2:-1])
            molecOb.fromToplogy(a_itpmolecObj.atomTypeList)

            molecs_in_system[molec_id] = molecOb
            count += 1

        delattr(a_itpmolecObj, 'lines_in_grofile') # this attribute is not needed anymore
        delattr(a_itpmolecObj, 'atomTypeList') # this attribute is not needed anymore

        return molecs_in_system

# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

    @staticmethod
    def get_bonded_interactions(a_molecs0_dict, a_itpmolec0, a_molecs1_dict, a_atomtypesH1,
                                a_atomtypesAlkali1, a_dist_thres, a_angle_thres):
        '''
        AdditionalFunctions.get_bonded_interactions is a method that get the hydrogen bonding interactions.
        It returns a list of H_bond objects (see H_bond.py for more information).
        This funtion return a list of H_bond objects (see H_Bond.py for more information).
        '''

        hBonds0_true = []
        hBonds0_list = []
        atomtypes_not_to_consider = a_atomtypesH1 + a_atomtypesAlkali1

        for molec0Name, molec0Info in a_molecs0_dict.items():
            for bond in a_itpmolec0.bonds:
                h_bondObj0 = H_bond.fromBondAndMolecobj(bond, molec0Info)
                hBonds0_list.append(h_bondObj0)

        for molec1Name, molec1Info in a_molecs1_dict.items():
            xcoordY_list = molec1Info.xCoord
            ycoordY_list = molec1Info.yCoord
            zcoordY_list = molec1Info.zCoord

            for i in range(len(xcoordY_list)):

                for bond in hBonds0_list:
                    xcoordY = xcoordY_list[i]
                    ycoordY = ycoordY_list[i]
                    zcoordY = zcoordY_list[i]

                    x_diff = xcoordY - bond.xcoordH
                    y_diff = ycoordY - bond.ycoordH
                    z_diff = zcoordY - bond.zcoordH

                    distance_HtoY = math.sqrt(x_diff ** 2 + y_diff ** 2 + z_diff ** 2)

                    if distance_HtoY <= a_dist_thres:
                        x_diff = xcoordY - bond.xcoordX
                        y_diff = ycoordY - bond.ycoordX
                        z_diff = zcoordY - bond.zcoordX

                        distance_XtoY = math.sqrt(x_diff ** 2 + y_diff ** 2 + z_diff ** 2)
                        distances = sorted([distance_XtoY, bond.distance_HtoX, distance_HtoY])
                        angle = math.degrees(math.acos((distances[0]**2 + distances[1]**2 - distances[2]**2) / (2 * distances[0] * distances[1])))

                        if angle > a_angle_thres:
                            Y_atomtype = molec1Info.atomtype[i]

                            if Y_atomtype not in atomtypes_not_to_consider:
                                bond.xcoordY = xcoordY_list[i]
                                bond.ycoordY = ycoordY_list[i]
                                bond.zcoordY = zcoordY_list[i]

                                bond.Y_atomType = molec1Info.atomtype[i]
                                bond.Y_atomNum = molec1Info.atomNum[i]
                                bond.Y_atomName = molec1Info.atomName[i]
                                print(bond.X_atomType, '-', bond.H_atomType, '---', bond.Y_atomType)

                                bond.distance_HtoY = distance_HtoY
                                bond.distance_XtoY = distance_XtoY
                                bond.is_H_bond = True
                                hBonds0_true.append(bond)

        return hBonds0_true

# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

    @staticmethod
    def chunks(lst, n):
        """AdditionalFunctions.chunks yields successive n-sized chunks from lst."""
        for i in range(0, len(lst), n):
            yield lst[i:i + n]


    @staticmethod
    def for_gmx_index(a_gmxCall, a_grofile):
        '''
        AdditionalFunctions.for_gmx_index call gromacs to get the index.ndx file
        '''
        print('\nCalling gmx make_ndx to create the index.ndx file\n')
        commandline = '(echo "q\n") | {} make_ndx -f {}'.format(a_gmxCall, a_grofile)

        print('\nUsing the coomandline: {}\n'.format(commandline))
        gmx_outputFile = open('index_output.txt', 'w')
        proc = subprocess.Popen(commandline, stderr=subprocess.STDOUT, stdout=gmx_outputFile, shell=True)
        streamdata = proc.communicate()
        outproc = proc.returncode
        gmx_outputFile.close()

        if outproc != 0:
            raise Exception('gmx make_ndx in commandline {} from for_gmx_index() failed'.format(commandline))

        print('\ngmx make_ndx finished succesfully!\n')

        return 'index_output.txt' # returning the name of the stdout and stderr from gmx make_ndx. It will be needed further

# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

    @staticmethod
    def for_gmx_density(a_gmxCall, a_grofile, a_xtcName, a_tprName, a_outdens_xvg,
                        a_zBoxLen, a_molec0_index, a_molec1_index, a_h_bond_index):
        '''
        AdditionalFunctions.for_gmx_density call gromacs to get the density profile using gmx density.
        The output file will have the name used for the argument a_outdens_xvg.
        '''

        print('\nCalling gmx density to create the {}  file\n'.format(a_outdens_xvg))
        commandline = '(echo "{}\n"; echo "{}\n"; echo "{}\n"; echo "q\n") | {} density -s {} -f {} -o {} -sl {} -d z -ng 3 -dens mass -n index.ndx'.format(
            a_molec0_index, a_molec1_index, a_h_bond_index, a_gmxCall, a_tprName, a_xtcName, a_outdens_xvg, str(a_zBoxLen))

        print('\nUsing the coomandline: {}\n'.format(commandline))
        gmx_outputFile = open('density_output.txt', 'w')
        proc = subprocess.Popen(commandline, stderr=subprocess.STDOUT, stdout=gmx_outputFile, shell=True)
        streamdata = proc.communicate()
        outproc = proc.returncode
        gmx_outputFile.close()

        if outproc != 0:
            raise Exception('gmx density in commandline {} from for_gmx_index() failed'.format(commandline))

        print('\ngmx density finished succesfully!\n')

        return None

# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

    @staticmethod
    def get_gmx_index_density_files(a_destinationFilePath, a_all_hbonds, a_molec0_residue, a_molec1_residue,
                                    a_zBoxLen, a_gmx_call, a_grofileName, a_xtc_filename,
                                    a_tpr_filename, a_outdens_xvg):
        '''
        AdditionalFunctions.get_gmx_index_density_files is a method that call AdditionalFunctions.for_gmx_index
        and AdditionalFunctions.for_gmx_density.
        '''

        # Getting the index.ndx file using gmx make_ndx
        # ------------------------------

        columns_in_ndxFile = 15 # index files from gromacs has by default 15 columns
        ndx_gmx_format = '%6d' # every digit in the column must comply wih this format (see gromacs manual)

        # Getting ready the atomnums in the all_hbonds to be added to the gmx index.ndx file
        # ------------------------------

        atom_nrs = [] # to save all the atomnums that are present in the H_bond objects
        for hbonds in a_all_hbonds: # hbonds is a list of H_bond objects
            for hbondObj in hbonds:
                atom_nrs.extend([hbondObj.H_atomNum, hbondObj.X_atomNum, hbondObj.Y_atomNum])

        atom_nrs = list(AdditionalFunctions.chunks(atom_nrs, columns_in_ndxFile)) # atom_nrs list must be divided into chunks of 15 numbers to comply with the index format for gmx
        # Example: atom_nrs = [x,x,x,x,x,x,x,x,x,x, ...........]
        # after using the chunk function, atom_nrs = [[15 x numbers here], [15 x numbers here], [15 x numbers here], ....]
        # where all x are different numbers

        # Now that the atomnums are in chucks of 15 numbers, the gmx index format must be applied to all numbers
        # ------------------------------

        hbond_atomnums = [] # list to save the atomnums from h_binds that will be added to the index.ndx gromacs file
        for element in atom_nrs:
            formatted = ''.join(ndx_gmx_format % n for n in element) + '\n'
            hbond_atomnums.append(formatted)

        hbond_atomnums = ['[ h_bond_atoms ]\n'] + hbond_atomnums # additing the tittle for this new group in the index file (see gromacs manual for more info).

        # Getting the final index.ndx file
        # ------------------------------

        index_out_filename = AdditionalFunctions.for_gmx_index(a_gmx_call, a_grofileName) # creating the index.ndx file using gmx make_ndx and the name of the stderr and stdout

        ndx_fileObj = FileUtil.FileReader(a_destinationFilePath, 'index.ndx') # reading the index.ndx file created in the previous code
        ndx_fileObj.file = ndx_fileObj.file + hbond_atomnums # adding the h_bonds atomnum to the index file (all this is required to calculate the density)
        ndx_fileObj.FileWriter(a_destinationFilePath, 'index.ndx') # priting the new index.ndx file
        del(ndx_fileObj)

        print('\nThe atom_nrs that are involved in the hydrogen bonds were added at the end of the index.ndx file \n')

        # The stderr and stdout from gmx make_ndx is needed to get the indexes for gmx density input
        # ------------------------------

        # There is a bug very well known in gromacs. regards index files. Sometimes the residues are repeated.
        # When that happens, the residue name cannot be used for gmx density. To avoid problem, the index referred to the
        # residues wil be used instead.
        out_indexfile = FileUtil.FileReader(a_destinationFilePath, index_out_filename).file # reading the stderr and stdout file of gmx make_ndx

        molec0_index = [] # if the residue0 is repeated we wil need to append it to get the first index
        molec1_index = [] # if the residue1 is repeated we wil need to append it to get the first index

        for i, line in enumerate(out_indexfile):
            if line.find(a_molec0_residue) >= 0:
                molec0_index.append(line)
            if line.find(a_molec1_residue) >= 0:
                molec1_index.append(line)

        molec0_index = molec0_index[0].split()[0] # taking the first index for residue0
        molec1_index = molec1_index[0].split()[0] # taking the first index for residue1
        h_bond_index = 'h_bond_atoms' # this was the same name used in the line 274 and it is at the end of the new index.ndx file

        # -------------------------------------------------------------------------------------------------
        # Getting the density file using gmx density
        # -------------------------------------------------------------------------------------------------

        zBoxLen = int(a_zBoxLen*10) # a_zBoxLen must be multiplied by 10 to be more slices and it must an integer (see gromacs manual)
        AdditionalFunctions.for_gmx_density(a_gmx_call, a_grofileName, a_xtc_filename,
                                            a_tpr_filename, a_outdens_xvg, zBoxLen,
                                            molec0_index, molec1_index, h_bond_index)


        return None




# ==================================================================================================
# JsonConfig class
# ==================================================================================================

class JsonConfig(Enum):
    '''JsonConfig is a class to read the Json file that has all the imput variables for the program'''

    # --------------------------------------------------
    # Reading the Json file
    # --------------------------------------------------
    with open('INPUT_USER_H_BOND_CALCS.json', 'r') as fjson:
        config = json.load(fjson)

    print('\nReading input variable from INPUT_USER_H_BOND_CALCS.json file\n')
    # --------------------------------------------------
    # For molecule0
    # --------------------------------------------------
    molec0_itp = config['MOLECULE0']['GMX_ITP_FILE'] # name of the itp file of molec0
    molec0_H_atomTypes = config['MOLECULE0']['HYDROGEN_ATOMTYPES'] # list of hydrogen atomtypes in molec0 (should be in the topology description)
    molec0_AlkaliAtomtypes = config['MOLECULE0']['ALKALI_METALS_ATOMTYPES'] # list of hydrogen atomtypes in molec0 (should be in the topology description)

    # --------------------------------------------------
    # For molecule1
    # --------------------------------------------------
    molec1_itp = config['MOLECULE1']['GMX_ITP_FILE'] # name of the itp file of molec1
    molec1_H_atomTypes = config['MOLECULE1']['HYDROGEN_ATOMTYPES'] # list of hydrogen atomtypes in molec1 (should be in the topology description)
    molec1_AlkaliAtomtypes = config['MOLECULE1']['ALKALI_METALS_ATOMTYPES'] # list of hydrogen atomtypes in molec1 (should be in the topology description)

    # --------------------------------------------------
    # Criteria to consider an h-bond interaction
    # --------------------------------------------------
    dist_thres = config['H_BOND_CONSTANTS']['DISTANCE_THRESHOLD'] # maximun distance required to be considered a hydrogen bond interaction
    angle_thres = config['H_BOND_CONSTANTS']['ANGLE_THRESHOLD'] # minimum angle required to be considered a hydrogen bond interaction

    # --------------------------------------------------
    # Name for Gromacs files
    # --------------------------------------------------
    gmx_call = config['GMX_FILES']['GMX_CALL'] # sourcing gmx executable
    xtc_filename = config['GMX_FILES']['XTC_FILENAME'] # name of the xtc file (see Additionals.py for more info)
    tpr_filename = config['GMX_FILES']['TPR_FILENAME'] # name of the tpr file (see Additionals.py for more info)
    outdens_xvg = config['GMX_FILES']['OUTPUT_GMX_DENSITY_FILENAME'] # name of the output file when doing gmx density commandline (see Additionals.py for more info)















