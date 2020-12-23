#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@Author: Andres Vodopivec
@Date: 2020-12-08
@Filename: Topology.py
@Last modified by: Andres Vodopivec
@Last modified time: 2020-12-08

Information:
------------
This file contains a Class that stores the topology of each type of molecule
present in the GRO file.

'''

from FileManager import FileUtil



class Topology():
    '''
    Class that stores the topology of each type of molecule present in the grofile.
    For information about the attributes, please refer to Gromacs manual.
    Specifically ipt file (input topology file).
    '''

    def __init__(self, moleculetype, residue, atoms_per_molec, atomtype_atomNum_dict, bonds, atomTypeList):
        self.moleculetype = moleculetype
        self.residue = residue
        self.atoms_per_molec = atoms_per_molec
        self.atomtype_atomNum_dict = atomtype_atomNum_dict
        self.bonds = bonds
        self.atomTypeList = atomTypeList

# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

    @classmethod
    def fromItpFile(cls, a_destinationFilePath, a_itpfileName, a_atomtypesH):
        '''
        This is the class constructor. It only takes the a_itpfile and a_atomtypesH.
        a) a_itpfile is the name of the gromacs itp file for a given molecule
        b) a_atomtypesH is a list of H atoms in the itp file
        c) a_destinationFilePath is where the file a_itpfileName is located
        '''
        # -------------------------------------
        # Verifying the input format is correct
        # -------------------------------------
        FileUtil.CheckFuncInput(a_destinationFilePath, str, 'Topology.fromItpFile')
        FileUtil.CheckFuncInput(a_itpfileName, str, 'Topology.fromItpFile')
        FileUtil.CheckFuncInput(a_atomtypesH, list, 'Topology.fromItpFile')

        # -------------------------------------
        # Reading the itp file
        # -------------------------------------
        itpFileOb = FileUtil.FileReader(a_destinationFilePath, a_itpfileName)
        itpfile = itpFileOb.file

        # -------------------------------------
        # Getting atoms from the itpfile
        # -------------------------------------
        # [ moleculetype ], [ atoms ], [ bonds ], and [ pairs ] are part of itp file format in Gromacs
        for index, line in enumerate(itpfile):
            if line.find('[ moleculetype ]') >= 0:
                molname_index = index + 2
            if line.find('[ atoms ]') >= 0:
                atoms_index_top = index + 2
            if line.find('[ bonds ]') >= 0 or line.find('[ pairs ]') >= 0:
                atoms_index_bottom = index - 1
                break

        moleculetype = itpfile[molname_index].split()[0]
        atoms_section = itpfile[atoms_index_top : atoms_index_bottom]
        atoms_per_molec = len(atoms_section)
        residue = atoms_section[0].split()[3]

        atomtype_atomNum_dict = {}
        atomTypeList = []

        for line in atoms_section:
            parts = line.split()
            key_atomtype = parts[1].strip()
            atomTypeList.append(key_atomtype)
            atom_nr = int(parts[0])

            if key_atomtype not in atomtype_atomNum_dict.keys():
                atomtype_atomNum_dict[key_atomtype] = []
                atomtype_atomNum_dict[key_atomtype].append(atom_nr)
            else:
                atomtype_atomNum_dict[key_atomtype].append(atom_nr)

        # -------------------------------------
        # Getting index of atoms bonded to Hydrogens from the itpfile
        # -------------------------------------
        # specifically looking in the bonds section
        # H is the hydrogen atom
        # X is the atom chemically bonded to the H atom
        for index, line in enumerate(itpfile):
            if line.find('[ bonds ]') >= 0:
                bonds_index_top = index + 2
            if line.find('[ angles ]') >= 0:
                bonds_index_bottom = index - 1
                break


        bonds_section = itpfile[bonds_index_top : bonds_index_bottom]
        atomsH_nrL = []

        for key_atomH in a_atomtypesH:
            atomsH_nrL.extend(atomtype_atomNum_dict[key_atomH])

        atomsH_nrL = sorted(atomsH_nrL)
        bonds = []

        for atomH_nr in atomsH_nrL:
            atomH_line = atoms_section[atomH_nr - 1]
            partsH = atomH_line.split()
            atyomtypeH = partsH[1]
            atomH = partsH[4]

            for line in bonds_section:
                # bonds sections is a list where each line has two number. Those are atomnumbers that corresponds to two atoms that are bonded.
                # the H atom can be any of those numbers, thus we need to check both.
                parts = line.split()
                atomnum1 = int(parts[0])
                atomnum2 = int(parts[1])
                bond_info = {}

                if atomnum1 == atomH_nr: # if the hydrogen is atomnum1
                    bond_info['atomNumH'] = atomnum1
                    bond_info['atomTypeH'] = atyomtypeH
                    bond_info['atomH'] = atomH

                    # Getting AtomX
                    # -------------
                    atomX_line = atoms_section[atomnum2 - 1] # atomnums start from 1 and python is starting from 0. Thus we take one unit.
                    partsX = atomX_line.split()
                    atyomtypeX = partsX[1]
                    atomX = partsX[4]
                    bond_info['atomNumX'] = atomnum2
                    bond_info['atomTypeX'] = atyomtypeX
                    bond_info['atomX'] = atomX
                    bonds.append(bond_info)

                if atomnum2 == atomH_nr: # if the hydrogen is atomnum2
                    bond_info['atomNumH'] = atomnum2
                    bond_info['atomTypeH'] = atyomtypeH
                    bond_info['atomH'] = atomH

                    # Getting AtomX
                    # -------------
                    atomX_line = atoms_section[atomnum1 - 1] # atomnums start from 1 and python is starting from 0. Thus we take one unit.
                    partsX = atomX_line.split()
                    atyomtypeX = partsX[1]
                    atomX = partsX[4]
                    bond_info['atomNumX'] = atomnum1
                    bond_info['atomTypeX'] = atyomtypeX
                    bond_info['atomX'] = atomX
                    bonds.append(bond_info)

                if atomnum1 > atomH_nr and atomnum2 > atomH_nr: # if both atomnum are higher that atomH_nr, then there is no need to continue the search for that H atom.
                    break

        return cls(moleculetype, residue, atoms_per_molec, atomtype_atomNum_dict, bonds, atomTypeList) # these 5 variable are basic data from a ipt file


# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

    def fromGroFile(self, a_destinationFilePath, a_grofileName):
        '''
        Tolopogy.fromGroFile takes the a_destinationFilePath, a_topolfileName,
        and a_atomtypesH.
        a) a_grofileName is the name of the grofile for the given system
        b) a_destinationFilePath is where the file a_topolfileName is located
        The return of this function are the lines_in_grofile and number_molecs
        present inf the grofile.
        '''
        # -------------------------------------
        # Verifying the input format is correct
        # -------------------------------------
        FileUtil.CheckFuncInput(a_destinationFilePath, str, 'Topology.fromGroFile')
        FileUtil.CheckFuncInput(a_grofileName, str, 'Topology.fromGroFile')

        # -------------------------------------
        # Reading the GRO file
        # -------------------------------------
        groFileOb = FileUtil.FileReader(a_destinationFilePath, a_grofileName)
        grofile = groFileOb.file[2:-1]
        lines_in_grofile = [] # taking the index of the lines from the GRofile that correspond to the molecule described in the itp file.

        for index, line in enumerate(grofile):
            residue = line[5:10].strip()

            if residue == self.residue:
                lines_in_grofile.append(index)

        self.lines_in_grofile = lines_in_grofile
        self.number_molecs = int(len(self.lines_in_grofile) / self.atoms_per_molec) # number of molecules in the GROfile that correspond to the topology given.

        return None

