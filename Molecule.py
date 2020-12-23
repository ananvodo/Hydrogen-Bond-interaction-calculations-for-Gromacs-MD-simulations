#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@Author: Andres Vodopivec
@Date: 2020-12-08
@Filename: Molecule.py
@Last modified by: Andres Vodopivec
@Last modified time: 2020-12-08

Information:
------------
This file contains a Class that stores the information of each type of molecule
present in the grofile.

'''

from FileManager import FileUtil
from enum import Enum

# ==================================================================================================
# ENUM CLASS FOR GROFILE PARSING
# ==================================================================================================

class Gro(Enum):
    '''
    Enum for the start and end of the components in the grofile
    Please refer to GRO file format in the Groamcas manual for more information
    S: stand for start
    E: stand for end
    the number belong to positions in the grofile lines.
    '''
    RESNUM_S = 0
    RESNUM_E = 5

    RESTYPE_S = 5
    RESTYPE_E = 10

    ATONAME_S = 10
    ATONAME_E = 15

    ATONUM_S = 15
    ATONUM_E = 20

    XCOORD_S = 20
    XCOORD_E = 28

    YCOORD_S = 28
    YCOORD_E = 36

    ZCOORD_S = 36
    ZCOORD_E = 44


# ==================================================================================================
# CLASS FOR MOLECULE COMPONENTS
# ==================================================================================================

class Molecule():
    '''Class that stores all data from each molecule in the grofile'''

    def __init__(self, residNum, residType, atomName, atomNum, xCoord, yCoord,
                 zCoord):
        """
        This is the Constructor to instanciate the object for Molecule() class.
        For information about the attributes, please refer to Gromacs manual.
        """
        self.residNum = residNum
        self.residType = residType
        self.atomName = atomName
        self.atomNum = atomNum
        self.xCoord = xCoord
        self.yCoord = yCoord
        self.zCoord = zCoord


# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------


    @classmethod
    def fromList(cls, a_indexesList, a_grofileList):
        '''
        Molecule.fromList is the class constructor. It only takes the a_grofileList.
        a_grofileList is a piece from the grofile that corresponds to the
        molecule.
        This method is to about the molecules in the grofile.
        a_indexesList are the indexes to be selected from the supplied grofile (a_grofileList)
        '''
        # -------------------------------------
        # Verifying the input format is correct
        # -------------------------------------
        FileUtil.CheckFuncInput(a_grofileList, list, 'Topology.fromFileUtil')

        # -------------------------------------
        # splitting the line into the desired components
        # -------------------------------------
        residNum = []
        residType = []
        atomName = []
        atomNum = []
        xCoord = []
        yCoord = []
        zCoord = []

        for line in a_indexesList:
            residNum.append(int(a_grofileList[line][Gro.RESNUM_S.value : Gro.RESNUM_E.value].strip()))
            residType.append(a_grofileList[line][Gro.RESTYPE_S.value : Gro.RESTYPE_E.value].strip())
            atomName.append(a_grofileList[line][Gro.ATONAME_S.value : Gro.ATONAME_E.value].strip())
            atomNum.append(int(a_grofileList[line][Gro.ATONUM_S.value : Gro.ATONUM_E.value].strip()))
            xCoord.append(float(a_grofileList[line][Gro.XCOORD_S.value : Gro.XCOORD_E.value].strip()))
            yCoord.append(float(a_grofileList[line][Gro.YCOORD_S.value : Gro.YCOORD_E.value].strip()))
            zCoord.append(float(a_grofileList[line][Gro.ZCOORD_S.value : Gro.ZCOORD_E.value].strip()))

        return cls(residNum, residType, atomName, atomNum, xCoord, yCoord,
                   zCoord)


# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

    def fromToplogy(self, a_atomtypeList):
        '''
        Molecule.fromToplogy method is to create an attibute with all the atomtypes.
        The atomtypes are not present in the GRO file. They are present in the
        itp file for the specific molecule.
        For more information, please refer to Gromacs manual.
        '''
        self.atomtype = a_atomtypeList

        return None







































