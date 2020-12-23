#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@Author: Andres Vodopivec
@Date: 2020-12-08
@Filename: H_bonds.py
@Last modified by: Andres Vodopivec
@Last modified time: 2020-12-08

Information:
------------
This file contains a Class that stores the Hydrogen bonds present in the
grofile

'''

import math

class H_bond():
    '''
    Class to store all the h-bond information.
    There are three atoms involved in a h-bond interaction, here represented as H, X and Y.
    Where H is the hydrogen atom, X is the atom chemically bonded to the hydrogen, and
    Y is the atom interacting the hydrogen through hydrogen bond interactions (dipole-diople attraction).

    Representation:
    X-H---Y
    where:
        - represents chemical bond
        --- represent dipole-dipole interaction (hydrogen bond interaction)

    For more information about all the attributes, please refer to the Gromacs manual for itp files and grofiles.
    '''
    def __init__(self, H_atomType, H_atomNum, H_atomName,
                 X_atomType, X_atomNum, X_atomName,
                 xcoordH, ycoordH, zcoordH,
                 xcoordX, ycoordX, zcoordX,
                 distance_HtoX):

        self.H_atomType = H_atomType
        self.H_atomNum = H_atomNum
        self.H_atomName = H_atomName

        self.X_atomType = X_atomType
        self.X_atomNum = X_atomNum
        self.X_atomName = X_atomName

        self.Y_atomType = None
        self.Y_atomNum = None
        self.Y_atomName = None

        self.xcoordH = xcoordH
        self.ycoordH = ycoordH
        self.zcoordH = zcoordH

        self.xcoordX = xcoordX
        self.ycoordX = ycoordX
        self.zcoordX = zcoordX

        self.xcoordY = None
        self.ycoordY = None
        self.zcoordY = None

        self.distance_HtoX = distance_HtoX
        self.distance_HtoY = None
        self.distance_XtoY = None
        self.is_H_bond = False

# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------


    @classmethod
    def fromBondAndMolecobj(cls, a_bond_info_dict, a_molecObj):
        indexH = a_bond_info_dict['atomNumH'] - 1 # atomnums start from 1 and python is starting from 0. Thus we take one unit.
        indexX = a_bond_info_dict['atomNumX'] - 1 # atomnums start from 1 and python is starting from 0. Thus we take one unit.

        H_atomType = a_bond_info_dict['atomTypeH']
        H_atomNum = a_molecObj.atomNum[indexH]
        H_atomName = a_molecObj.atomName[indexH]
        xcoordH = a_molecObj.xCoord[indexH]
        ycoordH = a_molecObj.yCoord[indexH]
        zcoordH = a_molecObj.zCoord[indexH]

        X_atomType = a_bond_info_dict['atomTypeX']
        X_atomNum = a_molecObj.atomNum[indexX]
        X_atomName = a_molecObj.atomName[indexX]
        xcoordX = a_molecObj.xCoord[indexX]
        ycoordX = a_molecObj.yCoord[indexX]
        zcoordX = a_molecObj.zCoord[indexX]

        x_delta = xcoordH - xcoordX
        y_delta = ycoordH - ycoordX
        z_delta = zcoordH - zcoordX

        distance_HtoX = math.sqrt(x_delta ** 2 + y_delta ** 2 + z_delta ** 2)

        return cls(H_atomType, H_atomNum, H_atomName,
                   X_atomType,X_atomNum, X_atomName,
                   xcoordH, ycoordH, zcoordH,
                   xcoordX, ycoordX, zcoordX,
                   distance_HtoX)

        return None






