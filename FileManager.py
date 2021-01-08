#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@Author: Andres Vodopivec
@Date: 2020-12-08
@Filename: FileManager.py
@Last modified by: Andres Vodopivec
@Last modified time: 2020-12-08

Information:
------------
This file has the class FileUtilBase that is used for read and write files as
list of str only

'''

from abc import ABC, abstractmethod
import os


# ==================================================================================================
# PlotFileUtilBaseManager abstract class
# ==================================================================================================

class FileUtilBase(ABC):

    def __init__(self, file):
        self.file = file

    # file.getter
    @property
    def file(self):
        return self.__file


    @file.setter
    def file(self, file):
        if isinstance(file, list):

            if all([isinstance(element, str) for element in file]):
                self.__file = file
            else:
                raise ValueError('All elements in file supplied for FileUtilBase are not str')

        else:
            raise ValueError('\nfile supplied for FileUtilBase is not a list\n')

# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

    @classmethod
    @abstractmethod
    def FileReader(cls, a_destinationFilePath, a_fileName):
        """
        FileUtilBase.FileReader is an abstractmethod to read a file from
        /destinationFilePath/fileName to be implemented in the child classes.
        """
        pass

# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

    def FileWriter(self, a_destinationFilePath, a_fileName):
        """
        FileUtilBase.FileWriter is to print an object to /destinationFilePath/fileName.
        """

        print('\nWriting the new file with name ' + a_fileName +
              ' in directory: ' + a_destinationFilePath + '\n')

        # Confirming filepath exists
        # ------------------------
        workdir = FileUtilBase.MakePathDir('FileUtilBase.FileWriter', a_destinationFilePath)
        filepath = os.path.join(workdir, a_fileName)

        # Prinitng file attribute
        # ------------------------
        with open(filepath, 'w') as fout:
            fout.writelines(self.file)

        return None

# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

    @staticmethod
    def CheckFuncInput(a_inputVar, a_inputType, a_classDotFuncName):
        '''
        FileUtilBase.CheckFuncInput is a static method to check that input to any
        function in other classes has the correct type (i.e. string, int, ...)
        '''
        if isinstance(a_inputVar, a_inputType):
            pass
        else:
            raise ValueError('\nArgument {} in {} is not a/an {}\n'.format(
                    a_inputVar,
                    a_classDotFuncName,
                    str(a_inputType).split()[1]))

        return None

# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

    @staticmethod
    def MakePathDir(a_classDotFuncName, a_destinationFilePath,
                    a_subFolder=None):
        '''
        FileUtilBase.MakePathDir is a static method to create the path for any
        file that is required as input for any function.
        This funtions looks that the directory provided is string format and
        that it exists.
        '''
        FileUtilBase.CheckFuncInput(a_destinationFilePath, str,
                                a_classDotFuncName)

        if a_subFolder is None:
            workPath = a_destinationFilePath
            # print("\nUsing {} in path {}\n".format(
            #         a_classDotFuncName, workPath))
        else:
            FileUtilBase.CheckFuncInput(a_subFolder, str, a_classDotFuncName)
            workPath = os.path.join(a_destinationFilePath, a_subFolder)
            # print("\nUsing {} in path {}\n".format(
            #         a_subFolder, workPath))

        if not os.path.exists(workPath):
            raise FileNotFoundError('{} did not found path {}'.format(
                    a_classDotFuncName, workPath))

        return workPath


# ==================================================================================================
# FileUtil class
# ==================================================================================================

class FileUtil(FileUtilBase):

    def __init__(self, file, **kwargs): # **kwargs serves as an initialization list. Here is not required here because file is the only attribute and there are no more initialized attributes
        super().__init__(file, **kwargs)

# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

    @classmethod
    def FileReader(cls, a_destinationFilePath, a_fileName):
        """
        FileUtil.FileReader is a Constructor that takes a path
        (/destinationFilePath/fileName) to instanciate the object.
        """

        print('\nReading new file with name ' + a_fileName +
              ' in directory: ' + a_destinationFilePath + '\n')


        filepath = os.path.join(a_destinationFilePath, a_fileName)

        # Confirming filepath exists
        # ------------------------
        filepath = FileUtilBase.MakePathDir('FileUtilBase.FileReader', filepath)

        # Reading the file
        # ------------------------
        file = []

        with open(filepath, 'r') as fin:
            file = fin.readlines()

        return cls(file)

# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

    def FileWriter(self, a_destinationFilePath, a_fileName): # This is not required, but it is just for example how to use a parent class implementation
        super().FileWriter(a_destinationFilePath, a_fileName)



# x = FileUtil.FileReader(os.getcwd(), 'test.gro')

# x.FileWriter(os.getcwd(), 'xx.gro')

