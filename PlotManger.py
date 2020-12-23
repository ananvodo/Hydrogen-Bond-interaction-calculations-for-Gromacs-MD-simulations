#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@Author: Andres Vodopivec
@Date: 2020-12-08
@Filename: PlotManager.py
@Last modified by: Andres Vodopivec
@Last modified time: 2020-12-08

Information:
------------
This file has the class PlotManager that is used for plotting in matplotlib.

'''

from abc import ABC, abstractmethod
import os
os.environ['QT_QPA_PLATFORM']='offscreen'
import matplotlib
if os.environ.get('DISPLAY', '') == '':
    print('no display found. Using non-interactive Agg backend')
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import itertools
import numpy as np
from FileManager import FileUtil
import math

# ==================================================================================================
# PlotManager abstract class
# ==================================================================================================


class PlotManager(ABC):
    '''PlotManager is an abstract class that is used for plotting in matplotlib'''

# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

    @abstractmethod
    def plot(self, a_destinationFilePath):
        '''PlotManager.plot is an abstractmethod for plotting'''
        pass


# ==================================================================================================
# HbondsPlotter
# ==================================================================================================
# This class is only to plot h_bond information.

class HbondsPlotter(PlotManager):
    '''This class is only to plotting h_bond information.'''

    def __init__(self, allBonds_list,
                 plotname='hbonds', xlabel='Simulation time (ns)',
                 ylabel='Number of Hydrogen Bonds',
                 MaxNLocatorX_int=True, MaxNLocatorY_int=True,
                 legend_loc='best', legend_position = 'outside'):
        self.allBonds_list = allBonds_list # list of lists. Every list inside the list are composed by H_bond objects. See H_bond.py for more info.
        self.plotname = plotname
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.MaxNLocatorX_int = MaxNLocatorX_int # to use integer numbers in x-axis ticks
        self.MaxNLocatorY_int = MaxNLocatorY_int # to use integer numbers in y-axis ticks
        self.legend_loc = legend_loc # see matplotlib documentation for all possible locations for the legend.
        self.legend_position = legend_position # legend to be inside or outside the plotting area


# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------


    def plot(self, a_destinationFilePath, a_frameTime):
        '''
        HbondsPlotter.plot is an abstractmethod for plotting hydrogen bonding information
        a_destinationFilePath is the destination where the plot are to be saved.
        a_frameTime is the time per frame (grofile) (i.e. every frame is taken at every 3ns of simulation)
        '''
        print('\nGetting all information ready to plot the h-bonds per frame\n')

        # -----------------------------------------------------
        # Function with constant parameters for matplolib
        # -----------------------------------------------------
        def plot_parameters(a_plotname):
            '''
            HbondsPlotter.plot.plot_parameters that contains all the constant information
            regards matplotlib
            '''
            ax.yaxis.set_major_locator(MaxNLocator(integer=self.MaxNLocatorY_int)) # to use integer numbers in y-axis
            ax.xaxis.set_major_locator(MaxNLocator(integer=self.MaxNLocatorX_int)) # to use integer numbers in x-axis
            plt.xlabel(self.xlabel) # tittle for x-axis
            plt.ylabel(self.ylabel) # tittle for y-axis

            if self.legend_position == 'outside': # legend is to be located outside the plotting area
                ax.legend(loc=self.legend_loc, fontsize=7, bbox_to_anchor=(1,1),fancybox=True, shadow=True)
            else: # legend is to be located inside the plotting area
                ax.legend(loc=self.legend_loc, fontsize=7)

            savepath = os.path.join(a_destinationFilePath, a_plotname)
            plt.savefig(savepath, dpi=1800, bbox_inches="tight")
            # plt.show()

            return None

        # -----------------------------------------------------
        # Getting the needed information to do the plots
        # -----------------------------------------------------

        num_grofiles = len(self.allBonds_list) # the total number of grofiles that are being analyzed
        classifier_group = {} # disctionary to store all the different types of bonds in all grofiles
        # Example of classifier_group:
        # classifier_group {'C-H---O': {
        #                               'C-H---O': [x....],
        #                               'CC3162 - HCA1 --- OY1': [y.....],
        #                               'CC3161 - HCA2 --- OY2': [y.....]
        #                               }
        #                   'O-H---O': {
        #                               'O-H---O': [x....],
        #                               'OC3162 - HCA1 --- OY1': [y.....],
        #                               'OC3161 - HCA2 --- OY2': [y.....]
        #                               }
        #                   }
        # Where: x.... is the total h-bonds element in  allBonds_list for the specific group. If there are 7 elements, then x... shoud have seven numbers
        # Where: y.... is the total h-bonds element in  allBonds_list for the specific sub_group. If there are 7 elements, then x... shoud have seven numbers

        classifier_group['Total'] = [0 for i in range(num_grofiles)] # initializing to zero all the total number of bonds per grofile

        for index, hbonds in enumerate(self.allBonds_list):
            count = 0 # to count total number of bonds in the current grofile bonding info (or element in allBonds_list)

            for hbond in hbonds: # hbonds is an object of the H_bond class (see H_bond.py file for more information)
                count += 1
                key = hbond.X_atomType + '-' + hbond.H_atomType + '---' + hbond.Y_atomType # the key will look like this: CC3162 - HCA1 --- OY1
                key_group = hbond.X_atomType[0] + '-' + hbond.H_atomType[0] + '---' + hbond.Y_atomType[0] # this key will look like tis: C-H---O

                if key_group not in classifier_group.keys(): # if the group_key is not present
                    classifier_group[key_group] = {} # create the empty disctionary
                    classifier_group[key_group][key_group] = [0 for i in range(num_grofiles)] # initialize all values to zero.


                if key not in classifier_group[key_group].keys():# if the key is not present
                    classifier_group[key_group][key] = [0 for i in range(num_grofiles)] # initialize all values to zero.

                classifier_group[key_group][key][index] += 1 # if the key_group and key are present, then we update count
                classifier_group[key_group][key_group][index] += 1 # if the key_group and key are present, then we update count
                classifier_group['Total'][index] = count # updating the total count for the group. This is to get the total number of h-bons per group

        marker = itertools.cycle((',', '+', '.', 'o', '*')) # to get markers using an iterator
        NUM_COLORS = len(classifier_group.keys()) # get the different colors neeed
        x_group = [i * a_frameTime for i in range(len(self.allBonds_list))] # getting the x-axis values. a_frameTime is a constant provided to the function

        # -----------------------------------------------------
        # Plotting the h-bonds for the groups
        # -----------------------------------------------------

        print('\nPlotting the summary of h-bonds per frame\n')
        ax = plt.figure().gca()
        cm = plt.get_cmap('terrain') # slecting the matplolib color set
        ax.set_prop_cycle('color', [cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)]) # setting the cycle to select the colors per line in the plot

        for group in classifier_group.keys():

            if group != 'Total':
                y_group = classifier_group[group][group]
                plt.plot(x_group, y_group, label=group, marker=next(marker), markersize=3)
            else:
                y_group = classifier_group[group]
                plt.plot(x_group, y_group, label=group, marker=next(marker), markersize=3)

        plotname = self.plotname + '_summary.png'
        plot_parameters(plotname) # function implemented at the begining of this function. It is just a bunch of matplotlib constansts.

        # -----------------------------------------------------
        # Plotting the h-bonds different h-bonds inside each group
        # -----------------------------------------------------

        print('\nPlotting the h-bonds per group in all frames\n')

        for key_group, group in classifier_group.items():

            if key_group == 'Total':
                continue

            ax = plt.figure().gca()
            cm = plt.get_cmap('terrain')
            NUM_COLORS = len(group.keys())
            ax.set_prop_cycle('color', [cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)])

            for key, y_subset in group.items():
                plt.plot(x_group, y_subset, label=key, marker=next(marker), markersize=3)

            plotname = key_group + '_' + self.plotname + '.png'
            plot_parameters(plotname) # function implemented at the begining of this function. It is just a bunch of matplotlib constansts.

        return None


# ==================================================================================================
# XvgPlotter
# ==================================================================================================
# This class is only to plot h_bond information.

class XvgPlotter(PlotManager):
    '''This class is only to plotting h_bond information.'''

    def __init__(self, x_values, features, title, xlabel, legends, ylabel, plot_filename):
        self.x_values = x_values
        self.features = features
        self.title = title
        self.xlabel = xlabel
        self.legends = legends
        self.ylabel = ylabel
        self.plot_filename = plot_filename

# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

    @classmethod
    def from_XVGfile(cls, a_destinationFilePath, a_fileName):
        '''
        XvgPlotter.from_XVGfile is a class constructor that takes the 2 arguments:
        a_destinationFilePath: the location where the xvg file is located.
        a_fileName: the name of the xvg file.
        '''

        file = FileUtil.FileReader(a_destinationFilePath, a_fileName).file
        legends = []

        for line in file:
            if line.find('title') >= 0: # Getting the title
                title = line.split('"')[1]
            elif line.find('xaxis') >= 0: # getting x-axis label
                xlabel = line.split('"')[1]
            elif line.find('yaxis') >= 0: # getting y-axis label
                ylabel = line.split('"')[1]
            elif line.find('@ s') >= 0: # getting the legends
                legends.append(line.split('"')[1])

        # Getting the dats to plot
        #-------------------------
        data = [line for line in file if (line.find('#') < 0 and line.find('@') < 0)]
        line_count = len(data)
        feature_num = len(data[0].split()) - 1  # number of features variables (x-axis is the first column so we deduct 1 unit)
        x_values = np.zeros(line_count, dtype=float)
        features = np.zeros((feature_num, line_count), dtype=float).T

        for index, line in enumerate(data):
            parts = line.split()
            x_values[index] = round(float(parts[0]), 4)

            for var in range(feature_num):
                new_part = var + 1 # we start from 1 beacause zero is x_values
                features[index][var] = round(float(parts[new_part]), 4)

        features = features.T
        plot_filename = a_fileName.split('.')[0] + '.png'

        return cls(x_values, features, title, xlabel, legends, ylabel, plot_filename)


# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

    def featureset_scaling(self):
        '''
        XvgPlotter.featureset_scaling is a method to scale de density values to be between 0 and 1.
        This is needed because some molecules have a higher density than others, and also the h-bond
        density values will be less than the molecules.
        Scaling all from 0 to 1, we are sure that we will be able to see the b-bonds location better.
        '''

        for i in range(len(self.features)):
            self.features[i] = self.features[i] / (self.features[i].max(axis=0) + np.spacing(0))

        return None


# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

    def plot(self):
        # Density profile information has always many zero values at the begining and at the end
        # of the data. We will only select the values are are higher than zero to concentrate
        # the plotting are to the part of interest.

        # Getting the highest index in the entire featureset that has a non-zero value
        # ------------------------------
        # This is to determine the right limint in plot.xlim
        highest_nonzero_index = 0
        for i in range(len(self.features)):
            highest_nonzero_i = np.where(self.features[i] != 0)[0][-1]

            if highest_nonzero_index < highest_nonzero_i:
                highest_nonzero_index = highest_nonzero_i

        xlim_right = math.ceil(self.x_values[highest_nonzero_index]) # rounding to the ones place (see python documentation for more info)

        # Getting the lowest index in the entire featureset that has a non-zero value
        # ------------------------------
        # This is to determine the left limint in plot.xlim
        lowest_nonzero_index = highest_nonzero_index
        for i in range(len(self.features)):
            lowest_nonzero_i = np.where(self.features[i] != 0)[0][0]

            if lowest_nonzero_index > lowest_nonzero_i:
                lowest_nonzero_index = lowest_nonzero_i

        xlim_left = int(self.x_values[lowest_nonzero_index]) # rounding to the int value


        # Plotting the information
        # ------------------------------

        fig, ax = plt.subplots()

        for var in range(len(self.features)):
            plt.plot(self.x_values, self.features[var], label=self.legends[var])

        plt.legend(loc='best', fontsize=9)
        plt.xlabel(self.xlabel)
        plt.ylabel(self.ylabel)
        plt.xlim(left=xlim_left, right=xlim_right)
        plt.grid(axis='y', linestyle='--', linewidth=0.5)
        savepath = os.path.join(os.getcwd(), self.plot_filename)
        plt.savefig(savepath, dpi=1800, bbox_inches="tight")

        return None



if __name__ == "__main__":
    xbgPlotter_Obj = XvgPlotter.from_XVGfile(os.getcwd(), 'example.xvg')
    xbgPlotter_Obj.featureset_scaling()
    xbgPlotter_Obj.plot()























