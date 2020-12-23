# Hydrogen-Bonds-calculations-for-Gromac-MD-simulations

@Author: Andres Vodopivec

@Date: 2020-12-08

This program can be used to calculate quantify and plot the location of hydrogen bonds in any Molecular Dynamic's system compatible with Gromacs5.x or higher, independent of the the force-field used, as long as it is an Atomistic force-field.

The purpose of this program is to quantify and plot the hydrogen bond interaction between two different species of molecules.
The two different species of molecules are referred as: molec0 and molec1.

### PROGRAM REQUIREMENTS:

**a)** This program requires: FileManager.py, Topology.py, Additionals.py, Molecule.py, H_bond.py, PlotManger.py, and INPUT_USER_H_BOND_CALCS.json.

All entries must be made in the INPUT_USER_H_BOND_CALCS.json file.

**b)** Entries can be made to the command line. Please use ```python h_bonds_calc.py --help``` for usage information.

**c)**  The GRO files to calculate the hydrogen bonds. This grofiles should be frames from a simulation.

The following command line can be used to extract your simulations frames:
```
gmx_mpi trjconv -sep -o frame.gro -skip 1 -f sim.xtc -s sim.tpr -pbc whole. 
```

GMX will generate a series of frameX.gro file, where X is a number from 0 to how many frames are present in your xtc file. 

If a list is supplied using the -g flag in the python commnand line, it is important that your GRO files are named using alpha-numeric (i.e. mygro0.gro, mygro1.gro, mygro2.gro, mygro3.gro, ....); otherwise, the code will fail.
    

**d)**  The two itp files (one for each of the two molecule species). We refer in them in the program as molec0_itp and molec1_itp.

**e)** The xtc and tpr files from which the GRO files were extracted

**f)** All the files mentioned in this section MUST be located in the same directory where the Main_hbonds.py is located.


### UNDERSTANDING PROGRAM NOTATION AND CONSIDERATIONS:

To consider a hydrogen bond, the following considerations must be taken into account:

a)  The hydrogen from a molecule species has to be chemically bonded to another atom in the same molecule. This is represented by X-H. Where H is the hydrogen atom, X is the atom chemically bonded to the hydrogen atom, and "-" represents the chemical bond.

b)  The third atom has to come the other molecule species, and it is represented by Y. Therefore, the notation X-H---Y is the Y atom interacting with the hydrogen from the other molecule species, and "---" represent the dipole-dipole interaction.

c)  The Y atom must be within a maximum distance threshold from the Hydrogen atom only (designated as dist_thres).

d)  The configuration X-H---Y should comply with the cosine rule. In this cosine rule we need all there distances (XtoH, XtoY, and HtoY), and the angle that Y forms with X-H should have a minimum threshold (designated as angle_thres) in order to be considered a hydrogen bond.

e)  Hydrogen bonds cannot take place if Y atom is an Alkali (i.e. sodium, potassium, or calcium) or another hydrogen atom.


### IMPORTANT:

a)  All files must be compatible with Gromacs in order for this program to work.

b)  Please be sure to read carefully all the help info and comments to clear any doubts.

c)  The hydrogen bond criteria in this program are supported in the publication:

    Experimental and Computational Studies of Choline Chloride-Based Deep Eutectic Solvents
    Sasha L. Perkins, Paul Painter, and Coray M. Colina
    Journal of Chemical & Engineering Data 2014 59 (11), 3652-3662
    DOI: 10.1021/je500520h


### PROGRAM ADDITIONALS:

This program also calculates and plot a density profile plot between molec0, molec1, and the locaiton of the hydrogen bonds. Gromacs must be installed in order to determine the density profile.

If there is no interest in doing this extra analysis, just uncomment lines from the sections: 'Getting the density profile using gmx density' and 'Doing the density profile plot'.


### USEFUL TIP:

If you do not want to carry so many files, you can create a python executable using pyinstaller. Make sure all the python files required are in the same directory and run:
```pyinstaller --onefile Main_hbonds.py ```

The executable will be located inside the dist folder. This excutable is all you need and you can run it in any directory, just make sure you have the gromacs files and Json file in the same directory.

#### Please feel free to use this program. If you find anything useful in this program or decide to use it for your research, please fork it, click on the star, and add it in your references or citations.

