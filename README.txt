The script was written in MATLAB (2020b), version 2019a and higher are recommended for the correct use. It calculates the ring current geometrical factor (RCGF) for the provided molecule and defined ring current path. It can also estimates the ring current susceptibility I/B (in nA/T) if the NMR data are provided. 

It is executed by running the script run_this_script.m. 

Besides this, 6 more functions are included as separate scripts, which all need to be in the working directory (or the path of their locations must be specified). The functions are:

area_finder.m – this evaluates the net cross-section of the provided molecule in respect to the x,y,z axes
geometry_plot.m – provides simple visualization of the provided geometry
RCM.m – calculates RCGF
unit3d.m – contains evaluated integral from equation (S9) 
output_box.m – generates the prompts showing RCGF and asking the user for more instructions
fit_and_plot – fits a linear function to the RCGF and provided experimental NMR data 


To calculate the RCGF, three files from a user must be provided (the prompt window will ask to select the files):
     The geometry, in the form of the .xyz file (using conventional layout of the xyz files with first two lines containing the number of atoms and a title, both lines to be omitted).
     Specification of the atom numbers (relating to the numbering in the .xyz file) for which the RCGF should be calculated as a comma-delimited .csv file. Each column defines new type of the atoms. If more numbers are in one column, the RCGF is calculated for all of them and then averaged for the column. 
     Connectivity matrix specifying from which segments is the RCM constructed as a comma-delimited .csv file. Every row corresponds to the one segment, one bond. The row contains three numbers: weight of the ring current (= 1 if the only option for the current to flow), number of the atom from which the segment starts and number of the atom at which the segment ends. 

For the fit of the linear function between the RCGF and experimental data, one more file containing the experimental chemical shift differences (aromatic-nonaromatic) needs to be provided. This is a simple single column .csv file with numbers on the n-th row corresponding to the n-th column in the file specifying the atom numbers for which the RCGF is calculated. 



If used in publication please cite as: doi:10.xxxx/xxxxxxx

GitHub page:
https://github.com/mjirasek

Questions regarding the script (how to run / found bugs) adress to m.jirasek@seznam.cz
Alternatively, Prof. Harry Anderson can be reached on harry.anderson@chem.ox.ac.uk

Michael Jirasek, PhD.
