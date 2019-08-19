# DiSolver.py

The following is a brief description of the input/output structure. For a full description of the code see our paper: LINK.


## Structure of input files

### in.csv
This file contains the main input parameters as well as switches for turning automated calibration and conversion of FEC on/off.
Do not delete any lines from the file.

![In file 1](/Images/incsv.PNG)

### flows.csv
This file contains the depths and flow rates of the fractures. If used in forward mode, the two columns *Upper depth limit* and *Lower depth limit* are not required. Positive flows are inflows and negative flows are outflows. If using the code for FFEC logging (i.e. a pumped borehole), add the depth and flow rate of the pump to this file, as you would a fracture.

![In file 2](/Images/flowscsv.PNG)

### initialcondition.csv
This file contains the depth vs. concentration data for the initial state. DiSolve.py will interpolate the data onto a 1D grid with the spatial discretization given in *in.csv*. If these data are FEC data rather than concentration data (e.g. in kg/m<sup>3</sup>), set "*Initial condition given as concentration (1) or as FEC (2)*" to 2 in *in.csv*. 

![In file 3](/Images/incon.PNG)

### measuredprofiles.csv

This file contains the measured profiles at the output times defined in *in.csv*. This is a required input for inversion modelling but optional for forward modelling.

![In file 4](/Images/measured.PNG)

## Output

The modelled depth vs. concentration data can be found in *Output\profiles.csv*. If used in inversion mode, *Output.csv* will contain the optimized output parameters.

The modelled data and measured data (if given) are plotted in *profiles.csv*.

![In file 5](/Output/profiles.png)