# DISOLV

The following is a brief description of the input/output structure. For a full description of the code see our paper: LINK.

## Installation

DISOLV can be installed for Python 2.7 or Python 3 using pip:

```
pip install disolv
```

It can be used within other Python code with:

```python
import disolv
disolv.run('input_dir', 'output_dir', calibrate=True, convertFEC=True,
           method='SLSQP')
```

## Structure of input files

### in.csv
This file contains the main input parameters. Do not delete any lines from the file.

![In file 1](/Images/incsv.PNG)

### flows.csv
This file contains the depths and flow rates of the fractures. If used in forward mode, the two columns *Upper depth limit* and *Lower depth limit* are not required. Positive flows are inflows and negative flows are outflows. If using the code for FFEC logging (i.e. a pumped borehole), add the depth and flow rate of the pump to this file, as you would a fracture.

![In file 2](/Images/flowscsv.PNG)

### initialcondition.csv
This file contains the depth vs. concentration data for the initial state. DISOLV will interpolate the data onto a 1D grid with the spatial discretization given in *in.csv*. If these data are FEC data rather than concentration data (e.g. in kg/m<sup>3</sup>), DISOLV can convert them to concentration data by setting the argument *convertFEC* to 'True'. 

![In file 3](/Images/incon.PNG)

### measuredprofiles.csv

This file contains the measured profiles at the output times defined in *in.csv*. This is a required input for inversion modelling but optional for forward modelling.

![In file 4](/Images/measured.PNG)

## Running the model

DISOLV can be imported and run as follows:

```python
    import disolv
    disolv.run('input_dir', 'output_dir', calibrate=False, convertFEC=False)
```   
    
The first and second arguments are the file paths to the input and output directories. *Calibrate* refers to whether the model is being run in forward ('False') or inverse ('True') mode, and *convertFEC* indicates whether the initial condition has been given in fluid electrical conductivity (μS cm<sup>−1</sup>) and must be converted to concentration (in kg m<sup>−3</sup>) (‘True’) or whether it has been given as a concentration (‘False’).

If DISOLV is run in inverse mode, the optimization method can be chosen in the final argument:

```python
    disolv.run('input_dir', 'output_dir', calibrate=True, convertFEC=False, method='SLSQP')
```
## Output

The modelled depth vs. concentration data can be found in *Output\profiles.csv*. If used in inversion mode, *Output.csv* will contain the optimized output parameters.

The modelled data and measured data (if given) are plotted in *profiles.csv*.

![In file 5](/Output/Berambadi.png)
