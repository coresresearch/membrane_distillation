# membrane_distillation
1D and Quasi-2D modeling tool for Direct Contact Membrane Distillation

The code here implements a 1D model of pore-scale transport processes in a Direct Contact Membrane Distillation (DCMD) water filtration system.

The system of equations integrates temporal conservation equations for total mass, element, and energy, over a sufficiently long time to estimate steady-state behavior as a funciton of feed temperature and distillate temperature.

All necessary inputs are specified in the header of `membrane_runner.py`.  This is also the file you run to carry out the simulation.
After running the simulation, copy the relevant inputs to the header of `MD_sim_plotting.py` and run that file to plot the predicted flux as a function of temperatures and tortuosity, and compare to co-current bench-scale experiments.
The model can take the nominal feed and distillate temperatures as the boundary conditions, or can perform a quasi-2D simulation if given a surface temperature profile for each side of the membrane.  By way of example, a surface profile for feed and distillate inlet temperatures of 30 and 20 C, respectively, is provided.

The simulation requires a number of Python dependencies, including [Cantera](https://cantera.org), [numpy](http://www.numpy.org/), [pandas](https://pandas.pydata.org/), and [scipy](https://www.scipy.org/).  For all but Cantera, installation can be done using simple [Conda](https://conda.io/) commands:

  ```conda install <package name>```
  
For Cantera installation, please see [Cantera's installation instructions](https://cantera.org/install/index.html).

If you are using Conda, I would recommend setting up a single environment with all required packages:

```conda create --name modeling --channel cantera/label/dev cantera numpy scipy matplotlib pandas```

This creates an environment `modeling` from which you can run the model.  
