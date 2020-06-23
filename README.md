# membrane_distillation
1D and Quasi-2D modeling tool for Direct Contact Membrane Distillation

The code here implements a 1D model of pore-scale transport processes in a Direct Contact Membrane Distillation (DCMD) water filtration system.

The system of equations integrates temporal conservation equations for mass, element, and energy, over a sufficiently long time to estimate steady-state behavior as a funciton of feed temperature and distillate temperature.

## Installing

Installing via [Conda](https://conda.io/) is highly recommended.  

The simulation requires a number of Python dependencies, including [Cantera](https://cantera.org), [numpy](http://www.numpy.org/), [ruamel.yaml](https://pypi.org/project/ruamel.yaml/), and [scipy](https://www.scipy.org/).  For all but Cantera, installation can be done using simple [Conda](https://conda.io/) commands:

> ```conda install <package name>```
  
For Cantera installation, please see [Cantera's installation instructions](https://cantera.org/install/index.html).

If you are using Conda, I would recommend setting up a single environment with all required packages (it is strongly advised that you do not install packages into your `base` conda environment):

> ```conda create --name membrane --channel cantera/label/dev cantera numpy scipy matplotlib pandas```

This creates an environment `modeling` from which you can run the model (you do not have to call it modeling). Activate this environment by executing the command:

> ```conda activate modeling```

The current version of the code also requires the ruamel.yaml module. Once you have activated the environment (`modeling`, or whatever name you gave it), install `ruamel.yaml` to this environment, via Conda:

> ```conda install -channel conda-forge ruamel.yaml```

(the `-channel` flag indicates to download from a specific conda channel.  Above, we downloaded and installed from the `cantera/label/dev` channel, to get the development version of cantera.  Here, we download and install from `conda-forge`, which hosts the `ruamel.yaml` package).

Once you have your Conda environment set up, you need to download the model files. If you use `git`, this is relatively simple to do via cloning from the command line:

> ```git clone https://github.com/coresresearch/membrane_distillation.git```

You can also click "Clone or download," above, followed by "Download ZIP".

Next cd into the folder where you downloaded the model:

> ```cd membrane_distillation```

## Running the Model

Once you have installed all necessary dependencies and downloaded the model files, there are two necessary steps:

1. Specify inputs in `membrane_distillation_inputs.yaml`
2. Call/run the model file `membrane_distillation_1d_model.py`

### Specifying inputs.

All necessary inputs are specified in the file `membrane_distillation_inputs.yaml`.  This includes all simulation and membrane parameters, as well as the cantera thermo-kinetic input parameters.  

Selected parameters can also be provided when the model is called:

* `save`: The folder and file where outputs are saved 
* `membrane`: The set of membrane properties in `membrane_distillation_inputs.yaml` to use for the simulation (the user must provide a string matching one of the name strings in the input file).
* `tau_g`: The membrane tortuosity factor.
* `transport`: The transport model (a string matching one of the implemented models. Current options are `Fick` for linear Fickean diffusion and `DGM` for the Dusty Gas Model).
* `feed_temp`: The membrane surface temperature at the (hot) feed water interface [degrees C].
* `permeate_temp`: The membrane surface temperature at the (cold) permeate water interface [degrees C].

In the case where certain parameters are given in both the input file and from the command line, _only the values from the command line will be used_.

### Run the Model File

There are two ways to run the model file.  

1. If you just want to run a single simulation for a given set of conditions, you can run `membrane_distillation_1d_model.py`, either from the command line or from an IDE such as `PyCharm`, `Spyder`, etc.  

    To use inputs from the input file only, simply call the model file.  From the command line, this looks like:  

    > ```python membrane_distillation_1d_model.py```

    If you want to overwrite some of the paramters in the input file, these can be provided via command line keywords.  For example, to save to a folder/file named `outputs/test.csv` (all output files are automatically written to the `outputs` folder), for the `DGM` transport model and a tortuosity factor of 1.3, you would run:

     ```python
     python membrane_distillation_1d_model.py --save='test' --transport='DGM' --tau_g=1.3
    ```

    (from within an IDE, you would mostly likely drop the `python` prefix and `.py` file suffix in these commands).

2. Running fits or repeated models:  If you want to run a series of models, or perform a quasi-2D simulation for a given surface temperature profile for each side of the membrane, `membrane_distillation_1d_model.py` can be imported as a module, and run from a python script file.  To run the same model as shown in the block above:

```python
from membrane_distillation_1d_model import membrane_distillation_1d_model
membrane_distillation_1d_model(save='test', transport='DGM', tau_g=1.3)
```

In this way, your script can easily loop over numerous simulations. The data is saved as you go, which you can then process and analyze after running the simulations. A demonstration of such a fitting routine under construction in the file `membrane_distillation_1d_fitting.py`.

## Using, Citing, and Contributing Code.

This code is licensed under a permissive, [BSD 3-clause license](https://github.com/coresresearch/membrane_distillation/blob/main/LICENSE).  You are free to re-purpose, edit, and use this model as you see fit, so long as you follow the license guidelines.

If you publish anything based on this work, we would appreciate using the following citatioon: 

```
@misc{CORES_membrane_distillation
    author = {Spencer D. Gilleon and Steven C. DeCaluwe},
    year = 2020,
    title = {CORES Research Group Membrane Distillation Model}
    url = {https://github.com/coresresearch/membrane_distillation},
}
```
Finally, if there are changes, improvements, and/or suggestions, please let us know!  You can open an issue, using the "Issues" tab above.  Or better yet, fork and clone this repository, commit your changes using Git, and open up a Pull Request!  We would gladly welcome your contributions!


