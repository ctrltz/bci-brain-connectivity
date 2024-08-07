## BCI & Brain Connectivity

This repository contains the code for the multiverse analysis of SNR and EEG-based functional connectivity in the context of sensorimotor BCI training. Please refer to the accompanying paper for details and cite it if you use the code from this repository:

> Kapralov, N., Jamshidi Idaji, M., Stephani, T., Studenova, A., Vidaurre, C., Ros., T., Villringer, A., Nikulin, V., 2023. Sensorimotor brain-computer interface performance depends on signal-to-noise ratio but not connectivity of the mu rhythm in a multiverse analysis of longitudinal data. bioRxiv. doi:[10.1101/2023.09.30.558407](https://doi.org/10.1101/2023.09.30.558407)

### Prerequisites

1. MATLAB (tested with R2022b) and the following toolboxes/files (can also be installed in the `toolboxes` folder using the `make setup_toolboxes` command in Linux):

 * [[link](https://github.com/bbci/bbci_public)] Berlin Brain-Computer Interface (BBCI) toolbox. **NOTE**: `toolboxes/BTB/BTB.mat` is provided for the BBCI toolbox to run. 
 * [[link](https://sccn.ucsd.edu/eeglab/download.php)] EEGLAB 2021.0. **NOTE**: The toolbox was manually set to [use double precision](https://eeglab.org/tutorials/misc/EEGLAB_option_menu.html) in order to avoid problems when analyzing data outside of EEGLAB.
 * [[link](https://research.ics.aalto.fi/ica/fastica/code/FastICA_2.5.zip)] FastICA v2.5 (only required if performing preprocessing from scratch)
 * [[link](https://www.parralab.org/nyhead/sa_nyhead.mat)] Pre-computed New York Head model
 * [[link](https://github.com/fooof-tools/fooof_mat)] MATLAB Wrapper for FOOOF 
 * [[link](https://github.com/jadref/tprod)] tprod toolbox. **NOTE**: apparently, it is problematic to compile this toolbox in the newer versions of MATLAB - one could either do it in previous versions or try to use the compiled files from this repository.
 
2. Python with FOOOF installed (`pip install fooof` or use `requirements.txt`). **NOTE**: keep in mind the [compatibility](https://www.mathworks.com/support/requirements/python-compatibility.html) between different versions of MATLAB and Python.

3. R (tested with version 4.2.2) and a number of packages specified in the environment file (`stats/renv.lock`). `renv::restore()` can be used to install the required packages, more details are available [here](https://rstudio.github.io/renv/articles/collaborating.html).

4. Local LaTeX environment for compiling the final PDF of the paper and supplementary material.

5. [Optional] For the `Makefile` to be useful, a Linux-based OS is required with the following utilities: `wget`, `unzip`, `curl`.

### Folder Structure

* `assets` - images that are included in some of the figures. Some of them are generated through MATLAB scripts
* `data` - folder with all the data (added to .gitignore)
  * `raw` - original dataset
  * `aux` - preprocessing info
  * `preproc`
    * `task1` - preprocessed data
  * `derivatives` - intermediate results
* `paper` - everything that is required for LaTeX
  * `figures` - here the generated figures are collected
  * `numbers` - here the generated output is collected
  * `sections` - sections of the manuscript
* `results` - figures & output for TeX that is generated by all scripts (added to .gitignore)
* `scripts` - MATLAB scripts
* `stats` - R scripts for statistical analysis
* `toolboxes` - folder with MATLAB toolboxes

### Creating a Local Copy of the Project

1. Clone the `master` branch of the repository.

2. Set up all the prerequisites (see above).

3. Download the original dataset using `make download_raw_data` or create a symbolic link to the location of the data if it was downloaded earlier (scripts in this repository do not change the original files):

```
ln -s <location> ./data/raw
```

4. Download the preprocessing information using `make download_aux_data` or manually from the [OSF repository](https://osf.io/tcvyd/).

5. Specify the path to a Python executable with FOOOF installed in line 7 of the file `scripts/BCI_MI_analysis_main.m`.

6. Specify the commands for running MATLAB and R in the lines 13 (MATLAB_ALIAS) and 15 (R_ALIAS) of the `Makefile`, respectively.

7. Launch the analysis using `make all` on Linux or follow the steps specified in the `Makefile`. Files `scripts/BCI_MI_analysis_main.m` and `stats/main.R` are the main entry points for the analyses performed in MATLAB and R, respectively.

8. [Optional] Run `git status`. If there are no changes detected, the results were repeated.
