[![](https://readthedocs.org/projects/biobb-wf-pmx-tutorial/badge/?version=latest)](https://biobb-wf-pmx-tutorial.readthedocs.io/en/latest/?badge=latest)
[![](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/bioexcel/biobb_wf_pmx_tutorial/master?filepath=biobb_wf_pmx_tutorial%2Fnotebooks%2Fbiobb_wf_pmx_tutorial.ipynb)

# Mutation Free Energy Calculations using BioExcel Building Blocks (biobb)

***

**Based on the official [pmx tutorial](http://pmx.mpibpc.mpg.de/sardinia2018_tutorial1/index.html).**

***

This tutorial aims to illustrate how to compute a **fast-growth mutation free energy** calculation, step by step, using the BioExcel **Building Blocks library (biobb)**. The particular example used is the **Staphylococcal nuclease** protein (PDB code 1STN), a small, minimal protein, appropriate for a short tutorial.

The **non-equilibrium free energy calculation** protocol performs a **fast alchemical transition** in the direction **WT->Mut** and back **Mut->WT**. The two equilibrium trajectories needed for the tutorial, one for **Wild Type (WT)** and another for the **Mutated (Mut)** protein (Isoleucine 10 to Alanine -I10A-), have already been generated and are included in this example. We will name **WT as stateA** and **Mut as stateB**.

![](https://raw.githubusercontent.com/bioexcel/biobb_wf_pmx_tutorial/master/biobb_wf_pmx_tutorial/notebooks/schema.png)

The tutorial calculates the **free energy difference** in the folded state of a protein. Starting from **two 1ns-length independent equilibrium simulations** (WT and mutant), snapshots are selected to start **fast (50ps) transitions** driving the system in the **forward** (WT to mutant) and **reverse** (mutant to WT) directions, and the **work values** required to perform these transitions are collected. With these values, **Crooks Gaussian Intersection** (CGI), **Bennett Acceptance Ratio** (BAR) and **Jarzynski estimator** methods are used to calculate the **free energy difference** between the two states.

*Please note that for the sake of disk space this tutorial is using 1ns-length equilibrium trajectories, whereas in the [original example](http://pmx.mpibpc.mpg.de/sardinia2018_tutorial1/eq.mdp) the equilibrium trajectories used were obtained from 10ns-length simulations.*

***

## Settings

### Biobb modules used

* [biobb_pmx](https://github.com/bioexcel/biobb_pmx): Tools to setup and run Alchemical Free Energy calculations.
* [biobb_md](https://github.com/bioexcel/biobb_md): Tools to setup and run Molecular Dynamics simulations.
* [biobb_analysis](https://github.com/bioexcel/biobb_analysis): Tools to analyse Molecular Dynamics trajectories.

### Auxiliar libraries used

* [nb_conda_kernels](https://github.com/Anaconda-Platform/nb_conda_kernels): Enables a Jupyter Notebook or JupyterLab application in one conda environment to access kernels for Python, R, and other languages found in other environments.
* [os](https://docs.python.org/3/library/os.html): Python miscellaneous operating system interfaces
* [plotly](https://plot.ly/python/offline/): Python interactive graphing library integrated in Jupyter notebooks.

### Conda Installation and Launch

```console
git clone https://github.com/bioexcel/biobb_wf_pmx_tutorial.git
cd biobb_wf_pmx_tutorial
conda env create -f conda_env/environment.yml
conda activate biobb_wf_pmx_tutorial
conda install -y -c bioconda biobb_analysis==2.0.1
jupyter-notebook biobb_wf_pmx_tutorial/notebooks/biobb_wf_pmx_tutorial.ipynb
```

***

## Tutorial

Click here to [view tutorial in Read the Docs](https://biobb-wf-pmx-tutorial.readthedocs.io/en/latest/?badge=latest)

***

## Version
October 2019 Release

## Copyright & Licensing
This software has been developed in the [MMB group](http://mmb.irbbarcelona.org) at the [BSC](http://www.bsc.es/) & [IRB](https://www.irbbarcelona.org/) for the [European BioExcel](http://bioexcel.eu/), funded by the European Commission (EU H2020 [823830](http://cordis.europa.eu/projects/823830), EU H2020 [675728](http://cordis.europa.eu/projects/675728)).

* (c) 2015-2019 [Barcelona Supercomputing Center](https://www.bsc.es/)
* (c) 2015-2019 [Institute for Research in Biomedicine](https://www.irbbarcelona.org/)

Licensed under the
[Apache License 2.0](https://www.apache.org/licenses/LICENSE-2.0), see the file LICENSE for details.

![](https://bioexcel.eu/wp-content/uploads/2019/04/Bioexcell_logo_1080px_transp.png "Bioexcel")
