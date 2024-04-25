# Percolative supercurrent in superconductor--ferromagnetic-insulator bilayers

This repository contains code and data for the paper _Percolative supercurrent in superconductor--ferromagnetic-insulator bilayers_.

It is structured as follows: 
- `exp_data` contains the data collected in the measurements.  
- `bilayers_fittings` contains the code and results for the fitting procedure of the NIS data with the Usadel model.
- `network_model` contains the resistive network model for zero-bias conductance and resistivity.

## Experimental data

### NIS data
The `NIS_data` folder contains the tunneling spectroscopy data of the junctions. 
- The `MDC004.Sample1` subfolder contains data for _Sample 1_ in the manuscript.
  
  |Filename| Description                   |
  |--------|-------------------------------|
  |`15.dat`|60x50  μm$`^2`$,   downsweep   |
  |`16.dat`|60x50  μm$`^2`$,   upsweep     |
  |`24.dat`|6x100  μm$`^2`$,   downsweep   |
  |`25.dat`|6x100  μm$`^2`$,   upsweep     |
  |`39.dat`|6x20  μm$`^2`$,   downsweep    |
  |`40.dat`|6x20  μm$`^2`$,   upsweep      |
  |`51.dat`|25x5  μm$`^2`$,   upsweep      |
  |`68.dat`|15x50  μm$`^2`$,   upsweep     |
  |`71.dat`|60x20  μm$`^2`$,   upsweep     |
  |`77.dat`|60x100  μm$`^2`$,   upsweep    |
  |`79.dat`|25x100  μm$`^2`$,   upsweep    |
  |`81.dat`|15x100  μm$`^2`$,   upsweep    |
  |`83.dat`|60x6  μm$`^2`$,   upsweep      |
  |`91.dat`|25x50  μm$`^2`$,   upsweep     |
  |`93.dat`|6x50  μm$`^2`$,   upsweep      |
  |`95.dat`|25x20  μm$`^2`$,   upsweep     |
  |`97.dat`|15x20  μm$`^2`$,   upsweep     |
  

- The `MDC004.Sample2` subfolder contains data for _Sample 2_ in the manuscript.
  
  |Filename| Description                   |
  |--------|-------------------------------|
  |`28.dat`|10x10  μm$`^2`$,   downsweep   |
  |`29.dat`|10x10  μm$`^2`$,   upsweep     |
  |`41.dat`|10x20  μm$`^2`$,   upsweep     |
  |`45.dat`|6x100  μm$`^2`$,   upsweep     |
  |`47.dat`|5x5  μm$`^2`$,   upsweep       |
  |`49.dat`|10x5  μm$`^2`$,   upsweep      |
  |`96.dat`|10x10  μm$`^2`$,   downsweep   |
  |`97.dat`|25x50  μm$`^2`$,   upsweep     |
  |`105.dat`|5x10  μm$`^2`$,   upsweep     |
  

### BAR data
The `BAR_data` folder contains transport data referring to the four-terminal resistance measurements of the Al bars.

  |Filename| Description                                     |
  |--------|-------------------------------------------------|
  |`69.dat`|Sample 2: 2, 5, 10, 20, and 50 μm, downsweep     |
  |`70.dat`|Sample 2: 2, 5, 10, 20, and 50 μm, upsweep       |
  |`197.dat`|Sample 1: 5, 10, and 20  μm, upsweep            |
  |`198.dat`|Sample 1: 5, 10, and 20  μm, downsweep          |
  |`202.dat`|Sample 1: 2 and 40  μm, upsweep                 |
  |`203.dat`|Sample 1: 2 and 40  μm, downsweep               |
  

### AFM data
The folder `AFM_data` contains the Atomic Force Microscopy data shown in the Supplementary Material.

### EDS data
The folder `EDS_data` contains the energy-dispersive X-ray spectroscopy data shown in the Supplementary Material.

## Bilayers fittings
The `bilayers_fitting` folder contains code and the output of the Usadel fittings of the NIS spectroscopy curves.
- `map_fit_minimal.ipynb` is the code used for the minimal fitting (used for Sample 1 in the paper).
- `map_fit_complete.ipynb` is the code used for the complete fitting (used for Sample 2 in the paper).
- `gen_theory_data.ipynb` is the code that uses the fitting results to solve the Usadel model and show the predicted NIS tunneling conductance.

## Network model
The `network model` folder contains the code and the output of the resistive network model for the zero-bias conductance and the model for the percolative supercurrent.
- `network_model.py` is the code of the model.
- `network_model_Sample1.ipynb`is the code for Sample 1.  
- `network_model_Sample2.ipynb`is the code for Sample 2.
- `network_model_general.ipynb` can be used to run the network model for zero-bias conductance with arbitrary parameters.
