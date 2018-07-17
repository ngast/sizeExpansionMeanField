# Size Expansions of Mean Field Approximation : Transient and Steady-State Analysis

This repository contains what is needed to reproduce the paper "Size Expansions of Mean Field Approximation : Transient and Steady-State Analysis" accepted at the conference "Performance Evaluation 2018" to be held in Toulouse (France).

Authors :
* Nicolas Gast (Inria)
* Luca Bortolussi (University of Trieste)
* Mirco Tribastone (IMT Lucca)

The repository includes the code to reproduce all simulations. 

## Receipe

### How to make the figures

Typing "make figures" generates all figures for the paper. This requires to have python3 installed as well as jupyter-notebook.

In order to faster the generation of the above figures, many
simulations of the various stochastic models have been precomputed
offline. You can remove them by typing "make
remove_precomputed_simulations" in order to force "make figure" to
regenerate them. However this might take a while (probably several hours). 

### How to compile the paper 

Type "make". This generates the figures and compiles the paper. 
