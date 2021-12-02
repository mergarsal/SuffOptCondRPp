# Sufficient optimality condition for the relative pose between calibrated cameras
Sufficient optimality condition 
for the relative pose problem 
between two calibrated cameras. 
Refer to our paper [HERE](https://doi.org/10.1137/21M1397970) 
for more information.



**Authors:** 
[Mercedes Garcia-Salguero](https://mapir.isa.uma.es/mapirwebsite_wordpress/?p=1718), 
[Javier Gonzalez-Jimenez](https://mapir.isa.uma.es/mapirwebsite_wordpress/?p=1536)


**License:** [GNUv3](https://github.com/mergarsal/SuffOptCondRPp/blob/main/LICENSE)


If you use this code for your research, please cite:

```
@ARTICLE{,
    author = {Garcia-Salguero, Mercedes and Gonzalez-Jimenez, Javier},
     month = {jan},
     title = {A Sufficient Condition of Optimality for the Relative Pose Problem between Cameras},
   journal = {SIAM Journal on Imaging Sciences},
    volume = {14},
    number = {4},
      year = {2021},
       url = {https://doi.org/10.1137/21M1397970},
       doi = {10.1137/21M1397970},
     pages = {1617--1634}
}
```



# Dependences

This repository include *manopt* as submodule,
so you need to clone it with the `--recursive` option

If you have already cloned it, you can still set the submodules with
```
git submodule update --init --recursive
```

*manopt*: 
    web: https://www.manopt.org
    github: https://github.com/NicolasBoumal/manopt
    

# Example
1. Run ``setup.m`` to import all the dependences
2. Run the basic example: ``example_suff_cond`` 

# Contents
This repo contains the following implementations: 
1. On-manifold estimation of the essential matrix.
    certifiers/solveRPpManifold.m
    
2. Generic closed-form certifiers:
    certifiers/checkOptimalityCertifier
    certifiers/computeLagrangeMultipliersReducedRelaxation.m
    
3. Sufficient condition in folder *condition*

4. Some useful functions in folder *utils*
