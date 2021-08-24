# UltraFrost
UltraFrost is developed for pore-scale quantitative characterization of permafrost samples using ultrasonic waves 

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5159712.svg)](https://doi.org/10.5281/zenodo.5159712)
 
## Dependencies 
The hybrid inverse and multi-phase poromechanical transfer function solver (UltraFrost) is implemented in Fortran and Python. The attached wrapper (trans_vp.so) for Fortran code can only run in Linux system. The following dependencies are required in order to run WaveFrost : 

gfortran compiler:
```
$ sudo apt install gfortran-9
```

blas and lapack library: 

```
$ sudo apt-get install libblas-dev checkinstall 
$ sudo apt-get install libblas-doc checkinstall 
$ sudo apt-get install liblapacke-dev checkinstall 
$ sudo apt-get install liblapack-doc checkinstall
```

## Instruction
An example code (inversion.py) is given to show the estimation of properties of a permafrost sample. Users can simply run following command to load the lab measurement, define loss function and run inversion analysis for the estimation of the permafrost or frozen soil properties. 
```
$ python inversion.py
```

An example experimental transfer function data for permafrost samples can be found in the 'data' folder. To use a different dataset, users can change the file name accordingly in line 8 of inversion.py. 

```
da = np.loadtxt('data/measurement5.csv') 
```

The aim function (los_func) is defined based transfer function and first kind of P wave velocity. 

```
loss = np.abs(yt - y_data)
los_func = np.abs(vp-vp_m)+ 1000*np.sqrt(np.sum(loss**2))
```

The inversion analysis is performed based on the Neighborhood algorithm and detailed information can be found in https://github.com/keithfma/neighborhood. Users should run the 'inversion.py' multiple times with various random initial values to increase the robustness of the inversion results. 

