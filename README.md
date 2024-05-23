# ZircSat

[![Build Status](https://github.com/lcandiot/ZircSat.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/lcandiot/ZircSat.jl/actions/workflows/CI.yml?query=branch%3Amain)

Calculate the liquidus temperature and estimate the Zircon saturation temperature based on thermodynamic models for arc magmas.

🚧 Note: this is a preliminary version. 🚧

This package is built on [MAGEMin](https://github.com/ComputationalThermodynamics/MAGEMin), a Gibbs energy minimization tool. Calculation of the liquidus temperature, i.e. the temperature at which the solid fraction is near 0, is performed using Newton-Raphson iterations. This tool is built entirely in the [Julia](https://julialang.org) language.

## Installation

First, install Julia on your system following the official [Julia download instructions](https://julialang.org/downloads/) for your system. Upon successful installation of Julia, execute the following steps to install and use ZircSat:

1. Open a terminal and clone this repository to any destination on your machine

   ```
   git clone https://github.com/lcandiot/ZircSat.jl.git
   ```
2. Step into the newly created repository

   ```
   cd ZircSat.jl
   ```
3. Create a personal working directory named `/USER`

   ```
   mkdir -p USER
   ```
4. Copy the `/examples/ex_calculate_Tliq_Tsat.jl` and ` /data/ZircSat_test_MarxerUlmer2019.csv` files into your newly generated  `/USER `directory. Rename `/examples/ex_calculate_Tliq_Tsat.jl` to `my_Tliq_Tsat_calculator.jl` and `/data/ZircSat_test_MarxerUlmer2019.csv` to `/data/my_data.csv`

   ```
   cp -r examples/ex_calculate_Tliq_Tsat.jl data/ZircSat_test_MarxerUlmer2019.csv USER/
   mv USER/ex_calculate_Tliq_Tsat.jl USER/my_Tliq_Tsat_calculator.jl
   mv USER/ZircSat_test_MarxerUlmer2019.csv USER/my_data.csv
   ```
5. Start the Julia REPL

   ```
   julia --project
   ```
6. To test if the installation works activate the package mode and test the ZircSat package

   ```
   ]
   test
   ```
7. Exit the package environment (hitting backspace on your keyboard).
8. Open `/USER/my_Tliq_Tsat_calculator.jl` with a code editor of your choice and modifiy the data path line

   ```
   data_path = "./USER/my_data.csv"
   ```
9. Run `/USER/my_Tliq_Tsat_calculator.jl` by entering the following command in the Julia REPL

   ```
   include("USER/my_Tliq_Tsat_calculator.jl")
   ```
10. Congratulations! You can now modify these two files according to your needs and add your own data to `/USER/my_data.csv`


## Input / Output

The `calculate_liquidus_temperature()` function expects a data table stored as `.csv` file. This table should contain the test pressure in kbar as well as the major oxides `SiO2-Al2O3-CaO-MgO-FeO-Fe2O3-TiO2-K2O-Na2O-Cr2O3-H2O` which are required by MAGEMins igneous database. For now only this database has been tested. While the order of the columns does not matter, the package is case sensitive. It is therefore recommended to use the [ZircSat_test_MarxerUlmer2019.csv](https://github.com/lcandiot/ZircSat.jl/tree/main/data/ZircSat_test_MarxerUlmer2019.csv) as a template. Further, the Zirconium concentration in units ppm is required as `Zr` column in the data file.

Upon calculation of the liquidus temperature, a new `.csv` file will be written to the same location as the input file. The new file contains the original data plus an additional column in which the calculated liquidus temperature is stored. Although the precision of this prediction is $\approx$ 1 °C for the test cases, the accuracy of the predicted liquidus temperature is dependent on the accuracy of MAGEMin predicting the stable mineral phases correctly.

## Working with this repository

It is recommended to generate a `/USER` directory in the top level of this repository and store or files which are not meant to be published there. This directory should also be the place to work from and run scripts. Feel free to use the [examples](https://github.com/lcandiot/ZircSat.jl/tree/main/examples) as starting point for your own routines.
