# ZircSat

[![Build Status](https://github.com/lcandiot/ZircSat.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/lcandiot/ZircSat.jl/actions/workflows/CI.yml?query=branch%3Amain)

Calculate the liquidus temperature and estimate the Zircon saturation temperature based on thermodynamic models for arc magmas.

🚧 Note: this is a preliminary version. 🚧

This package is built on [MAGEMin](https://github.com/ComputationalThermodynamics/MAGEMin), a Gibbs energy minimization tool. Calculation of the liquidus temperature, i.e. the temperature at which the solid fraction is near 0, is performed using Newton-Raphson iterations. This tool is built entirely in the [Julia](https://julialang.org) language.

## Installation
First, install Julia on your system following the official [Julia download instructions](https://julialang.org/downloads/) for your system. Upon successful installation of Julia, execute the following steps to install and use ZircSat:

1. Open a terminal and clone this repository to any destination on your machine:
    ```
    git clone https://github.com/lcandiot/ZircSat.jl.git
    ```

2. In your local repositories top level, open a Julia REPL by typing
    ```
    julia --project
    ```

 3. Test if your installation is working. Start the package mode by typing `]` followed by `test`. You can exit the package mode by hitting backspace/delete on your keyboard.
 4. In the [examples](https://github.com/lcandiot/ZircSat.jl/tree/main/examples) directory you can find scripts that illustrate how to use this package. See how to calculate the liquidus temperature by copy-pasting 
    ```
    include("./examples/ex_calculate_liquidus_temperature.jl")
    ```
    in your REPL.

## Input / Output
The `calculate_liquidus_temperature()` function expects a data table stored as `.csv` file. This table should contain the test pressure in kbar as well as the major oxides `SiO2-Al2O3-CaO-MgO-FeO-Fe2O3-TiO2-K2O-Na2O-Cr2O3-H2O` which are required by MAGEMins igneous database. For now only this database has been tested.

Upon calculation of the liquidus temperature, a new `.csv` file will be written to the same location as the input file. The new file contains the original data plus an additional column in which the calculated liquidus temperature is stored. Although the precision of this prediction is $\approx$ 1 °C for the test cases, the accuracy of the predicted liquidus temperature is dependent on the accuracy of MAGEMin predicting the stable mineral phases correctly.
