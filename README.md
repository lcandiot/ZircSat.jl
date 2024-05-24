# ZircSat

[![Build Status](https://github.com/lcandiot/ZircSat.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/lcandiot/ZircSat.jl/actions/workflows/CI.yml?query=branch%3Amain)

Calculate the liquidus temperature and estimate the Zircon saturation temperature based on thermodynamic models for silicate magmas.

🚧 Note: this is a preliminary version. 🚧

This package is built on [MAGEMin](https://github.com/ComputationalThermodynamics/MAGEMin), a Gibbs energy minimization tool. 

Calculation of the liquidus temperature, i.e. the temperature at which the solid fraction is $\approx$ 0, is performed using Newton-Raphson iterations and the algorithm is built entirely in the [Julia](https://julialang.org) language.

## Installation

First, install Julia following the official [Julia download instructions](https://julialang.org/downloads/) for your system. Upon successful installation of Julia, execute the following steps to install and use ZircSat:

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
   You should now see the following appear in the terminal
   <img width="617" alt="Bildschirmfoto 2024-05-24 um 14 25 56" src="https://github.com/lcandiot/ZircSat.jl/assets/50524459/9ccae926-c4aa-426e-9ce8-47b6e4f02a01">

   The green `julia>` string in the lower left corner of the terminal indicates that you are now in the standard Julia mode in which you can perform calculations and execute Julia code.

6. To test if the installation works, activate the package mode and test the ZircSat package

   ```
   ]
   test
   ```
   <img width="1023" alt="Bildschirmfoto 2024-05-24 um 14 30 31" src="https://github.com/lcandiot/ZircSat.jl/assets/50524459/d301e266-89d1-475b-b6dc-036e33a57e82">
   
   Note that the string in the lower left corner has changed to `(ZircSat) pkg>` indicating that your REPL is now in Julia's package mode. You should now see that both tests have passed succesfully as indicated in the image above.

7. You can now exit the package environment by hitting `BACKSPACE` on your keyboard. You should now again see the green `julia>` string in the lower left corner of the terminal.

8. Open `/USER/my_Tliq_Tsat_calculator.jl` with a code editor of your choice. Alternatively, you can also open it in the active REPL by typing `;` followed by 
   ```
   vi USER/my_Tliq_Tsat_calculator.jl
   ```
   Now, modifiy the data path line to
   ```
   data_path = "./USER/my_data.csv"
   ```
   and close the editor hitting `ESC`, `:wq`, and `ENTER`. Hit `BACKSPACE` again to go back to the standard `julia>` mode. You should now again see the green `julia>` string in the lower left corner of the terminal.

9. Run `/USER/my_Tliq_Tsat_calculator.jl` by entering the following command in the Julia REPL
   ```
   include("USER/my_Tliq_Tsat_calculator.jl")
   ```

10. Congratulations! You can now modify these two files according to your needs and add your own data to `/USER/my_data.csv`


## Input / Output

The `calculate_liquidus_temperature()` function expects a data table stored as `.csv` file. It is best practice to provide a comma delimited file. This table should contain the test pressure in kbar as well as the major oxides `SiO2-Al2O3-CaO-MgO-FeO-TiO2-K2O-Na2O-Cr2O3-H2O` which are required by MAGEMins igneous database. Note that `FeO` here is total iron. If only `Fe2O3` is availbale this value can be converted as `FeO  = Fe2O3 / 1.1111`, if both are given the conversion is `FeO  = FeO + ( Fe2O3 / 1.1111)`. For now only the igneous database has been tested. While the ordering of the columns does not matter, the package is case sensitive. It is therefore recommended to use the [ZircSat_test_MarxerUlmer2019.csv](https://github.com/lcandiot/ZircSat.jl/tree/main/data/ZircSat_test_MarxerUlmer2019.csv) as a template. Further, the Zirconium concentration in units ppm is required as `Zr` column in the data file.

Upon calculation of the liquidus temperature, a new `.csv` file will be written to the same location as the input file. The new file contains the original data plus three additional columns for calculated liquidus temperature, `T_liq [C]`, Zircon saturation temperature, `T_sat [C]`, and the difference between these two temperatures, `T_diff [C]`. Although the precision of this prediction is $\approx$ 1 °C for the test cases, the accuracy of the predicted liquidus temperature is dependent on the accuracy of MAGEMin predicting the stable mineral phases correctly, which is in turn constrained by the underlying igenous database.

## Working with this repository

It is recommended to generate a `/USER` directory as described above. This directory should also be the place to work from and run scripts. Feel free to use the [examples](https://github.com/lcandiot/ZircSat.jl/tree/main/examples) as starting point for your own routines.
