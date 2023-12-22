# Rocket Simulation

## Project Description

This program simulates the trajectory of a rocket and tell us at what time we should launch the rocket in order to achieve the international space station.
It will
1. Calculate ISS coordinates ("*coord_skyfield*")
2. Find the parameters ("*simulation.c*") and plot the graphs ("*graphs.py*")
3. Determine launch time ("*simulation.c*")


## Project Structure
- "*coord_skyfield*" contains a function that returns the coordinates of ISS at the time of interest
- "*graphs.py*" plots the graphs for our parameters over time
- "*simulation.c*" contains the code to run the simulation

### Implementation details

Overview
- The code is compiled and run on C. It calls Python functions with the "system" function.
- Python is used to plot the graphs and find variables with the specific Python library Skyfield
- C is used to find the needed parameters and output of the program


## Instructions
The code is run first compile *simulation.c*
```{bash}
gcc simulation.c
```
And then run it

```{bash}
./a.exe
```



## Requirements
```{bash}
$python --version
Python 3.11.4

$gcc --version
gcc.exe (Rev7, Built by MSYS2 project) 13.1.0
```

You will also need to install the libraries found on the requirement.txt file



