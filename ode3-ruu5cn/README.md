John Grant and Gary Henderson's stuff:
The terminal velocity stuff is pretty easy to run - everything is under VT, and I added everything necessary to the makefile
For cross-checking the terminal velocity stuff - we see that terminal velocity ~ sqrt(m). This makes sense, because terminal velocity occurs when (delta v) = 0 -> F = 0 -> mg = a(v_t)^2 -> v_t= sqrt(mg/a)


[![Review Assignment Due Date](https://classroom.github.com/assets/deadline-readme-button-24ddc0f5d75046c5622901739e7c5dd533143b0c8e959d652212380cedb1ea36.svg)](https://classroom.github.com/a/1FVNobkF)
# ode3

To build the ODE library and example programs, simply type `make` in this top level ode3 directory.

Description of example programs:<br>

**RKnTest**: Solves a single 1st order ODE using the single equation RK4 solver and the ODE array solver

**RKnStep**: Solves a single 1st order ODE using the single equation RK4 solver and the ODE array solver<br>
Basic version or ODE array solver applied to projectile motion with a simple model of air resistance, force of air resistance = -kv^2<br>
At each step in the calculation time and x,t positions are printed.<br>
Optional parameters [default values]<br>
* -v initial_velocity [100] m/s
* -t angle_thera [45] degrees
* -m mass_of_projectile [10] kg
* -k coefficient_of_air_resistance [0.1] kg/m


**RKnDemo**: Solves for projectile motion with a simple model of air resistance, force of air resistance = -kv^2<br>
This program includes graphical output.  Detailed output is saved in TGraph objects in RKnDemo.root.  The file **RKnPlotDemo.py** shows how to access date in the TGraphs and can be used to generate additional plots.<br>
Optional parameters [default values]<br>
* -v initial_velocity [100] m/s
* -t angle_thera [45] degrees
* -m mass_of_projectile [10] kg
* -k coefficient_of_air_resistance [0.1] kg/m

**baseball1**:  Starter template for first baseball problem

**baseball2**:  Starter template for second baseball problem

**baseball_drag.ipynb**: this notebook describes the drag force equations used in the text.


