#Propagator

----
The following program **Propagator.py** is a program that takes in parameter a **TLE.txt** file in standard 2 or 3 lines and works with multiple satellite at once.
------------------
# How to use the program
First you'll need to provide any TLE file into data/INPUT folder

Then run 

> **Propagator.py**

You will need to select the TLE file and boundaries of the study
>Begining
>Increment in sec
>Ending

You'll next need to choose the output format
PVT = coordinates in [X Y Z Xdot Ydot Zdot]
LLA = coordinates in [Latitude Longitude Altitude]

N.B.: All epochs are given in VTS standard 

A File.txt will be created for each separate Satellite into data/OUTPUT folder
> If you run multiple studies, it might overwrite preexisting file !

Files using PVT coordinates can be used in VTS


# Improvement
-Adding multiple output options (GePro, Cesium, .csv, etc)
-Overwrite security
-Graphic interface
-Innaccurate input security