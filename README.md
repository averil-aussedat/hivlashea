# HIVLASHEA: High order methods for Vlasov-Poisson model for sheath

## Installation of dependencies ##

Install c compiler (compile the code):

    $ sudo apt-get install gcc

Install parallel compiler + hdf5 library (compile the code for parallel use + library for data visualisation):

    $ sudo apt-get install libopenmpi-dev openmpi-bin libhdf5-openmpi-dev

Install source code manager (to share the code among us):

    $ sudo apt-get install git

Install visualize it (to visualise the hdf5 data):

    https://visit-dav.github.io/visit-website/releases-as-tables/#latest

## Run the project ##

0. Compile the project:

    $ ./compile_hivlashea.sh

1. Create a yaml file with all the simulation parameters.

Edit e.g. ./yaml/cemracs2022_2sp.yaml

2. Run the project.

If you have compiled the project, you should have a hivlashea.out file in the build/ folder.

You can e.g. go in the build folder and run the project:

./hivlashea.out ../yaml/cemracs2022_2sp.yaml

3. Parallel run of the project (local).

You can e.g. go in the build folder and run the project:

    $ mpirun -np 1 ./hivlashea.out
    $ mpirun -np 2 ./hivlashea.out

4. Parallel run of the project (supercomputer):

It depends on the supercomputing center... Usually you have to:

a. Write another compilation script (just replace the library paths by the correct paths),
   in order to compile the project on the distant computer.

b. Write a "job" file that will be executed (syntax depending on the distant computer).

## How to use git (basics) ##

1. Get the code.

git clone https://github.com/averil-prost/hivlashea.git

2. After a modification:

    $git status

Look at the modified files. Add the one needed with

    $git add .
    
(if you need to add all of them)

    $git add FILE1 FILE2...
    
(if you only need to add some of them)

    $git commit -m "MESSAGE"
    $git push origin master

In case of a conflict, if you have never seen one, ask for someone who knows, and learn by doing.

## How to use visit (basics) ##

bouton ouvrir sur la gauche

groupement intelligent

cliquer sur le .xmf

cliquer sur ajouter, pseudo couleur, values

cliquer sur tracer

menu contrôles, camera, 2d view, full frame on

