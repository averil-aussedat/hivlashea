# HIVLASHEA: High order methods for Vlasov-Poisson model for sheath

## Installation of dependencies! ##

Install c compiler (compile the code):

    $ sudo apt-get install gcc

Install parallel compiler + hdf5 library (compile the code for parallel use + library for data visualisation):

    $ sudo apt-get install libopenmpi-dev openmpi-bin libhdf5-openmpi-dev

Install source code manager (to share the code among us):

    $ sudo apt-get install git

Install visualize it (to visualise the hdf5 data):

    https://visit-dav.github.io/visit-website/releases-as-tables/#latest

## Run the project ##

To compile the project:

    $ ./compile_hivlashea.sh

To update the yaml file with all the simulation parameters (you can copy it beforehand for another test case):

    $ gedit ./yaml/cemracs2022_2sp.yaml

To run the project: after compiling, you should have a hivlashea.out file in the build/ folder. You can e.g. go in the build folder and run the project:

    $ ./hivlashea.out ../yaml/cemracs2022_2sp.yaml

To make a parallel run of the project (on your local computer), you can e.g. go in the build folder and run the project:

    $ mpirun -np 1 ./hivlashea.out ../yaml/cemracs2022_2sp.yaml
    $ mpirun -np 2 ./hivlashea.out ../yaml/cemracs2022_2sp.yaml

To make a parallel run of the project (on a distant server), this time it is more complicated. It depends on the supercomputing center... Usually you have to:

a. Write another compilation script (just replace the library paths by the correct paths),
   in order to compile the project on the distant computer.

b. Write a "job" file that will be executed (syntax depending on the distant computer).

## How to use git (basics) ##

To get the code, you have two ways to do it, it might depend on your connection:

git clone https://github.com/averil-prost/hivlashea.git
git clone git@github.com:averil-prost/hivlashea.git

If you need to make any modification from your local copy to the distant project, first look at which files you modified:

    $ git status

If you are happy with all modification, you can add all of them with:

    $ git add .
    
If you need to add only some of them, type instead:

    $ git add FILE1 FILE2...

Then you need to commit your changes (with a small message explaining the changes) and push your commit on the distant server:
    
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

