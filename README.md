# ChemReactor
Molecular simulation software that captures the spatio-temporal and chemical dynamics of the system. Can simulate different polymerization mechanisms (dimerization, bimolecular reaction, step-growth and chain-growth polymerization...). A combination of stochastic kinetics and stochastic dynamics (Langevin dynamics). 

System requirements
-------------------

In order to compile, a C++ compiler must be installed: gcc/g++ on Linux or Mac OS.

Requirements:
- gcc/g++ version 4.8 or higher
- OS: Ubuntu 16.04 Xenial Xerus LTS or higher
      Mac OS 10.12 Sierra or higher

The software has been extensively tested on Ubuntu 18.04 and 20.04.



Installation guide
------------------

Software installation (compilation) proceeds within 30 s, depending on computer configuration. Software run time depends on the parameters such, most notably number of time steps. 

1. Download the zipped code on github and unpack it
2. Using command line, enter the unzipped directory; ChemReactor-main, e.g.

    cd ~/Downloads/ChemReactor-main

3A. In the command line, type the following to compile: 

    g++ ChemReactor.cpp -o ChemReactor
    
    In the command line, type the following to run the code:

    ./ChemReactor

                --  or  --
3B. In the command line, type the following command:
     
    source computer_run.sh 



Demo
----

The code is a physical simulation and it requires the input of physical parameters, which are located in the ChemReactor.cpp. If needed, these can be changed prior to the simulation run.

1. Input files: 
   No input files are needed. The initial state of each simulation is a well-mixed solution of monomers; achieved by initializing randomly the positions of monomers.

2. The code outputs the following data files:

   2.A. identity.dat : Chemical identity or monomer types that are added initially to the simulation setup are saved in this file. It's a single column file, with length depending on how many particles were added initially. The column consists of monomer labels that are integers, e.g. 1, 2, 3 for a three monomer type simulation.

   2.B number_of_deleted_particles.dat : records the number of deleted moieties at each time step, if a continuous stirred tank reactor option is used.

   2.C ID_dimer_writer.dat : Records the monomer types of monomers and polymers in the simulation at each time step. For each time step, a list of monomers and polymers is written. Monomers are marked with their integer label, while for polymers, several integers are inputted corresponding to the labels of their monomers. e.g. a trimer 1-2-1, where 1 and 2 are monomer types, and - represents chemical bonds, will be marked as "1 2 1" in this file.

   2.D chain_list.dat : Equivalent to ID_dimer_writer.dat, except the particle simulation index is recorded for monomers and polymers. In each time step and for a total number of particles N, particle simulation index are element of [1, ....N].

   2.E chain_number.dat : Records number of chains with length >= 2.

   2.F chain_lengths.dat Records number of individual chains throughout the simulation.

   2.G bound_ID.dat Each row represents a single time step of length N, where N is the total number of particles at that time step. The file records for each particle whether it is bound to at least one other monomer (1) or it is not bound to any monomers (0).

   2.H bond_number.dat Equivalent to bound_ID.dat except the number of bonds are recorded. Possible entries are: 0: no bonds, 1: connected to one other monomer, 2: connected to two other monomers.

   2.I traj.dcd : binary file containing the trajectory of the chain, i.e. the x,y,z coordinates for each outputted snapshot (controlled by parameter "factor" in the epigenetic.cpp file. factor=10 in this version). *.dcd is a standard trajectory file in molecular simulations, read by many/most visualisation softwares, e.g. VMD (Visual Molecular Dynamics).


Expected run time for a demo
----------------------------

Data can be generated within 1-2 min, however, in order to achieve good sampling and production-quality data, the run time is several hours or days. The appropriate run time will depend on simulation parameters.


Instructions for use of data
----------------------------

Output data is stored in *.dat or *dcd files that can be used for analysis purposes.

