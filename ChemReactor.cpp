#include <iomanip>
#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <iterator>
#include <string>
#include <cmath>
#include <fstream>
#include <random>
#include <array>
#include <sstream>
#include <stdlib.h>
#include <vector>
#include <sstream>
#include "gaussianNoise.hpp"
#include "initialize.hpp"
#include "forces.hpp"
#include "integrator.hpp"
#include "pbc.hpp"
#include "reverse_function.hpp"
#include "first_digit.hpp"
#include "find_position.hpp"
#include "forward_overlap.hpp"
#include "backward_overlap.hpp"
#include "find_end.hpp"
#include "write_ends.hpp"
#include "identities.hpp"
#include "polymerize.hpp"
#include "depolymerize.hpp"
#include "CFSTR_add_remove.hpp"
#include "dcd.h"
using namespace std;
std::string line;
ifstream infile;

void get_configuration(std::vector<double> &x_coor, std::vector<double> &y_coor, std::vector<double> &z_coor, int bead_number){
infile.open("coor_column.dat");// file containing numbers in 3 columns 
int ctr = 0;
while (std::getline(infile, line))
    {
        std::stringstream ss(line);
        double a, b, c;
        ss >> a;
        ss >> b;
        ss >> c;
	x_coor.push_back(a);
	y_coor.push_back(b);
	z_coor.push_back(c);
        ctr+=1;
    }
}

// for the code to work properly, please provide quantities in units that are written in comments. 
// see interator.cpp for more info. Rescaling factors have been introduced to accept units as denoted below.

long long int seed       = 1000;
double kbT               = 4.1418; //pN*nm
double sigma_lj          = 0.332; //nm
double epsilon_lj        = 0.10*4.1418; // pN*nm
double d0                = 1.01*sigma_lj; //nm
double harmonic_constant = 10.0*epsilon_lj/(pow(d0, 2)); // [harmonic_constant]= kcal/mol/nm^2
double dt                = 0.054;   // [dt] = fs
double gamma_coeff       = 0.00925926;
double m                 = 1.0; // Daltons, mass of 10 base pairs that are approx in the 3.18 nm of the DNA
double rcut              = pow(2.0,(1.0/6.0))*sigma_lj; // cutoff distance for Lennard-Jones
long long int N_steps    = 1500;
int maximum_particle_number = 16001;
double box_size          = 10.0;
int factor               = 50;
int addParticle          = 1;
double binding_cutoff    = 1.0*sigma_lj;
int step_growth          = 0;
int chain_growth         = 1;
int sum_bond_num         = 0;
int sum_row              = 0;
int chain_number         = 0;
/////////////////////////////////////////
// parameters for time-dependent rates //
/////////////////////////////////////////
int time_dependent       = 0; //is the rate time-dependent or not, if it is, then a piece-wise function dictates the rate
int period               = 10;
int ct_td                = 0; // counter associated to the time-dependence  
/////////////////////////////////////////////////////////////////
// parameters for continuous-flow stirred tank reactor (CFSTR) //
/////////////////////////////////////////////////////////////////
double CFSTR_remove_p = 0.075; // probability to remove particles
int CFSTR_add_number_1 = 20;             // added monomers of type 1 
int CFSTR_add_number_2 = 20;             // added monomers of type 2 
int bead_number  = CFSTR_add_number_1+CFSTR_add_number_2;
///////////////////////////////
// probabilities of reaction //
///////////////////////////////
double p_plus_min        = 0.8;
double p_plus_max        = 0.5;
double p_minus_min       = 0.000;
double p_minus_max       = 0.000;
double ranr;
double dist;
double pressure;
int next_index; 
int tmp_idx;
bool in_vector; 
int begin_index;
vector<double> vx, vy, vz;
vector<double> x_coor, y_coor, z_coor; 
vector<double> new_x_coor, new_y_coor, new_z_coor;
vector<double> vdW_x, vdW_y, vdW_z;
vector<double> new_vx, new_vy, new_vz;
vector<double> ranr_x, ranr_y, ranr_z;
vector<double> harter_x, harter_y, harter_z;
vector<double> x_coor_pbc_whole, y_coor_pbc_whole, z_coor_pbc_whole;
vector<int> bond_num;
vector<int> polymer_ends = {0, 0};
vector<int> sequence;
vector<int> collect_chain_lengths;
vector<string> chain_atom_ID_list;      
vector<string> to_delete; 
vector<string> to_add; 



int ctr=0;

double u1 = 0.5, u2 = 0.5, two_pi = 2.0*3.14159265358979323846;

//void *dcd = open_dcd_write("traj.dcd", "dcd", 1, dt, bead_number);

///////////////////////////
// 2D matrix allocations //
///////////////////////////
vector<vector<int>> matrix(bead_number, vector<int>(bead_number, 0));
vector<vector<int>> matrix_cpy(bead_number, vector<int>(bead_number, 0));
vector<vector<double>> matrix_harter_x(bead_number, vector<double>(bead_number, 0.0));
vector<vector<double>> matrix_harter_y(bead_number, vector<double>(bead_number, 0.0));
vector<vector<double>> matrix_harter_z(bead_number, vector<double>(bead_number, 0.0));
vector<vector<double>> matrix_bond_r(bead_number, vector<double>(bead_number, 0.0));
vector<vector<int>> mask_end(bead_number, vector<int>(bead_number, 1));
vector<vector<double>> probability_matrix_plus(bead_number, vector<double>(bead_number, 0));
vector<vector<double>> probability_matrix_minus(bead_number, vector<double>(bead_number, 0));
vector<double> proportions;
vector<int> identities;
int count_added_particles = 0;
int count_deleted_particles = 0;
int memoryfactor = 100;


int main(){ 


//////////////////////////////////////////////
// INITIALIZE PROPORTIONS AND PROBABILITIES //
//////////////////////////////////////////////
vector<double> p_plus_probs;
vector<double> p_minus_probs;
proportions.push_back(CFSTR_add_number_1);
proportions.push_back(CFSTR_add_number_2);
p_plus_probs.push_back(0.8);
p_plus_probs.push_back(0.8);
p_minus_probs.push_back(0.005);
p_minus_probs.push_back(0.005);

ofstream myfile1;
myfile1.open ("bound_ID.dat");
ofstream myfile2;
myfile2.open ("bond_number.dat");
ofstream myfile4;
myfile4.open ("chain_number.dat");
ofstream myfile5;
myfile5.open ("chain_lengths.dat");
ofstream myfile6;
myfile6.open ("chain_list.dat");
ofstream myfile7;
myfile7.open ("identity.dat");
ofstream myfile9;
myfile9.open ("ID_dimer_writer.dat");
ofstream myfile10;
myfile10.open ("number_of_deleted_particles.dat");


srand(seed);
init_velocities(vx, vy, vz, new_vx, new_vy, new_vz, bead_number, kbT, m);
init_positions(x_coor, y_coor, z_coor, new_x_coor, new_y_coor, new_z_coor, x_coor_pbc_whole, y_coor_pbc_whole, z_coor_pbc_whole, box_size, bead_number); 
init_ranr(ranr_x, ranr_y, ranr_z, bead_number);
initialize_arrays(harter_x, harter_y, harter_z, new_x_coor, new_y_coor, new_z_coor, x_coor_pbc_whole, y_coor_pbc_whole, z_coor_pbc_whole, bond_num, bead_number, chain_atom_ID_list, vdW_x, vdW_y, vdW_z, identities);
assign_identity(identities, bead_number, ranr, proportions, probability_matrix_plus, probability_matrix_minus, p_plus_probs, p_minus_probs);


for(int mm=0; mm<bead_number; ++mm){
    myfile7 << identities[mm] << endl; 
}

for(int k=0; k<N_steps; ++k){
    to_add.clear();
    to_delete.clear();
    count_added_particles = 0;
    count_deleted_particles = 0;

    if (k%1==0){
        cout << "step: "<< k<< " bead number: " << bead_number << endl;
    }

    /////////////////////////
    // TIME-DEPENDENT RATE //
    /////////////////////////
    if (time_dependent){
        if (ct_td < period){
            p_plus_probs[0]     = p_plus_min;
            p_plus_probs[1]     = p_plus_min;
            p_minus_probs[0]    = p_minus_min;
            p_minus_probs[1]    = p_minus_min;
        }
        if (ct_td >= period){
            p_plus_probs[0]     = p_plus_max;
            p_plus_probs[1]     = p_plus_max;
            p_minus_probs[0]    = p_minus_max;
            p_minus_probs[1]    = p_minus_max;
        }
        ct_td += 1;
        if (ct_td == period*2){
            ct_td = 0;
	}
	cout << p_plus_probs[0] << endl;
        append_probabilities_to_matrix(probability_matrix_plus, probability_matrix_minus, identities, bead_number, p_plus_probs, p_minus_probs);
    } // end of time-dependent condition

    if (bead_number>maximum_particle_number){
    break;
    }

    ////////////////////
    // POLYMERIZATION //
    ////////////////////
    no_adding_particles_polymerize_simple_with_multiple_probabilities(matrix, bond_num, bead_number, ranr, dist, binding_cutoff, probability_matrix_plus, x_coor, y_coor, z_coor, box_size, to_delete, sequence, chain_atom_ID_list, identities);
    //////////////////////
    // DEPOLYMERIZATION //
    //////////////////////
    no_removing_particles_depolymerize_simple_with_multiple_probabilities(x_coor, y_coor, z_coor,  bond_num, matrix, mask_end, sequence, ranr, step_growth, chain_growth, dist, box_size, probability_matrix_minus, begin_index, bead_number, in_vector, next_index, chain_atom_ID_list); 
    if (k%factor==0){ 
    //////////////////////////////////////////
    // CONTINUOUS-FLOW STIRRED TANK REACTOR //
    //////////////////////////////////////////
    // add particles
    CFSTR_add(CFSTR_add_number_1, CFSTR_add_number_2, matrix, bond_num, bead_number, ranr, dist, binding_cutoff, probability_matrix_plus, probability_matrix_minus, x_coor, y_coor, z_coor, box_size, to_delete, sequence, chain_atom_ID_list, identities, matrix_cpy, matrix_harter_x, matrix_harter_y, matrix_harter_z, matrix_bond_r, mask_end, vx, vy, vz, new_vx, new_vy, new_vz, x_coor_pbc_whole, y_coor_pbc_whole, z_coor_pbc_whole, new_x_coor, new_y_coor, new_z_coor, kbT, harter_x, harter_y, harter_z, m, proportions, p_plus_probs, p_minus_probs, count_added_particles, ranr_x, ranr_y, ranr_z);
    bead_number += count_added_particles;
    //remove particles
    CFSTR_remove(ranr, to_delete, count_added_particles, box_size, kbT, bead_number, m, x_coor, y_coor, z_coor, new_x_coor, new_y_coor, new_z_coor, x_coor_pbc_whole, y_coor_pbc_whole, z_coor_pbc_whole, vx, vy, vz , new_vx, new_vy, new_vz, harter_x, harter_y, harter_z, bond_num, chain_atom_ID_list, matrix_cpy, matrix_harter_x, matrix_harter_y, matrix_harter_z, matrix_bond_r, mask_end, matrix, probability_matrix_plus, probability_matrix_minus, identities, ranr_x, ranr_y, ranr_z, count_deleted_particles, CFSTR_remove_p);
    bead_number -= count_deleted_particles;
    myfile10 << count_deleted_particles << endl;
    //////////////////////////
    // PROBABILITY MATRICES //
    //////////////////////////
    append_probabilities_to_matrix(probability_matrix_plus, probability_matrix_minus, identities, bead_number, p_plus_probs, p_minus_probs);
    } // condition on how often one adds/removes particles

    if(k%factor==0){
        myfile6 << "timestep " << k << endl;
        myfile9 << "timestep " << k << endl;
        for(int r=0; r<chain_atom_ID_list.size(); ++r){
	    std::string s2;
            myfile6 << chain_atom_ID_list[r] << endl;
            std::istringstream is( chain_atom_ID_list[r] );
            int n;
            while( is >> n ) {
		 s2.append(std::to_string(identities[n]));
                 s2.append(1u,' ');
            }
            myfile9 << s2 << endl;
        }
    }
    to_delete.clear();
    // integration 1
    harmonic_force(x_coor, y_coor, z_coor, harter_x, harter_y, harter_z, matrix, mask_end, matrix_harter_x, matrix_harter_y, matrix_harter_z, matrix_bond_r, bead_number, harmonic_constant, d0, box_size);
    update_ranr(ranr_x, ranr_y, ranr_z, bead_number);
    integrate_pt1(x_coor, y_coor, z_coor, new_x_coor, new_y_coor, new_z_coor, vx, vy, vz, new_vx, new_vy, new_vz, ranr_x, ranr_y, ranr_z, harter_x, harter_y, harter_z, gamma_coeff, dt, m, kbT, bead_number);
    // integration 2
    harmonic_force(x_coor, y_coor, z_coor, harter_x, harter_y, harter_z, matrix, mask_end, matrix_harter_x, matrix_harter_y, matrix_harter_z, matrix_bond_r, bead_number, harmonic_constant, d0, box_size); 
    integrate_pt2(x_coor, y_coor, z_coor, new_x_coor, new_y_coor, new_z_coor, vx, vy, vz, new_vx, new_vy, new_vz, ranr_x, ranr_y, ranr_z, harter_x, harter_y, harter_z, gamma_coeff, dt, m, kbT, bead_number);
    if(k%factor == 0){
        string s1;
        string s2;
        for(int f=0; f<bead_number; ++f){
            s2.append(to_string(bond_num[f]));
            s2.append(1u,' ');
            if (bond_num[f] == 1 or bond_num[f] == 2){
                s1.append(to_string(1));
                s1.append(1u,' ');
            }
            if (bond_num[f] == 0 ){
                s1.append(to_string(0));
                s1.append(1u,' ');
            }
        }
        myfile1 << s1 << endl;
        myfile2 << s2 << endl;
         write_ends(sequence, polymer_ends, matrix, bond_num, bead_number, sum_row, chain_number, next_index, in_vector, myfile4, myfile5, k);
    }


}
          
myfile1.close();
myfile2.close();
myfile4.close();
myfile5.close();
myfile6.close();
myfile7.close();
myfile9.close();
myfile10.close();
return 0;
}
