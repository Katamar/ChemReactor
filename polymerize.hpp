void print_double_matrix(std::vector<std::vector<double>> &matrix){
    for ( const std::vector<double> &v : matrix ){
       for ( double x : v ) std::cout << x << ' ';
           std::cout << std::endl;
    }
}
void print_int_matrix(std::vector<std::vector<int>> &matrix){
    for ( const std::vector<int> &v : matrix ){
       for ( int x : v ) std::cout << x << ' ';
           std::cout << std::endl;
    }
}

void append_double_to_matrix(std::vector<std::vector<double>> &matrix, double number){

    // define row
    std::vector<double> vector1(matrix.size(), number);
    // add row
    matrix.push_back(vector1);
    // add column
    for (int k = 0 ; k < matrix.size() ; ++k){
        matrix[k].push_back(number);
    }

}

void append_zeros_double_to_matrix(std::vector<std::vector<double>> &matrix){

    // define row
    std::vector<double> vector1(matrix.size(), 0.0);
    // add row
    matrix.push_back(vector1);
    // add column
    for (int k = 0 ; k < matrix.size() ; ++k){
        matrix[k].push_back(0.0);
    }

}

void append_zeros_int_to_matrix(std::vector<std::vector<int>> &matrix){

    // define row
    std::vector<int> vector1(matrix.size(), 0);
    // add row
    matrix.push_back(vector1);
    // add column
    for (int k = 0 ; k < matrix.size() ; ++k){
        matrix[k].push_back(0);
    }

}

void append_ones_double_to_matrix(std::vector<std::vector<double>> &matrix){

    // define row
    std::vector<double> vector1(matrix.size(), 1.0);
    // add row
    matrix.push_back(vector1);
    // add column
    for (int k = 0 ; k < matrix.size() ; ++k){
        matrix[k].push_back(1.0);
    }

}

void append_ones_int_to_matrix(std::vector<std::vector<int>> &matrix){

    // define row
    std::vector<int> vector1(matrix.size(), 1);
    // add row
    matrix.push_back(vector1);
    // add column
    for (int k = 0 ; k < matrix.size() ; ++k){
        matrix[k].push_back(1);
    }

}


void add_element(double box_size, double kbT, int &bead_number, double m, std::vector<double> &x_coor, std::vector<double> &y_coor, std::vector<double> &z_coor, std::vector<double> &new_x_coor, std::vector<double> &new_y_coor, std::vector<double> &new_z_coor, std::vector<double> &x_coor_pbc_whole, std::vector<double> &y_coor_pbc_whole, std::vector<double> &z_coor_pbc_whole, std::vector<double> &vx, std::vector<double> &vy, std::vector<double> &vz , std::vector<double> &new_vx, std::vector<double> &new_vy, std::vector<double> &new_vz, std::vector<double> &harter_x, std::vector<double> &harter_y, std::vector<double> &harter_z, std::vector<int> &bond_num, std::vector<double> &ranr_x, std::vector<double> &ranr_y, std::vector<double> &ranr_z){
    x_coor.push_back(((double)rand()/(double)RAND_MAX)*box_size);
    y_coor.push_back(((double)rand()/(double)RAND_MAX)*box_size);
    z_coor.push_back(((double)rand()/(double)RAND_MAX)*box_size);
    new_x_coor.push_back(0.0);
    new_y_coor.push_back(0.0);
    new_z_coor.push_back(0.0);
    x_coor_pbc_whole.push_back(0.0);
    y_coor_pbc_whole.push_back(0.0);
    z_coor_pbc_whole.push_back(0.0);
    vx.push_back(generateGaussianNoise(0.0, 1.0, 0.5, 0.5, 6.283185307179586)*sqrt(kbT*0.6022/m));
    vy.push_back(generateGaussianNoise(0.0, 1.0, 0.5, 0.5, 6.283185307179586)*sqrt(kbT*0.6022/m));
    vz.push_back(generateGaussianNoise(0.0, 1.0, 0.5, 0.5, 6.283185307179586)*sqrt(kbT*0.6022/m));
    new_vx.push_back(0.0);
    new_vy.push_back(0.0);
    new_vz.push_back(0.0);
    harter_x.push_back(0.0);
    harter_y.push_back(0.0);
    harter_z.push_back(0.0);
    ranr_x.push_back(0.0);
    ranr_y.push_back(0.0);
    ranr_z.push_back(0.0);
    bond_num.push_back(0);
}

void add_particle(int particle_index, int count_added_particles, double box_size, double kbT, int &bead_number, double m, std::vector<double> &x_coor, std::vector<double> &y_coor, std::vector<double> &z_coor, std::vector<double> &new_x_coor, std::vector<double> &new_y_coor, std::vector<double> &new_z_coor, std::vector<double> &x_coor_pbc_whole, std::vector<double> &y_coor_pbc_whole, std::vector<double> &z_coor_pbc_whole, std::vector<double> &vx, std::vector<double> &vy, std::vector<double> &vz , std::vector<double> &new_vx, std::vector<double> &new_vy, std::vector<double> &new_vz, std::vector<double> &harter_x, std::vector<double> &harter_y, std::vector<double> &harter_z, std::vector<int> &bond_num, std::vector<int> &identities, std::vector<std::string> &chain_atom_ID_list, std::vector<std::vector<int>> &matrix_cpy, std::vector<std::vector<double>> &matrix_harter_x, std::vector<std::vector<double>> &matrix_harter_y, std::vector<std::vector<double>> &matrix_harter_z, std::vector<std::vector<double>> &matrix_bond_r, std::vector<std::vector<int>> &mask_end, std::vector<std::vector<int>> &matrix, std::vector<std::vector<double>>  &probability_matrix_plus, std::vector<std::vector<double>>  &probability_matrix_minus, std::vector<double> proportions, std::vector<double> p_plus_probs, std::vector<double> p_minus_probs, double ranr, std::vector<double> &ranr_x, std::vector<double> &ranr_y, std::vector<double> &ranr_z){
    identities.push_back(identities[particle_index]);
    chain_atom_ID_list.push_back(std::to_string(bead_number+count_added_particles-1)); 
    append_zeros_int_to_matrix(matrix);
    append_zeros_int_to_matrix(matrix_cpy);
    append_zeros_double_to_matrix(matrix_harter_x);
    append_zeros_double_to_matrix(matrix_harter_y);
    append_zeros_double_to_matrix(matrix_harter_z);
    append_zeros_double_to_matrix(matrix_bond_r);
    append_ones_int_to_matrix(mask_end);
    append_zeros_double_to_matrix(probability_matrix_plus);
    append_zeros_double_to_matrix(probability_matrix_minus);
    add_element(box_size, kbT, bead_number, m, x_coor, y_coor, z_coor, new_x_coor, new_y_coor, new_z_coor, x_coor_pbc_whole, y_coor_pbc_whole, z_coor_pbc_whole, vx, vy, vz, new_vx, new_vy, new_vz, harter_x, harter_y, harter_z, bond_num, ranr_x, ranr_y, ranr_z);
}


int polymerize_simple_with_multiple_probabilities(std::vector<std::vector<int>> &matrix, std::vector<int> &bond_num, int &bead_number, double ranr, double dist, double binding_cutoff, std::vector<std::vector<double>>  &probability_matrix_plus, std::vector<std::vector<double>>  &probability_matrix_minus, std::vector<double> &x_coor, std::vector<double> &y_coor, std::vector<double> &z_coor, double box_size, std::vector<std::string> &to_delete, std::vector<int> &sequence, std::vector<std::string> &chain_atom_ID_list, std::vector<int> &identities, std::vector<std::vector<int>> &matrix_cpy, std::vector<std::vector<double>> &matrix_harter_x, std::vector<std::vector<double>> &matrix_harter_y, std::vector<std::vector<double>> &matrix_harter_z, std::vector<std::vector<double>> &matrix_bond_r, std::vector<std::vector<int>> &mask_end, std::vector<double> &vx, std::vector<double> &vy, std::vector<double> &vz, std::vector<double> &new_vx, std::vector<double> &new_vy, std::vector<double> &new_vz, std::vector<double> &x_coor_pbc_whole, std::vector<double> &y_coor_pbc_whole, std::vector<double> &z_coor_pbc_whole, std::vector<double> &new_x_coor, std::vector<double> &new_y_coor, std::vector<double> &new_z_coor, double kbT, std::vector<double> &harter_x, std::vector<double> &harter_y, std::vector<double> &harter_z, double m, std::vector<double> proportions, std::vector<double> p_plus_probs, std::vector<double> p_minus_probs, int &count_added_particles, std::vector<double> &ranr_x, std::vector<double> &ranr_y, std::vector<double> &ranr_z){
    double p_plus_max = 0.8;
    double p_plus_min = 0.5;
    count_added_particles = 0;
    for(int s=0; s<bead_number; ++s){
        for(int t=s+1; t<bead_number; ++t){
            dist = sqrt(pow(((x_coor[s]-x_coor[t]) - round((x_coor[s]-x_coor[t])/box_size)*box_size), 2)+pow(((y_coor[s]-y_coor[t]) - round((y_coor[s]-y_coor[t])/box_size)*box_size), 2)+pow(((z_coor[s]-z_coor[t]) - round((z_coor[s]-z_coor[t])/box_size)*box_size), 2));
            if (dist <= binding_cutoff and (bond_num[s]==0 and bond_num[t]==0) and (matrix[s][t]==0 and matrix[t][s]==0) ){
		ranr = ((double)rand()/(double)RAND_MAX);
		if(ranr<=p_plus_max and identities[s] == identities[t]){
                    count_added_particles += 1;
                    add_particle(s, count_added_particles, box_size, kbT, bead_number, m, x_coor, y_coor, z_coor, new_x_coor, new_y_coor, new_z_coor, x_coor_pbc_whole, y_coor_pbc_whole, z_coor_pbc_whole, vx, vy, vz, new_vx, new_vy, new_vz, harter_x, harter_y, harter_z, bond_num, identities, chain_atom_ID_list, matrix_cpy, matrix_harter_x, matrix_harter_y, matrix_harter_z, matrix_bond_r, mask_end, matrix, probability_matrix_plus, probability_matrix_minus, proportions, p_plus_probs, p_minus_probs, ranr, ranr_x, ranr_y, ranr_z);
                    count_added_particles += 1;
                    add_particle(t, count_added_particles, box_size, kbT, bead_number, m, x_coor, y_coor, z_coor, new_x_coor, new_y_coor, new_z_coor, x_coor_pbc_whole, y_coor_pbc_whole, z_coor_pbc_whole, vx, vy, vz, new_vx, new_vy, new_vz, harter_x, harter_y, harter_z, bond_num, identities, chain_atom_ID_list, matrix_cpy, matrix_harter_x, matrix_harter_y, matrix_harter_z, matrix_bond_r, mask_end, matrix, probability_matrix_plus, probability_matrix_minus, proportions, p_plus_probs, p_minus_probs, ranr, ranr_x, ranr_y, ranr_z);
                    matrix[s][t] = 1;
                    matrix[t][s] = 1;
                    bond_num[s] = 1;
                    bond_num[t] = 1;
		    to_delete.push_back(std::to_string(s));
		    to_delete.push_back(std::to_string(t));
		    sequence.clear();
                    sequence.push_back(s);
                    sequence.push_back(t);
		    std::string s1;
                    for(int w=0; w<sequence.size(); ++w){
                        int len = sequence.size();
                        if (len<1) {break;}
                        s1.append(std::to_string(sequence[w]));
                        s1.append(1u,' ');
                    }
                    chain_atom_ID_list.push_back(s1);
		    forward_overlap(chain_atom_ID_list, to_delete);
		    }
		if(ranr<=p_plus_min and identities[s] != identities[t]){
                    count_added_particles += 1;
		    // add for s particle
                    add_particle(s, count_added_particles, box_size, kbT, bead_number, m, x_coor, y_coor, z_coor, new_x_coor, new_y_coor, new_z_coor, x_coor_pbc_whole, y_coor_pbc_whole, z_coor_pbc_whole, vx, vy, vz, new_vx, new_vy, new_vz, harter_x, harter_y, harter_z, bond_num, identities, chain_atom_ID_list, matrix_cpy, matrix_harter_x, matrix_harter_y, matrix_harter_z, matrix_bond_r, mask_end, matrix, probability_matrix_plus, probability_matrix_minus, proportions, p_plus_probs, p_minus_probs, ranr, ranr_x, ranr_y, ranr_z);
                    count_added_particles += 1;
		    // add for t particle
                    add_particle(t, count_added_particles, box_size, kbT, bead_number, m, x_coor, y_coor, z_coor, new_x_coor, new_y_coor, new_z_coor, x_coor_pbc_whole, y_coor_pbc_whole, z_coor_pbc_whole, vx, vy, vz, new_vx, new_vy, new_vz, harter_x, harter_y, harter_z, bond_num, identities, chain_atom_ID_list, matrix_cpy, matrix_harter_x, matrix_harter_y, matrix_harter_z, matrix_bond_r, mask_end, matrix, probability_matrix_plus, probability_matrix_minus, proportions, p_plus_probs, p_minus_probs, ranr, ranr_x, ranr_y, ranr_z);

                    matrix[s][t] = 1;
                    matrix[t][s] = 1;
                    bond_num[s] = 1;
                    bond_num[t] = 1;
		    to_delete.push_back(std::to_string(s));
		    to_delete.push_back(std::to_string(t));
		    sequence.clear();
                    sequence.push_back(s);
                    sequence.push_back(t);
		    std::string s1;
                    for(int w=0; w<sequence.size(); ++w){
                        int len = sequence.size();
                        if (len<1) {break;}
                        s1.append(std::to_string(sequence[w]));
                        s1.append(1u,' ');
                    }
                    chain_atom_ID_list.push_back(s1);
		    forward_overlap(chain_atom_ID_list, to_delete);
		    }
	    }
        } // pt s
    } // pt t      
    return count_added_particles;
}	



void no_adding_particles_polymerize_simple_with_multiple_probabilities(std::vector<std::vector<int>> &matrix, std::vector<int> &bond_num, int bead_number, double ranr, double dist, double binding_cutoff, std::vector<std::vector<double>>  probability_matrix_plus, std::vector<double> &x_coor, std::vector<double> &y_coor, std::vector<double> &z_coor, double box_size, std::vector<std::string> &to_delete, std::vector<int> &sequence, std::vector<std::string> &chain_atom_ID_list, std::vector<int> &identities){
    double p_plus_max = 0.8;
    double p_plus_min = 0.3;
    for(int s=0; s<bead_number; ++s){
        for(int t=s+1; t<bead_number; ++t){
            dist = sqrt(pow(((x_coor[s]-x_coor[t]) - round((x_coor[s]-x_coor[t])/box_size)*box_size), 2)+pow(((y_coor[s]-y_coor[t]) - round((y_coor[s]-y_coor[t])/box_size)*box_size), 2)+pow(((z_coor[s]-z_coor[t]) - round((z_coor[s]-z_coor[t])/box_size)*box_size), 2));
            if (dist <= binding_cutoff and (bond_num[s]==0 and bond_num[t]==0) and (matrix[s][t]==0 and matrix[t][s]==0) ){
		ranr = ((double)rand()/(double)RAND_MAX);
		if(ranr<=p_plus_max and identities[s] == identities[t]){
                    matrix[s][t] = 1;
                    matrix[t][s] = 1;
                    bond_num[s] = 1;
                    bond_num[t] = 1;
		    to_delete.push_back(std::to_string(s));
		    to_delete.push_back(std::to_string(t));
		    sequence.clear();
                    sequence.push_back(s);
                    sequence.push_back(t);
		    std::vector<std::string> v;
		    std::string s1;
                    for (int w=0; w<sequence.size(); ++w){
                        int len = sequence.size();
                        if (len<1) {break;}
                        s1.append(std::to_string(sequence[w]));
                        s1.append(1u,' ');
                    }
                    chain_atom_ID_list.push_back(s1);
		    forward_overlap(chain_atom_ID_list, to_delete);
	        }
		if(ranr<=p_plus_min and identities[s] != identities[t]){
                    matrix[s][t] = 1;
                    matrix[t][s] = 1;
                    bond_num[s] = 1;
                    bond_num[t] = 1;
		    to_delete.push_back(std::to_string(s));
		    to_delete.push_back(std::to_string(t));
		    sequence.clear();
                    sequence.push_back(s);
                    sequence.push_back(t);
		    std::vector<std::string> v;
		    std::string s1;
                    for (int w=0; w<sequence.size(); ++w){
                        int len = sequence.size();
                        if (len<1) {break;}
                        s1.append(std::to_string(sequence[w]));
                        s1.append(1u,' ');
                    }
                    chain_atom_ID_list.push_back(s1);
		    forward_overlap(chain_atom_ID_list, to_delete);
	        }
	    }
        } // pt s
    } // pt t 	
}



void polymerize_simple(std::vector<std::vector<int>> &matrix, std::vector<int> &bond_num, int bead_number, double ranr, double dist, double binding_cutoff, double p_plus, int &dimer_counter, std::vector<double> &x_coor, std::vector<double> &y_coor, std::vector<double> &z_coor, double box_size, std::vector<std::string> &to_delete, std::vector<int> &sequence, std::vector<std::string> &chain_atom_ID_list){
    for(int s=0; s<bead_number; ++s){
        for(int t=s+1; t<bead_number; ++t){
            dist = sqrt(pow(((x_coor[s]-x_coor[t]) - round((x_coor[s]-x_coor[t])/box_size)*box_size), 2)+pow(((y_coor[s]-y_coor[t]) - round((y_coor[s]-y_coor[t])/box_size)*box_size), 2)+pow(((z_coor[s]-z_coor[t]) - round((z_coor[s]-z_coor[t])/box_size)*box_size), 2));
            if (dist <= binding_cutoff and (bond_num[s]==0 and bond_num[t]==0) and (matrix[s][t]==0 and matrix[t][s]==0) ){
		ranr = ((double)rand()/(double)RAND_MAX);
		if(ranr<=p_plus){
		    dimer_counter += 1;
                    matrix[s][t] = 1;
                    matrix[t][s] = 1;
                    bond_num[s] = 1;
                    bond_num[t] = 1;
		    to_delete.push_back(std::to_string(s));
		    to_delete.push_back(std::to_string(t));
		    sequence.clear();
                    sequence.push_back(s);
                    sequence.push_back(t);
		    std::vector<std::string> v;
		    std::string s1;
                    for (int w=0; w<sequence.size(); ++w){
                        int len = sequence.size();
                        if (len<1) {break;}
                        s1.append(std::to_string(sequence[w]));
                        s1.append(1u,' ');
                    }
                    chain_atom_ID_list.push_back(s1);
		    forward_overlap(chain_atom_ID_list, to_delete);
		    }
	    }
        } // pt s
    } // pt t 	
}	


	
void polymerize(std::vector<int> &sequence, std::vector<int> &collect_chain_lengths, std::vector<int> &polymer_ends, std::vector<std::vector<int>> &matrix, std::vector<int> &bond_num, int bead_number, int sum_row, int chain_number, int next_index, bool in_vector, std::ofstream &myfile6, int k, double ranr, double dist, double box_size, std::vector<double> &x_coor, std::vector<double> &y_coor, std::vector<double> &z_coor, double binding_cutoff, double p_plus, int step_growth, int chain_growth, std::vector<std::vector<int>> &mask_end, int begin_index, std::vector<std::string> &chain_atom_ID_list, std::vector<std::string> &to_delete){
    for(int s=0; s<bead_number; ++s){
        for(int t=s+1; t<bead_number; ++t){
            ranr = ((double)rand()/(double)RAND_MAX);
            dist = sqrt(pow(((x_coor[s]-x_coor[t]) - round((x_coor[s]-x_coor[t])/box_size)*box_size), 2)+pow(((y_coor[s]-y_coor[t]) - round((y_coor[s]-y_coor[t])/box_size)*box_size), 2)+pow(((z_coor[s]-z_coor[t]) - round((z_coor[s]-z_coor[t])/box_size)*box_size), 2));
            // step growth
            if (dist <= binding_cutoff and (bond_num[s]<2 and bond_num[t]<2) and (ranr<=p_plus) and step_growth == 1 and (matrix[s][t]==0 and matrix[t][s]==0) and (mask_end[s][t]==1 and mask_end[t][s]==1)){
                matrix[s][t] = 1;
                matrix[t][s] = 1;
        	bond_num[s] += 1;
        	bond_num[t] += 1;
		// put two oligomers together, they will be bonded in R-s-t-R
		if ((bond_num[s]==2 and bond_num[t]==2)){
                    for(int d=0; d<bead_number; ++d){
                        mask_end[s][d] = 1;
                        mask_end[d][s] = 1;
                        mask_end[t][d] = 1;
                        mask_end[d][t] = 1;
	            }
            	    begin_index = s;
                    delete_end(sequence, matrix, mask_end, bond_num, begin_index, bead_number, in_vector, next_index);
            	    begin_index = t;
                    delete_end(sequence, matrix, mask_end, bond_num, begin_index, bead_number, in_vector, next_index);
		    begin_index = s;
		    find_middle_end(sequence, matrix, mask_end, bond_num, begin_index, bead_number, in_vector, next_index);
		}
                // add particle s on tip of the chain
		else if (bond_num[s]==1 and bond_num[t]==2){
		    to_delete.push_back(std::to_string(s));
                    for(int d=0; d<bead_number; ++d){
                        mask_end[s][d] = 1;
                        mask_end[d][s] = 1;
                        mask_end[t][d] = 1;
                        mask_end[d][t] = 1;
	            }
		    begin_index = s;
                    find_end(sequence, matrix, mask_end, bond_num, begin_index, bead_number, in_vector, next_index);  
	        }
	        // add particle t on tip of chain 	
		else if (bond_num[t]==1 and bond_num[s]==2){
		    to_delete.push_back(std::to_string(t));
                    for(int d=0; d<bead_number; ++d){
                        mask_end[s][d] = 1;
                        mask_end[d][s] = 1;
                        mask_end[t][d] = 1;
                        mask_end[d][t] = 1;
	            }
		    begin_index = t;
                    find_end(sequence, matrix, mask_end, bond_num, begin_index, bead_number, in_vector, next_index);  
	        }
		else if (bond_num[t]==1 and bond_num[s]==1){
		    to_delete.push_back(std::to_string(s));
		    to_delete.push_back(std::to_string(t));
		    sequence.clear();
                    sequence.push_back(s);
                    sequence.push_back(t);
	        }
	        /////////////////////////////////////
                // VISUALIZATION OF POLYMER GROWTH //
                /////////////////////////////////////
		std::vector<std::string> v;
		std::string s1;
                for (int w=0; w<sequence.size(); ++w){
                    int len = sequence.size();
                    if (len<1) {break;}
                    s1.append(std::to_string(sequence[w]));
                    s1.append(1u,' ');
                }
                chain_atom_ID_list.push_back(s1);
		forward_overlap(chain_atom_ID_list, to_delete);
		mask_end[sequence[0]][sequence[sequence.size()-1]]==1;
		mask_end[sequence[sequence.size()-1]][sequence[0]]==1;
            }

            //////////////////
            // chain growth //
	    //////////////////
            if (dist <= binding_cutoff and ((bond_num[s]==0 and bond_num[t]==1) or (bond_num[t]==0 and bond_num[s]==1) or (bond_num[s]==0 and bond_num[t]==0)) and (ranr<=p_plus) and chain_growth == 1){
                    matrix[s][t] = 1;
                    matrix[t][s] = 1;
        	    bond_num[s] += 1;
        	    bond_num[t] += 1;
		    if (bond_num[s]==1 and bond_num[t]==2){
		        begin_index = s;
                        find_end(sequence, matrix, mask_end, bond_num, begin_index, bead_number, in_vector, next_index);  
	            }	
		    else if (bond_num[t]==1 and bond_num[s]==2){
		        begin_index = t;
                        find_end(sequence, matrix, mask_end, bond_num, begin_index, bead_number, in_vector, next_index);  
	            }	
		    else if (bond_num[t]==1 and bond_num[s]==1){
		        to_delete.push_back(std::to_string(s));
		        to_delete.push_back(std::to_string(t));
		        sequence.clear();
                        sequence.push_back(s);
                        sequence.push_back(t);
	            }	
	            /////////////////////////////////////
                    // VISUALIZATION OF POLYMER GROWTH //
                    /////////////////////////////////////
		    std::vector<std::string> v;
		    std::string s1;
                    for (int w=0; w<sequence.size(); ++w){
                        int len = sequence.size();
                        if (len<1) {break;}
                        s1.append(std::to_string(sequence[w]));
                        s1.append(1u,' ');
                    }
                    chain_atom_ID_list.push_back(s1);
		    forward_overlap(chain_atom_ID_list, to_delete);
		    mask_end[sequence[0]][sequence[sequence.size()-1]]==1;
		    mask_end[sequence[sequence.size()-1]][sequence[0]]==1;
            } // end of chain growth subroutine 
        }
    }
}	

