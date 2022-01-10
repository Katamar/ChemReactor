

void delete_rowcolumn_int_matrix(int index, std::vector<std::vector<int>> &matrix){
    int rowIndex = index;
    int columnIndex = index;
    
    matrix.erase(matrix.begin()+rowIndex);
    
    std::for_each(matrix.begin(), matrix.end(), [&](std::vector<int>& row) {
        row.erase(std::next(row.begin(), columnIndex));
    });
}

void delete_rowcolumn_double_matrix(int index, std::vector<std::vector<double>> &matrix){
    int rowIndex = index;
    int columnIndex = index;
    
    matrix.erase(matrix.begin()+rowIndex);
    
    std::for_each(matrix.begin(), matrix.end(), [&](std::vector<double>& row) {
        row.erase(std::next(row.begin(), columnIndex));
    });
}

void delete_index_int_vector(int index, std::vector<int> &vct){
    vct.erase(vct.begin()+index);
}

void delete_index_double_vector(int index, std::vector<double> &vct){
    vct.erase(vct.begin()+index);
}



void delete_element(int index, double box_size, double kbT, int &bead_number, double m, std::vector<double> &x_coor, std::vector<double> &y_coor, std::vector<double> &z_coor, std::vector<double> &new_x_coor, std::vector<double> &new_y_coor, std::vector<double> &new_z_coor, std::vector<double> &x_coor_pbc_whole, std::vector<double> &y_coor_pbc_whole, std::vector<double> &z_coor_pbc_whole, std::vector<double> &vx, std::vector<double> &vy, std::vector<double> &vz , std::vector<double> &new_vx, std::vector<double> &new_vy, std::vector<double> &new_vz, std::vector<double> &harter_x, std::vector<double> &harter_y, std::vector<double> &harter_z, std::vector<int> &bond_num, std::vector<int> &identities, std::vector<double> &ranr_x, std::vector<double> &ranr_y, std::vector<double> &ranr_z){
    delete_index_double_vector(index, x_coor);
    delete_index_double_vector(index, y_coor);
    delete_index_double_vector(index, z_coor);
    delete_index_double_vector(index, new_x_coor);
    delete_index_double_vector(index, new_y_coor);
    delete_index_double_vector(index, new_z_coor);
    delete_index_double_vector(index, x_coor_pbc_whole);
    delete_index_double_vector(index, y_coor_pbc_whole);
    delete_index_double_vector(index, z_coor_pbc_whole);
    delete_index_double_vector(index, vx);
    delete_index_double_vector(index, vy);
    delete_index_double_vector(index, vz);
    delete_index_double_vector(index, new_vx);
    delete_index_double_vector(index, new_vy);
    delete_index_double_vector(index, new_vz);
    delete_index_double_vector(index, harter_x);
    delete_index_double_vector(index, harter_y);
    delete_index_double_vector(index, harter_z);
    delete_index_double_vector(index, ranr_x);
    delete_index_double_vector(index, ranr_y);
    delete_index_double_vector(index, ranr_z);
    delete_index_int_vector(index, bond_num);
    delete_index_int_vector(index, identities);
}

void delete_particle(int index, std::vector<std::string> &to_delete, int count_added_particles, double box_size, double kbT, int &bead_number, double m, std::vector<double> &x_coor, std::vector<double> &y_coor, std::vector<double> &z_coor, std::vector<double> &new_x_coor, std::vector<double> &new_y_coor, std::vector<double> &new_z_coor, std::vector<double> &x_coor_pbc_whole, std::vector<double> &y_coor_pbc_whole, std::vector<double> &z_coor_pbc_whole, std::vector<double> &vx, std::vector<double> &vy, std::vector<double> &vz , std::vector<double> &new_vx, std::vector<double> &new_vy, std::vector<double> &new_vz, std::vector<double> &harter_x, std::vector<double> &harter_y, std::vector<double> &harter_z, std::vector<int> &bond_num, std::vector<std::string> &chain_atom_ID_list, std::vector<std::vector<int>> &matrix_cpy, std::vector<std::vector<double>> &matrix_harter_x, std::vector<std::vector<double>> &matrix_harter_y, std::vector<std::vector<double>> &matrix_harter_z, std::vector<std::vector<double>> &matrix_bond_r, std::vector<std::vector<int>> &mask_end, std::vector<std::vector<int>> &matrix, std::vector<std::vector<double>>  &probability_matrix_plus, std::vector<std::vector<double>>  &probability_matrix_minus, std::vector<int> &identities, std::vector<double> &ranr_x, std::vector<double> &ranr_y, std::vector<double> &ranr_z){
    delete_rowcolumn_int_matrix(index, matrix);
    delete_rowcolumn_int_matrix(index, matrix_cpy);
    delete_rowcolumn_int_matrix(index, mask_end);
    delete_rowcolumn_double_matrix(index, matrix_harter_x);
    delete_rowcolumn_double_matrix(index, matrix_harter_y);
    delete_rowcolumn_double_matrix(index, matrix_harter_z);
    delete_rowcolumn_double_matrix(index, matrix_bond_r);
    delete_rowcolumn_double_matrix(index, probability_matrix_plus);
    delete_rowcolumn_double_matrix(index, probability_matrix_minus);
    delete_element(index, box_size, kbT, bead_number, m, x_coor, y_coor, z_coor, new_x_coor, new_y_coor, new_z_coor, x_coor_pbc_whole, y_coor_pbc_whole, z_coor_pbc_whole, vx, vy, vz , new_vx, new_vy, new_vz, harter_x, harter_y, harter_z, bond_num, identities, ranr_x, ranr_y, ranr_z);
}


int depolymerize_simple_with_multiple_probabilities(int addParticle, double ranr, std::vector<std::string> &to_delete, int count_added_particles, double box_size, double kbT, int &bead_number, double m, std::vector<double> &x_coor, std::vector<double> &y_coor, std::vector<double> &z_coor, std::vector<double> &new_x_coor, std::vector<double> &new_y_coor, std::vector<double> &new_z_coor, std::vector<double> &x_coor_pbc_whole, std::vector<double> &y_coor_pbc_whole, std::vector<double> &z_coor_pbc_whole, std::vector<double> &vx, std::vector<double> &vy, std::vector<double> &vz , std::vector<double> &new_vx, std::vector<double> &new_vy, std::vector<double> &new_vz, std::vector<double> &harter_x, std::vector<double> &harter_y, std::vector<double> &harter_z, std::vector<int> &bond_num, std::vector<std::string> &chain_atom_ID_list, std::vector<std::vector<int>> &matrix_cpy, std::vector<std::vector<double>> &matrix_harter_x, std::vector<std::vector<double>> &matrix_harter_y, std::vector<std::vector<double>> &matrix_harter_z, std::vector<std::vector<double>> &matrix_bond_r, std::vector<std::vector<int>> &mask_end, std::vector<std::vector<int>> &matrix, std::vector<std::vector<double>>  &probability_matrix_plus, std::vector<std::vector<double>>  &probability_matrix_minus, std::vector<int> &identities, std::vector<double> &ranr_x, std::vector<double> &ranr_y, std::vector<double> &ranr_z){
    int count_deleted_particles = 0;
    std::vector<int> accumulate_delete;
    accumulate_delete.clear();
    std::vector<int> vctr_id;
    std::vector<std::string> chain_list_tmp;
    chain_list_tmp.clear();
    for(int m=0; m<bead_number; ++m){
        vctr_id.push_back(m);
    }
    chain_list_tmp.clear();
    for(int s=0; s<bead_number; ++s){
        for(int t=s+1; t<bead_number; ++t){
            ranr = ((double)rand()/(double)RAND_MAX);
            if ((ranr<=probability_matrix_minus[s][t]) and ((bond_num[s]==1 and bond_num[t]==1) or (bond_num[s]==1 and bond_num[t]==2) or (bond_num[s]==2 and bond_num[t]==1)) and (matrix[s][t]==1 and matrix[t][s]==1) ){
                count_deleted_particles += 2;
                matrix[s][t] = 0;
                matrix[t][s] = 0;
            	bond_num[s] =0;
            	bond_num[t] =0;
		// CHECK THIS BACKWARD OVERLAP!!!!
                backward_overlap_delete(s, t, chain_atom_ID_list);
		accumulate_delete.push_back(s);
		accumulate_delete.push_back(t);
            }
        }//particle s
    }//particle t
    std::sort(accumulate_delete.begin(), accumulate_delete.end(), std::greater<int>());
    for (int i=0; i<accumulate_delete.size(); ++i){
        delete_index_int_vector(accumulate_delete[i], vctr_id);
    }

    std::cout << " " << std::endl;
    std::cout << " " << std::endl;
    std::cout << " " << std::endl;
    for (int w=0; w<chain_atom_ID_list.size(); ++w){
        std::string s2;
        std::istringstream is( chain_atom_ID_list[w] );
        int n;
	s2.clear();
        while( is >> n ) {
		std::cout << " n: " << n << std::endl;
            for (int i=0; i<vctr_id.size(); ++i){
		 if(vctr_id[i] == n){
	             std::cout << "vectr_id: " << vctr_id[i] << " n: " << n << " replace: " << i << " " << " with: " << n << std::endl;
                     s2.append(std::to_string(i));
                     s2.append(1u,' ');
		 }
	    }
        }
	chain_list_tmp.push_back(s2);
    }
    std::sort(accumulate_delete.rbegin(), accumulate_delete.rend());
    for(int i=0; i<accumulate_delete.size(); ++i){
        delete_particle(accumulate_delete[i], to_delete, count_added_particles, box_size, kbT, bead_number, m, x_coor, y_coor, z_coor, new_x_coor, new_y_coor, new_z_coor, x_coor_pbc_whole, y_coor_pbc_whole, z_coor_pbc_whole, vx, vy, vz, new_vx, new_vy, new_vz, harter_x, harter_y, harter_z, bond_num, chain_atom_ID_list, matrix_cpy, matrix_harter_x, matrix_harter_y, matrix_harter_z, matrix_bond_r, mask_end, matrix, probability_matrix_plus, probability_matrix_minus, identities, ranr_x, ranr_y, ranr_z);
    }

    chain_atom_ID_list = chain_list_tmp;

    return count_deleted_particles;
}


void no_removing_particles_depolymerize_simple_with_multiple_probabilities(std::vector<double> &x_coor, std::vector<double> &y_coor, std::vector<double> &z_coor,  std::vector<int> &bond_num, std::vector<std::vector<int>> &matrix, std::vector<std::vector<int>> &mask_end, std::vector<int> &sequence, double ranr, int step_growth, int chain_growth, double dist, double box_size, std::vector<std::vector<double>>  probability_matrix_minus, int begin_index, int bead_number, bool in_vector, int next_index, std::vector<std::string> &chain_atom_ID_list){ 

        for(int s=0; s<bead_number; ++s){
            for(int t=s+1; t<bead_number; ++t){
                ranr = ((double)rand()/(double)RAND_MAX);
                dist = sqrt(pow(((x_coor[s]-x_coor[t]) - round((x_coor[s]-x_coor[t])/box_size)*box_size), 2)+pow(((y_coor[s]-y_coor[t]) - round((y_coor[s]-y_coor[t])/box_size)*box_size), 2)+pow(((z_coor[s]-z_coor[t]) - round((z_coor[s]-z_coor[t])/box_size)*box_size), 2));
                // depolymerization anywhere!
                if ((ranr<=probability_matrix_minus[s][t]) and  (bond_num[s]==1 and bond_num[t]==1) and (matrix[s][t]==1 and matrix[t][s]==1)){
		    if (bond_num[s]==1 and bond_num[t]==1){
                        backward_overlap(s, t, chain_atom_ID_list);
                        matrix[s][t] = 0;
                        matrix[t][s] = 0;
            	        bond_num[s] -=1;
            	        bond_num[t] -=1;
            	    }
		    /////////////////////////////////////
                    // VISUALIZATION OF POLYMER GROWTH //
                    /////////////////////////////////////
	            std::string s1;
                    for (int w=0; w<sequence.size(); ++w){
                        int len = sequence.size();
                        if (len<1) {break;}
                        s1.append(std::to_string(sequence[w]));
                        s1.append(1u,'-');
                    }
                } // depoly condition
            } // pt t
        }  // pt s
}	

void depolymerize_simple(std::vector<int> &bond_num, std::vector<std::vector<int>> &matrix, double ranr, double p_minus, int bead_number, int &dimer_counter, std::vector<std::string> &chain_atom_ID_list){ 
    for(int s=0; s<bead_number; ++s){
        for(int t=s+1; t<bead_number; ++t){
            ranr = ((double)rand()/(double)RAND_MAX);
            if ((ranr<=p_minus) and ((bond_num[s]==1 and bond_num[t]==1) or (bond_num[s]==1 and bond_num[t]==2) or (bond_num[s]==2 and bond_num[t]==1)) and (matrix[s][t]==1 and matrix[t][s]==1) ){
                matrix[s][t] = 0;
                matrix[t][s] = 0;
            	bond_num[s] =0;
            	bond_num[t] =0;
		dimer_counter -=1;
                backward_overlap(s, t, chain_atom_ID_list);
            }
        }//particle s
    }//particle t
}

void depolymerize(std::vector<double> &x_coor, std::vector<double> &y_coor, std::vector<double> &z_coor,  std::vector<int> &bond_num, std::vector<std::vector<int>> &matrix, std::vector<std::vector<int>> &mask_end, std::vector<int> &sequence, double ranr, int step_growth, int chain_growth, double dist, double box_size, double p_minus, int begin_index, int bead_number, bool in_vector, int next_index, std::vector<std::string> &chain_atom_ID_list, std::vector<int> &identities, std::vector<double> &ranr_x, std::vector<double> &ranr_y, std::vector<double> &ranr_z){ 

        for(int s=0; s<bead_number; ++s){
            for(int t=s+1; t<bead_number; ++t){
                ranr = ((double)rand()/(double)RAND_MAX);
                dist = sqrt(pow(((x_coor[s]-x_coor[t]) - round((x_coor[s]-x_coor[t])/box_size)*box_size), 2)+pow(((y_coor[s]-y_coor[t]) - round((y_coor[s]-y_coor[t])/box_size)*box_size), 2)+pow(((z_coor[s]-z_coor[t]) - round((z_coor[s]-z_coor[t])/box_size)*box_size), 2));
                // depolymerization anywhere!
                if ((step_growth==1) and (ranr<=p_minus) and ((bond_num[s]==2 and bond_num[t]==2) or (bond_num[s]==2 and bond_num[t]==1) or (bond_num[t]==2 and bond_num[s]==1) or (bond_num[s]==1 and bond_num[t]==1)) and (matrix[s][t]==1 and matrix[t][s]==1) ){
            	    if (bond_num[s]==2 and bond_num[t]==2){
                        backward_overlap(s, t, chain_atom_ID_list);
		        begin_index = s;
		        delete_middle_end(sequence, matrix, mask_end, bond_num, begin_index, bead_number, in_vector, next_index);
            	        bond_num[s] -=1;
            	        bond_num[t] -=1;
		        matrix[s][t] = 0;
		        matrix[t][s] = 0;
		        begin_index = s;
		        find_middle_end(sequence, matrix, mask_end, bond_num, begin_index, bead_number, in_vector, next_index);
		        begin_index = t;
		        find_middle_end(sequence, matrix, mask_end, bond_num, begin_index, bead_number, in_vector, next_index);

	            }
		    else if (bond_num[s]==2 and bond_num[t]==1){
                        backward_overlap(s, t, chain_atom_ID_list);
            	        begin_index = t;
                        delete_end(sequence, matrix, mask_end, bond_num, begin_index, bead_number, in_vector, next_index);
                        matrix[s][t] = 0;
                        matrix[t][s] = 0;
            	        begin_index = s; 
            	        bond_num[s] -=1;
            	        bond_num[t] -=1;
                        find_end(sequence, matrix, mask_end, bond_num, begin_index, bead_number, in_vector, next_index);
            	    }
		    else if (bond_num[t]==2 and bond_num[s]==1){
                        backward_overlap(s, t, chain_atom_ID_list);
            	        begin_index = s;
                        delete_end(sequence, matrix, mask_end, bond_num, begin_index, bead_number, in_vector, next_index);
                        matrix[s][t] = 0;
                        matrix[t][s] = 0;
            	        begin_index = t; 
            	        bond_num[s] -=1;
            	        bond_num[t] -=1;
                        find_end(sequence, matrix, mask_end, bond_num, begin_index, bead_number, in_vector, next_index);
            	    }
		    else if (bond_num[s]==1 and bond_num[t]==1){
                        backward_overlap(s, t, chain_atom_ID_list);
                        matrix[s][t] = 0;
                        matrix[t][s] = 0;
            	        bond_num[s] -=1;
            	        bond_num[t] -=1;
            	    }
		    /////////////////////////////////////
                    // VISUALIZATION OF POLYMER GROWTH //
                    /////////////////////////////////////
	            std::string s1;
                    for (int w=0; w<sequence.size(); ++w){
                        int len = sequence.size();
                        if (len<1) {break;}
                        s1.append(std::to_string(sequence[w]));
                        s1.append(1u,'-');
                    }
                } // anywhere depoly condition



                // end_chain depolymerization
                if ((chain_growth==1) and (ranr<=p_minus) and ((bond_num[s]==2 and bond_num[t]==1) or (bond_num[t]==2 and bond_num[s]==1) or (bond_num[s]==1 and bond_num[t]==1)) and (matrix[s][t]==1 and matrix[t][s]==1) ){
            	    if (bond_num[s]==2){
                        backward_overlap(s, t, chain_atom_ID_list);
            	        // particle t is released, delete masked element between t and other end of polymer
            	        begin_index = t;
                        delete_end(sequence, matrix, mask_end, bond_num, begin_index, bead_number, in_vector, next_index);
                        matrix[s][t] = 0;
                        matrix[t][s] = 0;
            	        begin_index = s; 
            	        bond_num[s] -=1;
            	        bond_num[t] -=1;
                        find_end(sequence, matrix, mask_end, bond_num, begin_index, bead_number, in_vector, next_index);
            	    }
            	    if (bond_num[t]==2){
                        backward_overlap(s, t, chain_atom_ID_list);
            	        // particle t is released, delete masked element between t and other end of polymer
            	        begin_index = s;
                        delete_end(sequence, matrix, mask_end, bond_num, begin_index, bead_number, in_vector, next_index);
                        matrix[s][t] = 0;
                        matrix[t][s] = 0;
            	        begin_index = t; 
            	        bond_num[s] -=1;
            	        bond_num[t] -=1;
                        find_end(sequence, matrix, mask_end, bond_num, begin_index, bead_number, in_vector, next_index);
            	    }
            	    if (bond_num[s]==1 and bond_num[t]==1){
                        backward_overlap(s, t, chain_atom_ID_list);
                        matrix[s][t] = 0;
                        matrix[t][s] = 0;
            	        bond_num[s] =0;
            	        bond_num[t] =0;
            	    }
                } // end_chain depoly condition
            } // pt t
        }  // pt s
}	

int depolymerize_with_multiple_probabilities(int addParticles, int step_growth, int chain_growth, double ranr, std::vector<std::string> &to_delete, int count_added_particles, double kbT, int &bead_number, double m, std::vector<double> &x_coor, std::vector<double> &y_coor, std::vector<double> &z_coor, std::vector<double> &new_x_coor, std::vector<double> &new_y_coor, std::vector<double> &new_z_coor, std::vector<double> &x_coor_pbc_whole, std::vector<double> &y_coor_pbc_whole, std::vector<double> &z_coor_pbc_whole, std::vector<double> &vx, std::vector<double> &vy, std::vector<double> &vz , std::vector<double> &new_vx, std::vector<double> &new_vy, std::vector<double> &new_vz, std::vector<double> &harter_x, std::vector<double> &harter_y, std::vector<double> &harter_z, std::vector<int> &bond_num, std::vector<std::string> &chain_atom_ID_list, std::vector<std::vector<int>> &matrix_cpy, std::vector<std::vector<double>> &matrix_harter_x, std::vector<std::vector<double>> &matrix_harter_y, std::vector<std::vector<double>> &matrix_harter_z, std::vector<std::vector<double>> &matrix_bond_r, std::vector<std::vector<int>> &mask_end, std::vector<std::vector<int>> &matrix, std::vector<std::vector<double>>  &probability_matrix_plus, std::vector<std::vector<double>>  &probability_matrix_minus, std::vector<int> &identities, std::vector<int> &sequence, double dist, double box_size, int begin_index, bool in_vector, int next_index, std::vector<double> &ranr_x, std::vector<double> &ranr_y, std::vector<double> &ranr_z){


    int count_deleted_particles = 0;
    std::vector<int> accumulate_delete;
    accumulate_delete.clear();
    std::vector<int> vctr_id;
    std::vector<std::string> chain_list_tmp;
    chain_list_tmp.clear();
    for(int m=0; m<bead_number; ++m){
        vctr_id.push_back(m);
    }
    chain_list_tmp.clear();
        for(int s=0; s<bead_number; ++s){
            for(int t=s+1; t<bead_number; ++t){
                ranr = ((double)rand()/(double)RAND_MAX);
                dist = sqrt(pow(((x_coor[s]-x_coor[t]) - round((x_coor[s]-x_coor[t])/box_size)*box_size), 2)+pow(((y_coor[s]-y_coor[t]) - round((y_coor[s]-y_coor[t])/box_size)*box_size), 2)+pow(((z_coor[s]-z_coor[t]) - round((z_coor[s]-z_coor[t])/box_size)*box_size), 2));
                // depolymerization anywhere!
                if ((step_growth==1) and (ranr<=probability_matrix_minus[s][t]) and ((bond_num[s]==2 and bond_num[t]==2) or (bond_num[s]==2 and bond_num[t]==1) or (bond_num[t]==2 and bond_num[s]==1) or (bond_num[s]==1 and bond_num[t]==1)) and (matrix[s][t]==1 and matrix[t][s]==1) ){
            	    std::cout << "disconnecting particles: " << s << " and " << t << std::endl;
            	    if (bond_num[s]==2 and bond_num[t]==2){
		        std::cout << " breaking in the middle between particles s: " << s << " and " << t << std::endl;
                        step_growth_backward_overlap_delete(s, t, chain_atom_ID_list);
		        begin_index = s;
		        delete_middle_end(sequence, matrix, mask_end, bond_num, begin_index, bead_number, in_vector, next_index);
            	        bond_num[s] -=1;
            	        bond_num[t] -=1;
		        matrix[s][t] = 0;
		        matrix[t][s] = 0;
		        begin_index = s;
		        find_middle_end(sequence, matrix, mask_end, bond_num, begin_index, bead_number, in_vector, next_index);
		        begin_index = t;
		        find_middle_end(sequence, matrix, mask_end, bond_num, begin_index, bead_number, in_vector, next_index);

	            }
		    else if (bond_num[s]==2 and bond_num[t]==1){
		        std::cout << " breaking in the end, s remains the end, t dettaches, s: " << s << " t: " << t << std::endl;
                        step_growth_backward_overlap_delete(s, t, chain_atom_ID_list);
            	        begin_index = t;
                        delete_end(sequence, matrix, mask_end, bond_num, begin_index, bead_number, in_vector, next_index);
                        matrix[s][t] = 0;
                        matrix[t][s] = 0;
            	        begin_index = s; 
            	        bond_num[s] -=1;
            	        bond_num[t] -=1;
                        find_end(sequence, matrix, mask_end, bond_num, begin_index, bead_number, in_vector, next_index);
			if(addParticles){
			    count_deleted_particles += 1;
			    accumulate_delete.push_back(t);
                         }
            	    }
		    else if (bond_num[t]==2 and bond_num[s]==1){
                        step_growth_backward_overlap_delete(s, t, chain_atom_ID_list);
            	        begin_index = s;
                        delete_end(sequence, matrix, mask_end, bond_num, begin_index, bead_number, in_vector, next_index);
                        matrix[s][t] = 0;
                        matrix[t][s] = 0;
            	        begin_index = t; 
            	        bond_num[s] -=1;
            	        bond_num[t] -=1;
                        find_end(sequence, matrix, mask_end, bond_num, begin_index, bead_number, in_vector, next_index);
			if(addParticles){
			    count_deleted_particles += 1;
			    accumulate_delete.push_back(s);
                         }
            	    }
		    else if (bond_num[s]==1 and bond_num[t]==1){
                        step_growth_backward_overlap_delete(s, t, chain_atom_ID_list);
                        matrix[s][t] = 0;
                        matrix[t][s] = 0;
            	        bond_num[s] -=1;
            	        bond_num[t] -=1;
			if(addParticles){
			    count_deleted_particles += 2;
			    accumulate_delete.push_back(s);
			    accumulate_delete.push_back(t);
                         }
            	    }
		    /////////////////////////////////////
                    // VISUALIZATION OF POLYMER GROWTH //
                    /////////////////////////////////////
	            std::string s1;
                    for (int w=0; w<sequence.size(); ++w){
                        int len = sequence.size();
                        if (len<1) {break;}
                        s1.append(std::to_string(sequence[w]));
                        s1.append(1u,'-');
                    }

                } // anywhere depoly condition

                // end_chain depolymerization
                if ((chain_growth==1) and (ranr<=probability_matrix_minus[s][t]) and ((bond_num[s]==2 and bond_num[t]==1) or (bond_num[t]==2 and bond_num[s]==1) or (bond_num[s]==1 and bond_num[t]==1)) and (matrix[s][t]==1 and matrix[t][s]==1) ){
		    std::cout << " end chain depolymerization between s: " << s << " and t: " << t << std::endl;
            	    if (bond_num[s]==2){
                        chain_growth_backward_overlap_delete(s, t, chain_atom_ID_list);
            	        // particle t is released, delete masked element between t and other end of polymer
            	        begin_index = t;
                        delete_end(sequence, matrix, mask_end, bond_num, begin_index, bead_number, in_vector, next_index);
                        matrix[s][t] = 0;
                        matrix[t][s] = 0;
            	        begin_index = s; 
            	        bond_num[s] -=1;
            	        bond_num[t] -=1;
                        find_end(sequence, matrix, mask_end, bond_num, begin_index, bead_number, in_vector, next_index);
			if(addParticles){
			    count_deleted_particles += 1;
			    accumulate_delete.push_back(t);
                         }
            	    }
            	    if (bond_num[t]==2){
                        chain_growth_backward_overlap_delete(s, t, chain_atom_ID_list);
            	        // particle s is released, delete masked element between t and other end of polymer
            	        begin_index = s;
                        delete_end(sequence, matrix, mask_end, bond_num, begin_index, bead_number, in_vector, next_index);
                        matrix[s][t] = 0;
                        matrix[t][s] = 0;
            	        begin_index = t; 
            	        bond_num[s] -=1;
            	        bond_num[t] -=1;
                        find_end(sequence, matrix, mask_end, bond_num, begin_index, bead_number, in_vector, next_index);
			if(addParticles){
			    count_deleted_particles += 1;
			    accumulate_delete.push_back(s);
                         }
            	    }
            	    if (bond_num[s]==1 and bond_num[t]==1){
                        chain_growth_backward_overlap_delete(s, t, chain_atom_ID_list);
                        matrix[s][t] = 0;
                        matrix[t][s] = 0;
            	        bond_num[s] =0;
            	        bond_num[t] =0;
			if(addParticles){
			    count_deleted_particles += 2;
			    accumulate_delete.push_back(s);
			    accumulate_delete.push_back(t);
                         }
            	    }
                } // end_chain depoly condition
            } // pt t
        }  // pt s
    ///////////////////
    // DELETING PART //
    ///////////////////
    std::sort(accumulate_delete.begin(), accumulate_delete.end(), std::greater<int>());
    for (int i=0; i<accumulate_delete.size(); ++i){
        delete_index_int_vector(accumulate_delete[i], vctr_id);
    }

    for (int w=0; w<chain_atom_ID_list.size(); ++w){
        std::string s2;
        std::istringstream is( chain_atom_ID_list[w] );
        int n;
	s2.clear();
        while( is >> n ) {
		int printn = n;
            for (int i=0; i<vctr_id.size(); ++i){
		 if(vctr_id[i] == n){
                     s2.append(std::to_string(i));
                     s2.append(1u,' ');
		 }
	    }
        }
	chain_list_tmp.push_back(s2);
    }
    std::sort(accumulate_delete.rbegin(), accumulate_delete.rend());
    for(int i=0; i<accumulate_delete.size(); ++i){
        delete_particle(accumulate_delete[i], to_delete, count_added_particles, box_size, kbT, bead_number, m, x_coor, y_coor, z_coor, new_x_coor, new_y_coor, new_z_coor, x_coor_pbc_whole, y_coor_pbc_whole, z_coor_pbc_whole, vx, vy, vz, new_vx, new_vy, new_vz, harter_x, harter_y, harter_z, bond_num, chain_atom_ID_list, matrix_cpy, matrix_harter_x, matrix_harter_y, matrix_harter_z, matrix_bond_r, mask_end, matrix, probability_matrix_plus, probability_matrix_minus, identities, ranr_x, ranr_y, ranr_z);
    }
    chain_atom_ID_list = chain_list_tmp;
    return count_deleted_particles;
}	
