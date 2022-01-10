
void CFSTR_add(int &CFSTR_add_number_1, int &CFSTR_add_number_2, std::vector<std::vector<int>> &matrix, std::vector<int> &bond_num, int &bead_number, double ranr, double dist, double binding_cutoff, std::vector<std::vector<double>>  &probability_matrix_plus, std::vector<std::vector<double>>  &probability_matrix_minus, std::vector<double> &x_coor, std::vector<double> &y_coor, std::vector<double> &z_coor, double box_size, std::vector<std::string> &to_delete, std::vector<int> &sequence, std::vector<std::string> &chain_atom_ID_list, std::vector<int> &identities, std::vector<std::vector<int>> &matrix_cpy, std::vector<std::vector<double>> &matrix_harter_x, std::vector<std::vector<double>> &matrix_harter_y, std::vector<std::vector<double>> &matrix_harter_z, std::vector<std::vector<double>> &matrix_bond_r, std::vector<std::vector<int>> &mask_end, std::vector<double> &vx, std::vector<double> &vy, std::vector<double> &vz, std::vector<double> &new_vx, std::vector<double> &new_vy, std::vector<double> &new_vz, std::vector<double> &x_coor_pbc_whole, std::vector<double> &y_coor_pbc_whole, std::vector<double> &z_coor_pbc_whole, std::vector<double> &new_x_coor, std::vector<double> &new_y_coor, std::vector<double> &new_z_coor, double kbT, std::vector<double> &harter_x, std::vector<double> &harter_y, std::vector<double> &harter_z, double m, std::vector<double> proportions, std::vector<double> p_plus_probs, std::vector<double> p_minus_probs, int &count_added_particles, std::vector<double> &ranr_x, std::vector<double> &ranr_y, std::vector<double> &ranr_z){

int num1 = int(CFSTR_add_number_1);
int num2 = int(CFSTR_add_number_2);
// add particles of type 1
for(int k=0; k<num1; ++k){
    count_added_particles += 1;
    identities.push_back(1);
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

// add particles of type 2
for(int k=0; k<num2; ++k){
    count_added_particles += 1;
    identities.push_back(2);
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

}


//CFSTR_remove_monomer_number
void CFSTR_remove(double ranr, std::vector<std::string> &to_delete, int count_added_particles, double box_size, double kbT, int &bead_number, double m, std::vector<double> &x_coor, std::vector<double> &y_coor, std::vector<double> &z_coor, std::vector<double> &new_x_coor, std::vector<double> &new_y_coor, std::vector<double> &new_z_coor, std::vector<double> &x_coor_pbc_whole, std::vector<double> &y_coor_pbc_whole, std::vector<double> &z_coor_pbc_whole, std::vector<double> &vx, std::vector<double> &vy, std::vector<double> &vz , std::vector<double> &new_vx, std::vector<double> &new_vy, std::vector<double> &new_vz, std::vector<double> &harter_x, std::vector<double> &harter_y, std::vector<double> &harter_z, std::vector<int> &bond_num, std::vector<std::string> &chain_atom_ID_list, std::vector<std::vector<int>> &matrix_cpy, std::vector<std::vector<double>> &matrix_harter_x, std::vector<std::vector<double>> &matrix_harter_y, std::vector<std::vector<double>> &matrix_harter_z, std::vector<std::vector<double>> &matrix_bond_r, std::vector<std::vector<int>> &mask_end, std::vector<std::vector<int>> &matrix, std::vector<std::vector<double>>  &probability_matrix_plus, std::vector<std::vector<double>>  &probability_matrix_minus, std::vector<int> &identities, std::vector<double> &ranr_x, std::vector<double> &ranr_y, std::vector<double> &ranr_z, int &count_deleted_particles, double CFSTR_remove_p){
std::vector<int> sequence;
sequence.clear();
std::vector<int> unique_remove_list;
std::vector<int> vctr_id;
std::vector<std::string> chain_list_tmp;
chain_list_tmp.clear();
for(int m=0; m<bead_number; ++m){
    vctr_id.push_back(m);
}
int index;
int begin_index;
int next_index;
bool in_vector;
double ranrt;

for(int k=0; k<bead_number; ++k){
    ranrt = ((double)rand()/(double)RAND_MAX);
    if (ranrt <= CFSTR_remove_p){
        if(bond_num[k] == 0){
	    unique_remove_list.push_back(k);
        }
        if(bond_num[k] == 1){
            begin_index = k;
	    find_end(sequence, matrix, mask_end, bond_num, begin_index, bead_number, in_vector, next_index);
	    if(next_index>k){
	    unique_remove_list.push_back(k);
	    unique_remove_list.push_back(next_index);

            }
        }
    }

}


if(unique_remove_list.size()>0){
    std::sort(unique_remove_list.begin(), unique_remove_list.end(), std::greater<int>());
}


for (int i=0; i<unique_remove_list.size(); ++i){
    delete_index_int_vector(unique_remove_list[i], vctr_id);
}
for (int w=0; w<chain_atom_ID_list.size(); ++w){
    std::string s2;
    std::istringstream is( chain_atom_ID_list[w] );
    int n;
    s2.clear();
    while( is >> n ) {
        for (int i=0; i<vctr_id.size(); ++i){
    	    if(vctr_id[i] == n){
                    s2.append(std::to_string(i));
                    s2.append(1u,' ');
    	    }
        }
    }
    if(s2.empty() == 0){
        chain_list_tmp.push_back(s2);
    }
}
chain_atom_ID_list = chain_list_tmp;


// delete particles
if(unique_remove_list.size()>0){
    for(int k=0; k<unique_remove_list.size(); ++k){
        index = unique_remove_list[k];
	count_deleted_particles +=1;
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
        delete_rowcolumn_int_matrix(index, matrix);
        delete_rowcolumn_int_matrix(index, matrix_cpy);
        delete_rowcolumn_int_matrix(index, mask_end);
        delete_rowcolumn_double_matrix(index, matrix_harter_x);
        delete_rowcolumn_double_matrix(index, matrix_harter_y);
        delete_rowcolumn_double_matrix(index, matrix_harter_z);
        delete_rowcolumn_double_matrix(index, matrix_bond_r);
        delete_rowcolumn_double_matrix(index, probability_matrix_plus);
        delete_rowcolumn_double_matrix(index, probability_matrix_minus);
    }
}



}
