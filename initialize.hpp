void init_positions(std::vector<double> &x_coor, std::vector<double> &y_coor, std::vector<double> &z_coor, std::vector<double> &new_x_coor, std::vector<double> &new_y_coor, std::vector<double> &new_z_coor, std::vector<double> &x_coor_pbc_whole, std::vector<double> &y_coor_pbc_whole, std::vector<double> &z_coor_pbc_whole, double box_size, int bead_number){
    for (int i=0; i<bead_number; ++i){
        x_coor.push_back(((double)rand()/(double)RAND_MAX)*box_size);
        y_coor.push_back(((double)rand()/(double)RAND_MAX)*box_size);
        z_coor.push_back(((double)rand()/(double)RAND_MAX)*box_size);
	
    }
}

void init_velocities(std::vector<double> &vx, std::vector<double> &vy, std::vector<double> &vz, std::vector<double> &new_vx, std::vector<double> &new_vy, std::vector<double> &new_vz, int bead_number, double kbT, double m){
    for (int i=0; i<bead_number; ++i){
	vx.push_back(generateGaussianNoise(0.0, 1.0, 0.5, 0.5, 6.283185307179586)*sqrt(kbT*0.6022/m));
	vy.push_back(generateGaussianNoise(0.0, 1.0, 0.5, 0.5, 6.283185307179586)*sqrt(kbT*0.6022/m));
	vz.push_back(generateGaussianNoise(0.0, 1.0, 0.5, 0.5, 6.283185307179586)*sqrt(kbT*0.6022/m));
	new_vx.push_back(0.0);
	new_vy.push_back(0.0);
	new_vz.push_back(0.0);
    }
}

void init_ranr(std::vector<double> &ranr_x, std::vector<double> &ranr_y, std::vector<double> &ranr_z, int bead_number){
    for (int i=0; i<bead_number; ++i){
	ranr_x.push_back(0.0);
	ranr_y.push_back(0.0);
	ranr_z.push_back(0.0);
    }  
}

void update_ranr(std::vector<double> &ranr_x, std::vector<double> &ranr_y, std::vector<double> &ranr_z, int bead_number){
    for (int i=0; i<bead_number; ++i){
	ranr_x[i] = generateGaussianNoise(0.0, 1.0, 0.5, 0.5, 6.283185307179586);
	ranr_y[i] = generateGaussianNoise(0.0, 1.0, 0.5, 0.5, 6.283185307179586);
	ranr_z[i] = generateGaussianNoise(0.0, 1.0, 0.5, 0.5, 6.283185307179586);
    }
}

void initialize_arrays(std::vector<double> &harter_x, std::vector<double> &harter_y, std::vector<double> &harter_z, std::vector<double> &new_x_coor, std::vector<double> &new_y_coor, std::vector<double> &new_z_coor, std::vector<double> &x_coor_pbc_whole, std::vector<double> &y_coor_pbc_whole, std::vector<double> &z_coor_pbc_whole, std::vector<int> &bond_num, int bead_number, std::vector<std::string> &chain_atom_ID_list,  std::vector<double> &vdW_x, std::vector<double> &vdW_y, std::vector<double> &vdW_z, std::vector<int> &identities){
    for (int i=0; i<bead_number; ++i){
        harter_x.push_back(0.0);
        harter_y.push_back(0.0);
        harter_z.push_back(0.0);
        new_x_coor.push_back(0.0);
        new_y_coor.push_back(0.0);
        new_z_coor.push_back(0.0);
	vdW_x.push_back(0.0);
	vdW_y.push_back(0.0);
	vdW_z.push_back(0.0);
        x_coor_pbc_whole.push_back(0.0);
        y_coor_pbc_whole.push_back(0.0);
        z_coor_pbc_whole.push_back(0.0);
        bond_num.push_back(0);
	chain_atom_ID_list.push_back(std::to_string(i));
	identities.push_back(0);
    }
}

