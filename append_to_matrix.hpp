void append_double_to_matrix(std::vector<std::vector<double>> &matrix, double numeric_fill){

    // define row
    std::vector<int> vector1(matrix, numeric_fill);
    // add row
    matrix.push_back(vector1);
    // add column
    for (int k = 0 ; k < matrix.size() ; ++k){
        matrix[k].push_back(numeric_fill);
    }

}

void append_int_to_matrix(std::vector<std::vector<int>> &matrix, int numeric_fill){

    // define row
    std::vector<int> vector1(matrix, numeric_fill);
    // add row
    matrix.push_back(vector1);
    // add column
    for (int k = 0 ; k < matrix.size() ; ++k){
        matrix[k].push_back(numeric_fill);
    }

}

void add_particle(double kbT, double &bead_number, std::vector<double> &x_coor, std::vector<double> &y_coor, std::vector<double> &z_coor, std::vector<double> &new_x_coor, std::vector<double> &new_y_coor, std::vector<double> &new_z_coor, std::vector<double> &x_coor_pbc_whole, std::vector<double> &y_coor_pbc_whole, std::vector<double> &z_coor_pbc_whole, std::vector<double> &vx, std::vector<double> &vy, std::vector<double> &vz, std::vector<double> &new_vx, std::vector<double> &new_vy, std::vector<double> &new_vz, std::vector<double> &harter_x, std::vector<double> &harter_y, std::vector<double> &harter_z, std::vector<double> &x_coor_pbc_whole, std::vector<double> &y_coor_pbc_whole, std::vector<double> &z_coor_pbc_whole, std::vector<double> &bond_num, std::vector<std::string> &chain_atom_ID_list, std::vector<int> &identities){
    //////////////////////
    // ADDING PARTICLES //
    //////////////////////
    bead_number += 1;
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
    x_coor_pbc_whole.push_back(0.0);
    y_coor_pbc_whole.push_back(0.0);
    z_coor_pbc_whole.push_back(0.0);
    bond_num.push_back(0);
}
