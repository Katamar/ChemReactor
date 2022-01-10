void harmonic_force(std::vector<double> &x_coor, std::vector<double> &y_coor, std::vector<double> &z_coor, std::vector<double> &harter_x, std::vector<double> &harter_y, std::vector<double> &harter_z, std::vector<std::vector<int>> &matrix, std::vector<std::vector<int>> &mask_end, std::vector<std::vector<double>> &matrix_harter_x, std::vector<std::vector<double>> &matrix_harter_y, std::vector<std::vector<double>> &matrix_harter_z, std::vector<std::vector<double>> &matrix_bond_r, int bead_number, double harmonic_constant, double d0, double box_size){
    for(int j=0; j<bead_number; j++){
        harter_x[j] = 0;
        harter_y[j] = 0;
        harter_z[j] = 0;
    }
    for (auto &i : matrix_harter_x)
        std::fill(i.begin(), i.end(), 0);
    for (auto &i : matrix_harter_y)
        std::fill(i.begin(), i.end(), 0);
    for (auto &i : matrix_harter_z)
        std::fill(i.begin(), i.end(), 0);
    for (auto &i : matrix_bond_r)
        std::fill(i.begin(), i.end(), 0);
    for(int s=0; s<bead_number; ++s){
        for(int t=s+1; t<bead_number; ++t){
            matrix_harter_x[s][t] = (x_coor[s]-x_coor[t]) - round((x_coor[s]-x_coor[t])/box_size)*box_size;
            matrix_harter_y[s][t] = (y_coor[s]-y_coor[t]) - round((y_coor[s]-y_coor[t])/box_size)*box_size;
            matrix_harter_z[s][t] = (z_coor[s]-z_coor[t]) - round((z_coor[s]-z_coor[t])/box_size)*box_size;
            matrix_bond_r[s][t] = sqrt(pow(matrix_harter_x[s][t], 2)+pow(matrix_harter_y[s][t], 2)+pow(matrix_harter_z[s][t], 2));
            matrix_harter_x[s][t] *= (-1)*harmonic_constant*(matrix_bond_r[s][t]-d0)/matrix_bond_r[s][t]; 
            matrix_harter_y[s][t] *= (-1)*harmonic_constant*(matrix_bond_r[s][t]-d0)/matrix_bond_r[s][t]; 
            matrix_harter_z[s][t] *= (-1)*harmonic_constant*(matrix_bond_r[s][t]-d0)/matrix_bond_r[s][t];
            matrix_harter_x[t][s] = (-1)*matrix_harter_x[s][t]; 
            matrix_harter_y[t][s] = (-1)*matrix_harter_y[s][t]; 
            matrix_harter_z[t][s] = (-1)*matrix_harter_z[s][t]; 
            matrix_harter_x[s][t] *= matrix[s][t]; 
            matrix_harter_x[t][s] *= matrix[t][s];
            matrix_harter_y[s][t] *= matrix[s][t]; 
            matrix_harter_y[t][s] *= matrix[t][s];
            matrix_harter_z[s][t] *= matrix[s][t]; 
            matrix_harter_z[t][s] *= matrix[t][s];

            matrix_harter_x[s][t] *= mask_end[s][t]; 
            matrix_harter_x[t][s] *= mask_end[t][s];
            matrix_harter_y[s][t] *= mask_end[s][t]; 
            matrix_harter_y[t][s] *= mask_end[t][s];
            matrix_harter_z[s][t] *= mask_end[s][t]; 
            matrix_harter_z[t][s] *= mask_end[t][s];
        }
    } // end of double particle loop
    for(int i=0; i<bead_number; i++){
        for(int j=0; j<bead_number; j++){
            harter_x[i] += matrix_harter_x[i][j];
            harter_y[i] += matrix_harter_y[i][j];
            harter_z[i] += matrix_harter_z[i][j];
        }
    }
}


void lj_force(std::vector<double> &x_coor, std::vector<double> &y_coor, std::vector<double> &z_coor, std::vector<double> &vdW_x, std::vector<double> &vdW_y, std::vector<double> &vdW_z, int bead_number, double epsilon_lj, double sigma_lj, double rcut){ 

double rvector=0;

    for (int i=0; i<bead_number; ++i){
        vdW_x[i]=0;
        vdW_y[i]=0;
        vdW_z[i]=0;
        for (int j=0; j<bead_number; j++){
            if (j != i-1 && j != i+1 && j != i){
                rvector=sqrt(pow((x_coor[i]-x_coor[j]), 2)+pow((y_coor[i]-y_coor[j]), 2)+pow((z_coor[i]-z_coor[j]), 2));
                    if (rvector>rcut) { 
                        vdW_x[i]=0;
                        vdW_y[i]=0;
                        vdW_z[i]=0;
                    }
                    else{
                        vdW_x[i]+=24.0*epsilon_lj*(x_coor[i]-x_coor[j])*1/pow(rvector, 2)*(2*pow(sigma_lj, 12)*pow(rvector, -12)-pow(sigma_lj, 6)*pow(rvector, -6));
                        vdW_y[i]+=24.0*epsilon_lj*(y_coor[i]-y_coor[j])*1/pow(rvector, 2)*(2*pow(sigma_lj, 12)*pow(rvector, -12)-pow(sigma_lj, 6)*pow(rvector, -6));
                        vdW_z[i]+=24.0*epsilon_lj*(z_coor[i]-z_coor[j])*1/pow(rvector, 2)*(2*pow(sigma_lj, 12)*pow(rvector, -12)-pow(sigma_lj, 6)*pow(rvector, -6));
                    }
            }
        }
    }
}
