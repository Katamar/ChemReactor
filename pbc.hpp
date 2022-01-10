void pbc_write(std::vector<double> &x_coor, std::vector<double> &y_coor, std::vector<double> &z_coor, std::vector<double> &x_coor_pbc_whole, std::vector<double> &y_coor_pbc_whole, std::vector<double> &z_coor_pbc_whole, int bead_number, double box_size){
    for(int p=0; p<bead_number; ++p){
        x_coor_pbc_whole[p] = x_coor[p] - round(x_coor[p]/box_size)*box_size;
        y_coor_pbc_whole[p] = y_coor[p] - round(y_coor[p]/box_size)*box_size;
        z_coor_pbc_whole[p] = z_coor[p] - round(z_coor[p]/box_size)*box_size;
    }
}
