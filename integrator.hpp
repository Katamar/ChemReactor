void integrate_pt1(std::vector<double> &x_coor, std::vector<double> &y_coor, std::vector<double> &z_coor, std::vector<double> &new_x_coor, std::vector<double> &new_y_coor, std::vector<double> &new_z_coor, std::vector<double> &vx, std::vector<double> &vy, std::vector<double> &vz, std::vector<double> &new_vx, std::vector<double> &new_vy, std::vector<double> &new_vz,std::vector<double> &ranr_x, std::vector<double> &ranr_y, std::vector<double> &ranr_z, std::vector<double> &harter_x, std::vector<double> &harter_y, std::vector<double> &harter_z, double gamma_coeff, double dt, double m, double kbT, int bead_number){
for (int i=0; i<bead_number; ++i){
    new_vx[i] = vx[i] - gamma_coeff*dt/2.0*vx[i] + dt/(2.0*m)*0.6022*(harter_x[i]) + dt/(2.0*m)*sqrt((2.0*0.6022*kbT*gamma_coeff*m)/dt)*ranr_x[i];
    //new_vx[i] = vx[i] - gamma_coeff*dt/2.0*vx[i] + dt/(2.0*m)*0.6022*(harter_x[i]+vdW_x[i]) + dt/(2.0*m)*sqrt((2.0*0.6022*kbT*gamma_coeff*m)/dt)*ranr_x[i];
    //new_vx[i] = vx[i] - gamma_coeff*dt/2.0*vx[i] + dt/(2.0*m)*0.6022*(harter_x[i]+angle_x[i]+vdW_x[i]) + dt/(2.0*m)*sqrt((2.0*0.6022*kbT*gamma_coeff*m)/dt)*ranr_x[i];
    new_x_coor[i]=x_coor[i]+new_vx[i]*dt;
    x_coor[i]=new_x_coor[i];
    vx[i]=new_vx[i];
    
    new_vy[i] = vy[i] - gamma_coeff*dt/2.0*vy[i] + dt/(2.0*m)*0.6022*(harter_y[i]) + dt/(2.0*m)*sqrt((2.0*0.6022*kbT*gamma_coeff*m)/dt)*ranr_y[i];
    //new_vy[i] = vy[i] - gamma_coeff*dt/2.0*vy[i] + dt/(2.0*m)*0.6022*(harter_y[i]+vdW_y[i]) + dt/(2.0*m)*sqrt((2.0*0.6022*kbT*gamma_coeff*m)/dt)*ranr_y[i];
    //new_vy[i] = vy[i] - gamma_coeff*dt/2.0*vy[i] + dt/(2.0*m)*0.6022*(harter_y[i]+angle_y[i]+vdW_y[i]) + dt/(2.0*m)*sqrt((2.0*0.6022*kbT*gamma_coeff*m)/dt)*ranr_y[i];
    new_y_coor[i]=y_coor[i]+new_vy[i]*dt;
    y_coor[i]=new_y_coor[i];
    vy[i]=new_vy[i];

    new_vz[i] = vz[i] - gamma_coeff*dt/2.0*vz[i] + dt/(2.0*m)*0.6022*(harter_z[i]) + dt/(2.0*m)*sqrt((2.0*0.6022*kbT*gamma_coeff*m)/dt)*ranr_z[i];
    //new_vz[i] = vz[i] - gamma_coeff*dt/2.0*vz[i] + dt/(2.0*m)*0.6022*(harter_z[i]+vdW_z[i]) + dt/(2.0*m)*sqrt((2.0*0.6022*kbT*gamma_coeff*m)/dt)*ranr_z[i];
    //new_vz[i] = vz[i] - gamma_coeff*dt/2.0*vz[i] + dt/(2.0*m)*0.6022*(harter_z[i]+angle_z[i]+vdW_z[i]) + dt/(2.0*m)*sqrt((2.0*0.6022*kbT*gamma_coeff*m)/dt)*ranr_z[i];
    new_z_coor[i]=z_coor[i]+new_vz[i]*dt;
    z_coor[i]=new_z_coor[i];
    vz[i]=new_vz[i];
    }
}

void integrate_pt2(std::vector<double> &x_coor, std::vector<double> &y_coor, std::vector<double> &z_coor, std::vector<double> &new_x_coor, std::vector<double> &new_y_coor, std::vector<double> &new_z_coor, std::vector<double> &vx, std::vector<double> &vy, std::vector<double> &vz, std::vector<double> &new_vx, std::vector<double> &new_vy, std::vector<double> &new_vz,std::vector<double> &ranr_x, std::vector<double> &ranr_y, std::vector<double> &ranr_z, std::vector<double> &harter_x, std::vector<double> &harter_y, std::vector<double> &harter_z, double gamma_coeff, double dt, double m, double kbT, int bead_number){
    for (int i=0; i<bead_number; ++i){
        new_vx[i] = vx[i] - gamma_coeff*dt/2.0*vx[i] + dt/(2.0*m)*0.6022*(harter_x[i]) + dt/(2.0*m)*sqrt((2.0*0.6022*kbT*gamma_coeff*m)/dt)*ranr_x[i];
        //new_vx[i] = vx[i] - gamma_coeff*dt/2.0*vx[i] + dt/(2.0*m)*0.6022*(harter_x[i]+vdW_x[i]) + dt/(2.0*m)*sqrt((2.0*0.6022*kbT*gamma_coeff*m)/dt)*ranr_x[i];
        //new_vx[i] = vx[i] - gamma_coeff*dt/2.0*vx[i] + dt/(2.0*m)*0.6022*(harter_x[i]+angle_x[i]+vdW_x[i]) + dt/(2.0*m)*sqrt((2.0*0.6022*kbT*gamma_coeff*m)/dt)*ranr_x[i];
        vx[i]=new_vx[i];
    
        new_vy[i] = vy[i] - gamma_coeff*dt/2.0*vy[i] + dt/(2.0*m)*0.6022*(harter_y[i]) + dt/(2.0*m)*sqrt((2.0*0.6022*kbT*gamma_coeff*m)/dt)*ranr_y[i];
	//new_vy[i] = vy[i] - gamma_coeff*dt/2.0*vy[i] + dt/(2.0*m)*0.6022*(harter_y[i]+vdW_y[i]) + dt/(2.0*m)*sqrt((2.0*0.6022*kbT*gamma_coeff*m)/dt)*ranr_y[i];
        //new_vy[i] = vy[i] - gamma_coeff*dt/2.0*vy[i] + dt/(2.0*m)*0.6022*(harter_y[i]+angle_y[i]+vdW_y[i]) + dt/(2.0*m)*sqrt((2.0*0.6022*kbT*gamma_coeff*m)/dt)*ranr_y[i];
        vy[i]=new_vy[i];
    
        new_vz[i] = vz[i] - gamma_coeff*dt/2.0*vz[i] + dt/(2.0*m)*0.6022*(harter_z[i]) + dt/(2.0*m)*sqrt((2.0*0.6022*kbT*gamma_coeff*m)/dt)*ranr_z[i];
        //new_vz[i] = vz[i] - gamma_coeff*dt/2.0*vz[i] + dt/(2.0*m)*0.6022*(harter_z[i]+vdW_z[i]) + dt/(2.0*m)*sqrt((2.0*0.6022*kbT*gamma_coeff*m)/dt)*ranr_z[i];
        //new_vz[i] = vz[i] - gamma_coeff*dt/2.0*vz[i] + dt/(2.0*m)*0.6022*(harter_z[i]+angle_z[i]+vdW_z[i]) + dt/(2.0*m)*sqrt((2.0*0.6022*kbT*gamma_coeff*m)/dt)*ranr_z[i];
        vz[i]=new_vz[i];
    }
}
