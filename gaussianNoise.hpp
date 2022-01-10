double generateGaussianNoise(double mu, double sigma, double u1, double u2, double two_pi){
    u1 = ((double)rand()/(double)RAND_MAX);
    u2 = ((double)rand()/(double)RAND_MAX);
    if (u1==0) {u1=0.001;}
    if (u1==1) {u1=0.999;}
    return (sqrt(-2.0 * log(u1)) * cos(two_pi * u2) ) * sigma + mu;
}
