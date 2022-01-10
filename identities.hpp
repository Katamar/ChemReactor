void assign_identity(std::vector<int> &identities, int bead_number, double ranr, std::vector<double> proportions, std::vector<std::vector<double>>  &probability_matrix_plus, std::vector<std::vector<double>>  &probability_matrix_minus, std::vector<double> p_plus_probs, std::vector<double> p_minus_probs){
    std::vector<double> proportion_sums; 
    proportion_sums.push_back(0.0);
    double prop_ctr = 0;
    for(int i=0; i<proportions.size(); ++i){
	prop_ctr += proportions[i];
	proportion_sums.push_back(prop_ctr);
    }
    for(int i=0; i<bead_number; ++i){
        for(int j=0; j<proportion_sums.size()-1; ++j){
            if (i >= proportion_sums[j] and i < proportion_sums[j+1]){
		identities[i] = j+1;
	    }
	}
    }
    
    ////////////////////////////////
    // FORWARD PROBABILITY MATRIX //
    ////////////////////////////////

    int number_of_types = proportions.size();
    for(int s=0; s<bead_number; ++s){
        for(int t=s+1; t<bead_number; ++t){
            if (identities[s] == identities[t] and s != t){
	        probability_matrix_plus[s][t] = p_plus_probs[int(identities[s]-1)];
	        probability_matrix_plus[t][s] = p_plus_probs[int(identities[s]-1)];
	    }
	    else{
	        probability_matrix_plus[s][t] = p_plus_probs[int(identities[s]-1)];
	        probability_matrix_plus[t][s] = p_plus_probs[int(identities[s]-1)];
	    }
        }
    }
    for(int s=0; s<bead_number; ++s){
	probability_matrix_plus[s][s] = 0;
    }
    /////////////////////////////////
    // BACKWARD PROBABILITY MATRIX //
    /////////////////////////////////

    for(int s=0; s<bead_number; ++s){
        for(int t=s+1; t<bead_number; ++t){
            if (identities[s] == identities[t] and s != t){
	        probability_matrix_minus[s][t] = p_minus_probs[int(identities[s]-1)];
	        probability_matrix_minus[t][s] = p_minus_probs[int(identities[s]-1)];
	    }
	    else{
	        probability_matrix_minus[s][t] = p_minus_probs[int(identities[s]-1)];
	        probability_matrix_minus[t][s] = p_minus_probs[int(identities[s]-1)];
	    }
        }
    }
    for(int s=0; s<bead_number; ++s){
	probability_matrix_minus[s][s] = 0;
    }
}

void assign_probability_matrices(std::vector<int> &identities, int bead_number, double ranr, std::vector<double> proportions, std::vector<std::vector<double>>  &probability_matrix_plus, std::vector<std::vector<double>>  &probability_matrix_minus, std::vector<double> p_plus_probs, std::vector<double> p_minus_probs){
    std::vector<double> proportion_sums; 
    proportion_sums.push_back(0.0);
    double prop_ctr = 0;
    
    std::cout << "IDENTITIES: " << std::endl;
    for(int i=0; i<bead_number; ++i){
        std::cout << " particle: " << i << " : " <<  identities[i] << std::endl;
    }
    std::cout << "END IDENTITIES: " << std::endl;
    
    ////////////////////////////////
    // FORWARD PROBABILITY MATRIX //
    ////////////////////////////////

    int number_of_types = proportions.size();
    for(int s=0; s<bead_number; ++s){
        for(int t=s+1; t<bead_number; ++t){
            if (identities[s] == identities[t] and s != t){
	        probability_matrix_plus[s][t] = p_plus_probs[int(identities[s]-1)];
	        probability_matrix_plus[t][s] = p_plus_probs[int(identities[s]-1)];
	    }
	    else{
	        probability_matrix_plus[s][t] = p_plus_probs[p_plus_probs.size()-1];
	        probability_matrix_plus[t][s] = p_plus_probs[p_plus_probs.size()-1];
	    }
        }
    }
    for(int s=0; s<bead_number; ++s){
	probability_matrix_plus[s][s] = 0;
    }
    /////////////////////////////////
    // BACKWARD PROBABILITY MATRIX //
    /////////////////////////////////

    for(int s=0; s<bead_number; ++s){
        for(int t=s+1; t<bead_number; ++t){
            if (identities[s] == identities[t] and s != t){
	        probability_matrix_minus[s][t] = p_minus_probs[int(identities[s]-1)];
	        probability_matrix_minus[t][s] = p_minus_probs[int(identities[s]-1)];
	    }
	    else{
	        probability_matrix_minus[s][t] = p_minus_probs[p_minus_probs.size()-1];
	        probability_matrix_minus[t][s] = p_minus_probs[p_minus_probs.size()-1];
	    }
        }
    }
    for(int s=0; s<bead_number; ++s){
	probability_matrix_minus[s][s] = 0;
    }

}

////////////////////////////////////////////
// MODIFY PROB MATRICES AFTER POLY/DEPOLY //
////////////////////////////////////////////
void append_probabilities_to_matrix(std::vector<std::vector<double>> &probability_matrix_plus, std::vector<std::vector<double>> &probability_matrix_minus, std::vector<int> &identities, int bead_number, std::vector<double> p_plus_probs, std::vector<double> p_minus_probs){
    for (int k = 0 ; k < bead_number; ++k){
        for (int l = 0 ; l < bead_number; ++l){
         ////////////////////////////////////////////////
         // square root option for crossover particles //
         ////////////////////////////////////////////////
        probability_matrix_plus[l][k]  = sqrt(p_plus_probs[int(identities[k]-1)] * p_plus_probs[int(identities[l]-1)]);
        probability_matrix_plus[k][l]  = sqrt(p_plus_probs[int(identities[k]-1)] * p_plus_probs[int(identities[l]-1)]);
        probability_matrix_minus[l][k] = sqrt(p_minus_probs[int(identities[k]-1)] * p_minus_probs[int(identities[l]-1)]);
        probability_matrix_minus[k][l] = sqrt(p_minus_probs[int(identities[k]-1)] * p_minus_probs[int(identities[l]-1)]);
        }
    }
    for (int k = 0 ; k < bead_number; ++k){
        for (int l = 0 ; l < bead_number; ++l){
            if(k==l){
                probability_matrix_plus[l][k] = 0;
                probability_matrix_plus[k][l] = 0;
                probability_matrix_minus[l][k] = 0;
                probability_matrix_minus[k][l] = 0;
            }
        }
    }
}


