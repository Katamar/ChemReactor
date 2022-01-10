void write_ends(std::vector<int> &sequence, std::vector<int> &polymer_ends, std::vector<std::vector<int>> &matrix, std::vector<int> &bond_num, int bead_number, int sum_row, int chain_number, int next_index, bool in_vector, std::ofstream &myfile4, std::ofstream &myfile5, int k){	
    polymer_ends.clear();
    for(int i=0; i<bead_number; i++){
        sum_row = 0;
        for(int j=0; j<bead_number; j++){
            sum_row += matrix[i][j];
        }
	if (sum_row == 1) {
	    polymer_ends.push_back(i);
	}
    }
    chain_number = (polymer_ends.size())/2;
    myfile4 << chain_number << std::endl;
    
    for(int s=0; s<polymer_ends.size(); s++){
        sequence.clear();
	sequence.push_back(polymer_ends[s]);
        for(int j=0; j<bead_number; j++){
            if(matrix[polymer_ends[s]][j]==1){
	        sequence.push_back(j);
	        next_index = j;
	    }
        }
	while(bond_num[next_index]==2){
            for(int j=0; j<bead_number; j++){
		in_vector = std::find(sequence.begin(), sequence.end(), j) != sequence.end();
                if(matrix[next_index][j]==1 and in_vector != 1){
	    	    sequence.push_back(j);
		    next_index = j;
	            continue;
	        }
	    }
	}
        myfile5 << sequence.size() << std::endl;
        polymer_ends.erase(std::remove(polymer_ends.begin(), polymer_ends.end(), sequence[0]), polymer_ends.end());	
        polymer_ends.erase(std::remove(polymer_ends.begin(), polymer_ends.end(), sequence[sequence.size()-1]), polymer_ends.end());
        s=-1;
    }
}
