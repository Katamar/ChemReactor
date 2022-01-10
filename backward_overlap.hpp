void vector_to_string(std::vector<int> &vector, std::string &string){
    std::stringstream result;
    std::copy(std::begin(vector), std::end(vector), std::ostream_iterator<int>(result, " "));
    string = result.str();
}

void backward_overlap(int s, int t, std::vector<std::string> &chain_atom_ID_list){
    int index1;
    int index2;
    int n1;
    std::vector<int> myNumbers1;
    std::vector<std::string> tmp_chain_atom_ID_list;
    std::vector<int> part1;
    std::vector<int> part2;
    std::string l1;
    std::string l2;
    std::string l3;

    for(int i=0; i<chain_atom_ID_list.size(); ++i){
	if(chain_atom_ID_list[i].size()==1){
            tmp_chain_atom_ID_list.push_back(chain_atom_ID_list[i]);
            continue;
	}

	myNumbers1.clear();
	part1.clear();
	part2.clear();
        myNumbers1.clear();
        std::stringstream s1( chain_atom_ID_list[i] );
        while ( s1 >> n1 )
            myNumbers1.push_back( n1 );
        int both = 0;
        int in_string = 0;
        for(int k=0; k<myNumbers1.size(); ++k){
		if (myNumbers1[k]==s){
                    index1 = k;
		    in_string = 1;
		    both += 1;
		}
		if (myNumbers1[k]==t){
                    index2 = k;
		    in_string = 1;
		    both += 1;
		}
		if ( both==2 ){
                    break;
		}
	}
	if(in_string == 0) {
            tmp_chain_atom_ID_list.push_back(chain_atom_ID_list[i]);
            continue;
	}
	// s > t
	if (index1 > index2){
	    // part1: o -- t, part2: s -- o 
            for(int k=0; k<index1; ++k){
                part1.push_back(myNumbers1[k]);
	    }
            for(int k=index1; k<myNumbers1.size(); ++k){
                part2.push_back(myNumbers1[k]); 
	    }
	}
	// t > s
	else{
	    // part1: o -- s , part2: t -- o
            for(int k=0; k<index2; ++k){
                part1.push_back(myNumbers1[k]); 
	    }
            for(int k=index2; k<myNumbers1.size(); ++k){
                part2.push_back(myNumbers1[k]); 
	    }
	}

        vector_to_string(part1, l1); 
        vector_to_string(part2, l2); 
    }
    tmp_chain_atom_ID_list.push_back(l1);
    tmp_chain_atom_ID_list.push_back(l2);
    chain_atom_ID_list.clear();
    chain_atom_ID_list = tmp_chain_atom_ID_list; 
}

void backward_overlap_delete(int s, int t, std::vector<std::string> &chain_atom_ID_list){
    int index1;
    int index2;
    int n1;
    std::vector<int> myNumbers1;
    std::vector<std::string> tmp_chain_atom_ID_list;
    std::vector<int> part1;
    std::vector<int> part2;
    std::string l1;
    std::string l2;
    std::string l3;

    for(int i=0; i<chain_atom_ID_list.size(); ++i){
	if(chain_atom_ID_list[i].size()==1){
            tmp_chain_atom_ID_list.push_back(chain_atom_ID_list[i]);
            continue;
	}

	myNumbers1.clear();
	part1.clear();
	part2.clear();
        myNumbers1.clear();
        std::stringstream s1( chain_atom_ID_list[i] );
        while ( s1 >> n1 )
            myNumbers1.push_back( n1 );
        int both = 0;
        int in_string = 0;
        for(int k=0; k<myNumbers1.size(); ++k){
		if (myNumbers1[k]==s){
                    index1 = k;
		    in_string = 1;
		    both += 1;
		}
		if (myNumbers1[k]==t){
                    index2 = k;
		    in_string = 1;
		    both += 1;
		}
		if ( both==2 ){
                    break;
		}
	}
	if(in_string == 0) {
            tmp_chain_atom_ID_list.push_back(chain_atom_ID_list[i]);
            continue;
	}
	// s > t
	if (index1 > index2){
	    // part1: o -- t, part2: s -- o 
            for(int k=0; k<index1; ++k){
                part1.push_back(myNumbers1[k]);
	    }
            for(int k=index1; k<myNumbers1.size(); ++k){
                part2.push_back(myNumbers1[k]); 
	    }
	}
	// t > s
	else{
	    // part1: o -- s , part2: t -- o
            for(int k=0; k<index2; ++k){
                part1.push_back(myNumbers1[k]); 
	    }
            for(int k=index2; k<myNumbers1.size(); ++k){
                part2.push_back(myNumbers1[k]); 
	    }
	}

        vector_to_string(part1, l1); 
        vector_to_string(part2, l2); 
    }
    chain_atom_ID_list.clear();
    chain_atom_ID_list = tmp_chain_atom_ID_list; 
}

void chain_growth_backward_overlap_delete(int s, int t, std::vector<std::string> &chain_atom_ID_list){
    int index1;
    int index2;
    int n1;
    std::vector<int> myNumbers1;
    std::vector<std::string> tmp_chain_atom_ID_list;
    std::vector<int> part1;
    std::vector<int> part2;
    std::string l1;
    std::string l2;
    std::string l3;

    for(int i=0; i<chain_atom_ID_list.size(); ++i){
	if(chain_atom_ID_list[i].size()==1){
            tmp_chain_atom_ID_list.push_back(chain_atom_ID_list[i]);
            continue;
	}

	myNumbers1.clear();
	part1.clear();
	part2.clear();
        myNumbers1.clear();
        std::stringstream s1( chain_atom_ID_list[i] );
        while ( s1 >> n1 )
            myNumbers1.push_back( n1 );
        int both = 0;
        int in_string = 0;
        for(int k=0; k<myNumbers1.size(); ++k){
		if (myNumbers1[k]==s){
                    index1 = k;
		    in_string = 1;
		    both += 1;
		}
		if (myNumbers1[k]==t){
                    index2 = k;
		    in_string = 1;
		    both += 1;
		}
		if ( both==2 ){
                    break;
		}
	}
	if(in_string == 0) {
            tmp_chain_atom_ID_list.push_back(chain_atom_ID_list[i]);
            continue;
	}
	// s > t
	if (index1 > index2){
	    // part1: o -- t, part2: s -- o 
            for(int k=0; k<index1; ++k){
                part1.push_back(myNumbers1[k]);
	    }
            for(int k=index1; k<myNumbers1.size(); ++k){
                part2.push_back(myNumbers1[k]); 
	    }
	}
	// t > s
	else{
	    // part1: o -- s , part2: t -- o
            for(int k=0; k<index2; ++k){
                part1.push_back(myNumbers1[k]); 
	    }
            for(int k=index2; k<myNumbers1.size(); ++k){
                part2.push_back(myNumbers1[k]); 
	    }
	}

        vector_to_string(part1, l1); 
        vector_to_string(part2, l2); 
    }
    if( l1.size() > l2.size()){
      tmp_chain_atom_ID_list.push_back(l1);
    }
    if( l2.size() > l1.size()){
      tmp_chain_atom_ID_list.push_back(l2);
    }
    chain_atom_ID_list.clear();
    chain_atom_ID_list = tmp_chain_atom_ID_list; 
}

void step_growth_backward_overlap_delete(int s, int t, std::vector<std::string> &chain_atom_ID_list){
    int index1;
    int index2;
    int n1;
    std::vector<int> myNumbers1;
    std::vector<std::string> tmp_chain_atom_ID_list;
    std::vector<int> part1;
    std::vector<int> part2;
    std::string l1;
    std::string l2;
    std::string l3;

    for(int i=0; i<chain_atom_ID_list.size(); ++i){
	if(chain_atom_ID_list[i].size()==1){
            tmp_chain_atom_ID_list.push_back(chain_atom_ID_list[i]);
            continue;
	}

	myNumbers1.clear();
	part1.clear();
	part2.clear();
        myNumbers1.clear();
        std::stringstream s1( chain_atom_ID_list[i] );
        while ( s1 >> n1 )
            myNumbers1.push_back( n1 );
        int both = 0;
        int in_string = 0;
        for(int k=0; k<myNumbers1.size(); ++k){
		if (myNumbers1[k]==s){
                    index1 = k;
		    in_string = 1;
		    both += 1;
		}
		if (myNumbers1[k]==t){
                    index2 = k;
		    in_string = 1;
		    both += 1;
		}
		if ( both==2 ){
                    break;
		}
	}
	if(in_string == 0) {
            tmp_chain_atom_ID_list.push_back(chain_atom_ID_list[i]);
            continue;
	}
	// s > t
	if (index1 > index2){
	    // part1: o -- t, part2: s -- o 
            for(int k=0; k<index1; ++k){
                part1.push_back(myNumbers1[k]);
	    }
            for(int k=index1; k<myNumbers1.size(); ++k){
                part2.push_back(myNumbers1[k]); 
	    }
	}
	// t > s
	else{
	    // part1: o -- s , part2: t -- o
            for(int k=0; k<index2; ++k){
                part1.push_back(myNumbers1[k]); 
	    }
            for(int k=index2; k<myNumbers1.size(); ++k){
                part2.push_back(myNumbers1[k]); 
	    }
	}

        vector_to_string(part1, l1); 
        vector_to_string(part2, l2); 
    }
    
    std::istringstream is1( l1 );
    std::istringstream is2( l2 );
    int n11;
    int n22;
    int ct1 = 0;
    int ct2 = 0;
    while( is1 >> n11 ) {
        ct1 += 1;
    }
    while( is2 >> n22 ) {
        ct2 += 1;
    }
    

    if (ct1 > 1){
    tmp_chain_atom_ID_list.push_back(l1);
    }
    if (ct2 > 1){
    tmp_chain_atom_ID_list.push_back(l2);
    }
    chain_atom_ID_list.clear();
    chain_atom_ID_list = tmp_chain_atom_ID_list; 
}
