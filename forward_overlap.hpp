void forward_overlap(std::vector<std::string> &chain_atom_ID_list, std::vector<std::string> &to_delete){

size_t pos;
std::string A, B, C;
std::string s1;
std::string s2;
bool fd;

for(int k=0; k<chain_atom_ID_list.size(); ++k){
    int ctr = 0;
    for(int l=k+1; l<chain_atom_ID_list.size(); ++l){
        for(int d=0; d<2; ++d){
	    if(ctr==0){
	    s1 = chain_atom_ID_list[k];
	    s2 = chain_atom_ID_list[l];
            }
	    if(ctr==1){
	    s1 = chain_atom_ID_list[k];
	    s2 = chain_atom_ID_list[l];
            reverse_function(s2);
            }
	    fd = first_digit(s1, s2);
	    if(s1.length()>s2.length()){
                pos = s1.find_last_of(s2[0]);
	    }
	    else{
                pos = s2.find_last_of(s1[0]);
	    }
	    pos = find_position(s1, s2);
	    if(pos == -1){
	        break;
	    }
	    else if(s1==s2){
	        to_delete.push_back(s1);
	    }
            else {
                if (fd){
            //    std::cout << "case 1: strings start with the same digit! " << std::endl;
                   if (s1.length() > s2.length()){
                       C = s1;
	               to_delete.push_back(s2);
                   }
	           else if(s2.length() > s1.length()){
                       C = s2;
	               to_delete.push_back(s1);
                   }
	        }

	        else if (not fd){
		    // case 2: two string don't start with the same digit
	            if(s1.length() > s2.length()){
                        A = s1;
                        B = s2;
                    }
                    else{
                        A = s2;
                        B = s1;
                    }
                    pos = A.find_last_of(B[0]);
	            if(pos + B.size() < A.size()){
                    // case 2: smchain_atom_ID_list string is in the middle e.g. 1 4 and 9 6 1 4 3 8 
	                C = A;
	                if (s1.length()>s2.length()){
	                    to_delete.push_back(s2);
	                    }
	                else{
	                    to_delete.push_back(s1);
	                    }
                    }
	            else{
                    // case 3: smchain_atom_ID_list string in the end e.g. 1 4 and 9 6 1 4  
	    	        C = A;
	                if (s1.length()>s2.length()){
	                    to_delete.push_back(s2);
	                    }
	                else{
	                    to_delete.push_back(s1);
	                    }
	            }
	        } // strings that begin with different digits
	    } // there is some overlap between strings
	    ctr +=1 ;
        } // reverse string loop
    }  // k loop
}  // l loop

for(int m=0; m<to_delete.size(); ++m){
    chain_atom_ID_list.erase(std::remove(chain_atom_ID_list.begin(), chain_atom_ID_list.end(), to_delete[m]), chain_atom_ID_list.end());
}

} // end of function forward_overlap


