int find_position(std::string &s1, std::string &s2){
    int pos = -1;
    int number_c1;
    int number_c2;
    std::vector<int> c1;
    std::vector<int> c2;
    c1.clear();
    c2.clear();
    std::stringstream iss_c1(s1);
    std::stringstream iss_c2(s2);
    while ( iss_c1 >> number_c1)
        c1.push_back( number_c1);
    while ( iss_c2 >> number_c2)
        c2.push_back( number_c2);
    for(int d=0; d<c1.size(); ++d){
	iss_c1 << c1[d];
    }
    for(int d=0; d<c2.size(); ++d){
	iss_c2 << c2[d];
    }
    for(int d=0; d<c1.size(); ++d){
        for(int s=0; s<c2.size(); ++s){
	    if(c1[d]==c2[s]){
		pos = d;
		break;
	    }
	}
	if (pos != -1){
	    break;
        }
    }
    return pos;
}
