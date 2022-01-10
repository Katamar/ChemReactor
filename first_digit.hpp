bool first_digit(std::string &str1, std::string &str2){

    int n1=0;
    std::vector<int> myNumbers1;
    myNumbers1.clear();
    std::stringstream s1( str1 );
    
    while ( s1 >> n1 )
        myNumbers1.push_back( n1 );
    
    int n2=0;
    std::vector<int> myNumbers2;
    myNumbers2.clear();
    std::stringstream s2( str2 );
    
    while ( s2 >> n2 )
        myNumbers2.push_back( n2 );
    return myNumbers1[0] == myNumbers2[0]; 
}

