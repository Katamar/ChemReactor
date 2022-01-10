void reverse_function(std::string &test){

    int number;
    std::stringstream rev;
    std::vector<int> myNumbers;
    myNumbers.clear();
    std::stringstream iss( test );
    while ( iss >> number )
        myNumbers.push_back( number );
    std::reverse(std::begin(myNumbers), std::end(myNumbers));
    for(int d=0; d<myNumbers.size(); ++d){
	rev << myNumbers[d];
	rev << " ";
    }
    test = rev.str();
}

