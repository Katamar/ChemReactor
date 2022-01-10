void find_end(std::vector<int> &sequence, std::vector<std::vector<int>> &matrix, std::vector<std::vector<int>> &mask_end, std::vector<int> &bond_num, int begin_index, int bead_number, bool in_vector, int &next_index){   
      sequence.clear();
      sequence.push_back(begin_index);
      for(int j=0; j<bead_number; j++){
          if(matrix[begin_index][j]==1){
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
      if (sequence.size()>2){
          mask_end[sequence[0]][sequence[sequence.size()-1]] = 0;	    
          mask_end[sequence[sequence.size()-1]][sequence[0]] = 0;
      }	    

}

void find_middle_end(std::vector<int> &sequence, std::vector<std::vector<int>> &matrix, std::vector<std::vector<int>> &mask_end, std::vector<int> &bond_num, int begin_index, int bead_number, bool in_vector, int next_index){  
      sequence.clear();
      sequence.push_back(begin_index);
      for(int j=0; j<bead_number; j++){
          if(matrix[begin_index][j]==1){
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
      begin_index = sequence[sequence.size()-1];

      sequence.clear();
      sequence.push_back(begin_index);
      for(int j=0; j<bead_number; j++){
          if(matrix[begin_index][j]==1){
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

      if (sequence.size()>2){
          mask_end[sequence[0]][sequence[sequence.size()-1]] = 0;	    
          mask_end[sequence[sequence.size()-1]][sequence[0]] = 0;
      }	    
}

void delete_end(std::vector<int> &sequence, std::vector<std::vector<int>> &matrix, std::vector<std::vector<int>> &mask_end, std::vector<int> &bond_num, int begin_index, int bead_number, bool in_vector, int next_index){   
      sequence.clear();
      sequence.push_back(begin_index);
      for(int j=0; j<bead_number; j++){
          if(matrix[begin_index][j]==1){
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
      if (sequence.size()>2){
          mask_end[sequence[0]][sequence[sequence.size()-1]] = 1;	    
          mask_end[sequence[sequence.size()-1]][sequence[0]] = 1;
      }	    
}

void delete_middle_end(std::vector<int> &sequence, std::vector<std::vector<int>> &matrix, std::vector<std::vector<int>> &mask_end, std::vector<int> &bond_num, int begin_index, int bead_number, bool in_vector, int next_index){  
      sequence.clear();
      sequence.push_back(begin_index);
      for(int j=0; j<bead_number; j++){
          if(matrix[begin_index][j]==1){
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
      begin_index = sequence[sequence.size()-1];

      sequence.clear();
      sequence.push_back(begin_index);
      for(int j=0; j<bead_number; j++){
          if(matrix[begin_index][j]==1){
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

      if (sequence.size()>2){
          mask_end[sequence[0]][sequence[sequence.size()-1]] = 1;	    
          mask_end[sequence[sequence.size()-1]][sequence[0]] = 1;
      }
}
