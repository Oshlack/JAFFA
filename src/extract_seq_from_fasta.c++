// Copyright 2019 Nadia Davidson for Murdoch Childrens Research
// Institute Australia. This program is distributed under the GNU
// General Public License. We also ask that you cite this software in
// publications where you made use of it for any part of the data
// analysis.

/** extract_seq_from_fasta is a simple program to extract sequences of interest
 ** from a fasta file, outputing a new fasta. Essentially it's a fast grep.
 **
 ** Author: Nadia Davidson, nadia.davidson@mcri.edu.au
 ** Modified: 4 August 2014
 **/ 

#include <iostream>
#include <fstream>
#include <istream>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <stdlib.h>

using namespace std;

// the help information which is printed when a user puts in the wrong
// combination of command line options.
void print_usage(){
  cout << endl;
  cout << "extract_seq_from_fasta extracts the sequences from a fasta file using a list of sequence IDs" << endl;
  cout << endl;
  cout << "Usage: cat <fasta file> | extract_seq_from_fasta <id file>  >  <out file>" << endl;
  cout << "         <id>           - A list of sequence ids" << endl;
  cout << "         <fasta file>   - Input fasta" << endl; 
  cout << "         <out file>     - Filtered output fasta" << endl;
  cout << endl;
}

// the real stuff starts here.
int main(int argc, char **argv){

  vector<string> ids_to_keep;
  if(argc!=2){
    print_usage();
    exit(1);
  }

  ifstream file;
  //Open the id file and put the list into a vector
  file.open(argv[1]);
  if(!(file.good())){
    cout << "Unable to open file " << argv[1] << endl;
    exit(1);
  }
  string line;
  string id;
  while ( getline (file,line) ){
    istringstream line_stream(line);
    line_stream >> id;
    ids_to_keep.push_back(id);
  }
  file.close();

  //Now open the fasta file and report if the id is in the list.
  /**  file.open(argv[2]);
  if(!(file.good())){
    cout << "Unable to open file " << argv[1] << endl;
    exit(1);
    }**/
  bool print=false;
  while ( getline (cin,line) ){
    int start=line.find(">")+1;
    if(start==1){ //if this is the ID line...
      int end=line.find_first_of("\t\n ")-1;
      string id=line.substr(start,end);
      if(find(ids_to_keep.begin(), ids_to_keep.end(), id) != ids_to_keep.end()){
	print=true;
	cout << line << endl;
      } else {
	print=false;
      }
    } else {
      if(print)
	cout << line << endl; //output
    }
  }
  //  file.close();
  
  return(0);
}
