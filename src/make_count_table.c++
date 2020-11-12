// Copyright 2020 Nadia Davidson. This program is distributed under the GNU
// General Public License. We also ask that you cite this software in
// publications where you made use of it for any part of the data
// analysis.

/** 
 ** Take a table of trsnscript counts and convert to gene-level counts.
 ** Assumes reads are uniquely mapped to transcripts.
 **
 ** Author: Nadia Davidson
 ** Modified: 2020
 **/ 

#include <iostream>
#include <fstream>
#include <istream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <algorithm>
#include <regex>
#include <stdlib.h>
#include <algorithm>    

using namespace std;

// the help information which is printed when a user puts in the wrong
// combination of command line options.
void print_usage(){
  cerr << endl;
  cerr << "Usage: make_count_table <transcript count table> <ref_table> > <out table>" << endl;
  cerr << endl;
}

// Main does the I/O and call multi_gene for each read and its associated
// set of transcript alignments.
int main(int argc, char **argv){

  //wrong number of arguements. Print help.
  if(argc!=3){
    print_usage();
    exit(1);
  }

  //Get the gene names and positions.
  //This is taken from process_transcriptome_alignment_table
  //Should probably separate this out as a function to remove redundancy
  /** 
   ** Now read in the gene to transcript ID mapping 
   **/
  ifstream file; 
  file.open(argv[2]);
  if(!(file.good())){
    cerr << "Unable to open file " << argv[3] << endl;
    exit(1);
  }
  //assume the first line is the header
  //use this to work out which columns to get
  map<string,string > trans_gene_map;
  string line;
  getline (file,line);
  int n_name=0; int n_name2=0; //column index
  istringstream line_stream(line);
  for(int i=0; line_stream ; i++){
    string temp;
    line_stream >> temp;
    if(temp=="name") n_name=i;
    if(temp=="name2") n_name2=i;
  }
  //loop over transcripts and fill information into data structure.
  while ( getline (file,line) ){ 
    istringstream line_stream(line);
    string gene; string trans; 
    for(int i=0; line_stream ; i++){
      string temp;
      line_stream >> temp;
      if(i==n_name) trans=temp; 
      if(i==n_name2) gene=temp; 
    }
    trans_gene_map[trans]=gene;
  }
  file.close();
  cerr << "Done reading in transcript IDs" << endl;

  //Now read the countTable. Store everything as a map
  map<string,int> gene_counts_map;
  
  string filename=argv[1];
  file.open(filename);
  if(!(file.good())){
    cerr << "Unable to open file " << filename << endl;
    exit(1);
  } 
  /**********  Read all the alignments ****************/
  cerr << "Reading the input transcript counts file, "<< filename << endl;
  while(getline(file,line) ){
    istringstream line_stream(line);
    int counts;
    string transcript;
    line_stream >> counts;
    line_stream >> transcript;
      
    smatch m; //extract the gene id
    regex_search(transcript,m,regex("_(EN[^_]*)__")) ; //Assumed ENSEMBL annotation naming here
    string gene=trans_gene_map[m[1].str()];
    gene_counts_map[gene]+=counts;
  }
  file.close();

  //loop over count map and output
  map<string, int>::iterator it = gene_counts_map.begin();
  while(it!=gene_counts_map.end()){
    cout << it->first << "\t" << it->second << endl;
    it++;
  }

  cerr << "Finished." << endl;

  return(0);
}
