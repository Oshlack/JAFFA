// Copyright 2021 Nadia Davidson. This program is distributed under the GNU
// General Public License. We also ask that you cite this software in
// publications where you made use of it for any part of the data
// analysis.

/** make_3_gene_fusion_table checks for reads with multiple fusions
 **
 ** Author: Nadia Davidson
 ** Modified: July 2021
 **/

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <set>
#include <vector>
#include <unordered_map>
#include <stdlib.h>

using namespace std;

//struct for fusion information
struct fusion {
  string fusion_name ; 
  string classification ; //highest..
  int read_support ; //aggregate over breakpoints
};

int main(int argc, char **argv){

  if(argc!=4){
    cout << "Usage: make_3_gene_fusion_table <.summary> <.txt> [output readID file]" << endl;
    exit(1);
  }

  /**
   ** Read the list of fusions
   **/
  ifstream file;
  file.open(argv[1]);
  if(!(file.good())){
    cerr << "Unable to open file " << argv[1] << endl;
  }
  string line;
  unordered_map< string, fusion> fusion_map;
  // read in the table of candidate fusions
  getline (file,line); //skip the header line
  while ( getline (file,line) ){
    istringstream line_stream(line);
    string temp;
    for(int i=0; i<2; i++) line_stream >> temp; //skip 2 columns
    string reads_string; 
    line_stream >> reads_string;
    int reads = atoi(reads_string.c_str());
    for(int i=3; i<14; i++) line_stream >> temp; //skip 11 columns
    string fname;
    line_stream >> fname;
    for(int i=15; i<18; i++) line_stream >> temp; //skip 2 columns
    string classif;
    line_stream >> classif;
    if(fusion_map.find(fname)==fusion_map.end()){
      struct fusion newFusion; 
      newFusion.fusion_name=fname;
      newFusion.classification=classif;
      newFusion.read_support=reads;
      fusion_map[fname]=newFusion;
    } else {
      fusion_map[fname].read_support+=reads;
      if(classif=="HighConfidence" || 
	 (classif=="LowConfidence" && fusion_map[fname].classification!="HighConfidence"))
	  fusion_map[fname].classification=classif;
    }
  }
  file.close();

  //now read the .txt file
  file.open(argv[2]);
  if(!(file.good())){
    cerr << "Unable to open file " << argv[2] << endl;
  }

  unordered_map< string, set<string> > fusion_reads;
  // read in the table of fusion reads
  while ( getline (file,line) ){
    istringstream line_stream(line);
    string read_id;
    line_stream >> read_id;
    string temp;
    for(int i=0; i<2; i++) line_stream >> temp; //skip 2 columns
    string fname;
    line_stream >> fname;
    if(fusion_map.find(fname)!=fusion_map.end())
      fusion_reads[read_id].insert(fname);
  }

  //loop through the reads, only keeping ones with multiple fusions.
  ofstream readOutFile;
  readOutFile.open (argv[3]);
  unordered_map< string , set< string> >::iterator itr = fusion_reads.begin();
  unordered_map< string , set< string> >::iterator itrEnd = fusion_reads.end();
  unordered_map<string,int> tfusion_counts;
  while(itr!=itrEnd){
    set< string> fs = itr->second;
    if(fs.size()>1){
      set< string >::iterator itrF = fs.begin();
      set< string >::iterator itrFEnd = fs.end();
      vector<string> fgenes;
      while(itrF!=itrFEnd){ 
	string fus_name=*itrF; //split the fusion names into their genes
	string start_gene = fus_name.substr(0, fus_name.find(":"));
	string end_gene = fus_name.substr(fus_name.find(":")+1,fus_name.size());
	if(fgenes.size()==0){ // add first set of genes
	  fgenes.push_back(start_gene) ; 
	  fgenes.push_back(end_gene) ;
	} else {
	  if(fgenes[0]==end_gene){ //work out the position of this fusion relative to others
	    fgenes.insert(fgenes.begin(),start_gene); //shift all up..
	  } else if (fgenes[fgenes.size()-1]==start_gene) {
	    fgenes.push_back(end_gene);
	  } else {
	    cerr << "read id: " << itr->first << " does not conform to 3 gene fusion structure." << endl;
	    fgenes.clear();
	    break;
	  }
	}
	itrF++;
      }
      if(fgenes.size()>=3){ 
	string tfname=fgenes[0];
	for(int i=1; i<fgenes.size() ; i++)
	  tfname+=":"+fgenes[i];
	readOutFile << itr->first << "\t" << tfname << endl;
	tfusion_counts[tfname]++;
      }
    }
    itr++;
  }
  readOutFile.close();

  //now output the 3 fusion results
  unordered_map< string, int>::iterator itrTF=tfusion_counts.begin();
  //print the header:
  cout << "Fusion\tReads\tConstituentFusionReads\tConstituentFusionClassifications\tClassification" << endl;
  while(itrTF!=tfusion_counts.end()){
    string tfname=itrTF->first;
    cout << tfname << "\t" << itrTF->second << "\t";

    vector<string> gene;
    stringstream ss(tfname);
    string g; //tokensize the 3 gene fusion name
    while(getline(ss, g,':')){
      gene.push_back(g);
    } //output the number of reads for each constituent fusion
    cout << fusion_map[gene[0]+":"+gene[1]].read_support;
    for(int i=2; i < gene.size() ; i++){
      cout << ":" << fusion_map[gene[i-1]+":"+gene[i]].read_support;
    }
    cout << "\t"; 
    //output the classification for each constituent fusion
    string classification=fusion_map[gene[0]+":"+gene[1]].classification;
    cout << classification;
    for(int i=2; i < gene.size() ; i++){
      string temp_class=fusion_map[gene[i-1]+":"+gene[i]].classification;
      cout << ":" << temp_class;
      if(classification=="HighConfidence" ||
         (classification=="LowConfidence" && temp_class!="HighConfidence"))
	classification=temp_class;
    }
    cout << "\t" << classification << endl;  
    itrTF++;
  }

}
