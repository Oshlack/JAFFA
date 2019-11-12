// Copyright 2019 Nadia Davidson for Murdoch Childrens Research
// Institute Australia. This program is distributed under the GNU
// General Public License. We also ask that you cite this software in
// publications where you made use of it for any part of the data
// analysis.

/** make_simple_read_table is a simple program to count the number of
 ** discordant pairs that support a fusions.
 **
 **
 ** Author: Nadia Davidson
 ** Modified: November 2019
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

using namespace std;

// the help information which is printed when a user puts in the wrong
// combination of command line options.
void print_usage(){
  cout << endl;
  cout << "Usage: cat <bam file> | make_simple_read_table <candidate_fusion_table> <trans_gene_map> >  <out table>" << endl;
  cout << endl;
}

struct fusion_candidate {
  string read ;
  int break_min ;
  int break_max ;
  pair<string,string> fusion ;
};

// the real stuff starts here.
int main(int argc, char **argv){

  if(argc!=3){
    print_usage();
    exit(1);
  }

  /**
   ** Read the list of candidate fusions
   **/
  ifstream file;
  file.open(argv[1]);
  if(!(file.good())){
    cerr << "Unable to open file " << argv[1] << endl;
    exit(1);
  }
  string line;
  vector<fusion_candidate>  candidate_table;
  vector< pair<string,string> > fusion_list;
  // read in the table of candidate fusions
  while ( getline (file,line) ){ 
    fusion_candidate new_cand;
    istringstream line_stream(line);
    line_stream >> new_cand.read;
    line_stream >> new_cand.break_min;
    line_stream >> new_cand.break_max;
    string temp;
    line_stream >> temp;
    smatch m; //extract the gene id and sort alphabetically
    regex_search(temp,m,regex("(.*):(.*)"));
    if(m[1].str()<m[2].str())
      new_cand.fusion=make_pair(m[1].str(), m[2].str());
    else
      new_cand.fusion=make_pair(m[2].str(), m[1].str());
    //sort
    fusion_list.push_back(new_cand.fusion);
    candidate_table.push_back(new_cand);
  }
  file.close();
  // remove duplicates in the fusion list
  sort( fusion_list.begin(), fusion_list.end() );
  fusion_list.erase( unique( fusion_list.begin(), fusion_list.end() ), fusion_list.end() );
  cerr << "Done reading in candidate fusions" << endl;

  /** 
   ** Now read in the gene to transcript ID mapping 
   **/
  file.open(argv[2]);
  if(!(file.good())){
    cerr << "Unable to open file " << argv[2] << endl;
    exit(1);
  }
  //assume the first line is the header
  getline (file,line);
  int n_name=0; int n_name2=0;
  istringstream line_stream(line);
  for(int i=0; line_stream ; i++){
    string temp;
    line_stream >> temp;
    if(temp=="name") n_name=i;
    if(temp=="name2") n_name2=i;
  }
  map<string, vector<string> > gene_trans_map ;
  while ( getline (file,line) ){
    istringstream line_stream(line);
    string gene;
    string trans;
    for(int i=0; line_stream ; i++){
      string temp;
      line_stream >> temp;
      if(i==n_name) trans=temp;
      if(i==n_name2) gene=temp;
    }
    gene_trans_map[gene].push_back(trans);
  }
  file.close();
  cerr << "Done reading in transcript IDs" << endl;

  /** 
   ** Now read in the mapped reads
   **/
  map<string,vector<string> > trans_read_map;
  vector<string> trans_names;
  while ( getline (cin,line) ){
    string temp1;
    string temp2;
    istringstream line_stream(line);
    line_stream >> temp1;
    line_stream >> temp2; line_stream >> temp2;
    trans_read_map[temp2].push_back(temp1);
  }
  // Loop over all the trans_read_map and fix the names:
  map<string,vector<string> > trans_read_map_fixed;
  map<string,vector<string> >::iterator tr_itr=trans_read_map.begin();
  for(;tr_itr!=trans_read_map.end() ; tr_itr++){
    smatch m;
    regex_search(tr_itr->first,m,regex("ENST\\d{11}(.\\d+)?"));
    trans_read_map_fixed[m[0].str()]=tr_itr->second;
  }
  trans_read_map.clear();
  cerr << "Done reading in bam file" << endl;

  /**
   ** now loop over each gene and get all the reads aligning to it
   **/
  map<string, pair<vector < string> , vector< string > > > gene_reads;
  map<string,vector<string> >::iterator gt_itr=gene_trans_map.begin();
  //loop over genes and get transcripts id
  for(;gt_itr!=gene_trans_map.end(); gt_itr++){
    vector<string> reads;
    //for all transcripts look up the reads
    for(int t=0; t < gt_itr->second.size(); t++){
       vector<string> reads_to_add = trans_read_map_fixed[gt_itr->second.at(t)];
      reads.insert(reads.end(), reads_to_add.begin(), reads_to_add.end());
    }
    // remove redundance
    sort( reads.begin(), reads.end() );
    reads.erase( unique( reads.begin(), reads.end() ), reads.end() );   
    //now loop over the reads and separate the read start and ends
    for(int r=0; r<reads.size(); r++){
      smatch m;
      regex_search(reads.at(r),m,regex("(.*)/([12])$"));
      if(m[2].str()=="1")
	gene_reads[gt_itr->first].first.push_back(m[1].str());
      else
	gene_reads[gt_itr->first].second.push_back(m[1].str());
    }
  }
  trans_read_map_fixed.clear();
  gene_trans_map.clear();
  cerr << "Done assigning reads to genes" << endl;

  /**
   ** Calculate the counts for each fusion
   **/
  //loop over the fusion list
  map<pair<string,string >, int> spanning_reads;
  for(int f=0; f<fusion_list.size(); f++){
    //check from intersection of read ids.
    vector<string> g1_r1=gene_reads[fusion_list.at(f).first].first;
    vector<string> g1_r2=gene_reads[fusion_list.at(f).first].second;
    vector<string> g2_r1=gene_reads[fusion_list.at(f).second].first;
    vector<string> g2_r2=gene_reads[fusion_list.at(f).second].second;
    //    cout << fusion_list.at(f) << endl;
    int total=0;
    for(vector<string>::iterator i = g1_r1.begin(); i!=g1_r1.end(); ++i){
      if (find(g2_r2.begin(), g2_r2.end(), *i) != g2_r2.end()){
	//cout << *i << endl;
	total++;
      }
    }
    for(vector<string>::iterator i = g2_r1.begin(); i!=g2_r1.end(); ++i){
      if (find(g1_r2.begin(), g1_r2.end(), *i) != g1_r2.end()){
	//cout << *i << endl;
	total++;
      }
    }
    spanning_reads[fusion_list.at(f)]=total;
    //cout << total << endl;
  }
  cerr << "Done calculating spanning pairs" << endl;

  cout << "transcript\tbreak_min\tbreak_max\tfusion_genes\tspanning_pairs\tspanning_reads" << endl;
  vector<fusion_candidate>::iterator cand_itr=candidate_table.begin();
  for(;cand_itr!=candidate_table.end(); cand_itr++){
    cout << cand_itr->read << "\t" ;
    cout << cand_itr->break_min << "\t" ;
    cout << cand_itr->break_max << "\t" ;
    cout << cand_itr->fusion.first << ":" << cand_itr->fusion.second << "\t" ;
    cout << spanning_reads[cand_itr->fusion] << "\t" ;
    cout << "1" << endl;
  }

  return(0);
}
