// Copyright 2019 Nadia Davidson. This program is distributed under the GNU
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
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <regex>
#include <stdlib.h>

using namespace std;

// the help information which is printed when a user puts in the wrong
// combination of command line options.
void print_usage(){
  cerr << endl;
  cerr << "Usage: cat <bam file> | cut -f1-3 | make_simple_read_table <candidate_fusion_table> <trans_gene_map> >  <out table>" << endl;
  cerr << endl;
}

//struct for fusion information
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
  vector< string > genes_of_interest;
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
    genes_of_interest.push_back(m[1].str());
    genes_of_interest.push_back(m[2].str());
  }
  file.close();
  // remove duplicates in the fusion list
  sort( fusion_list.begin(), fusion_list.end() );
  fusion_list.erase( unique( fusion_list.begin(), fusion_list.end() ), fusion_list.end() );
  sort( genes_of_interest.begin(), genes_of_interest.end() );
  genes_of_interest.erase( unique( genes_of_interest.begin(), genes_of_interest.end() ), genes_of_interest.end() );
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
  unordered_map<string, vector<string> > gene_trans_map ;
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
    //only include genes in candidate fusions
    if(find(genes_of_interest.begin(), genes_of_interest.end(), gene) != genes_of_interest.end())
      gene_trans_map[gene].push_back(trans);
  }
  file.close();
  cerr << "Done reading in transcript IDs" << endl;

  /** 
   ** Now read in the mapped reads
   **/
  unordered_map< string, vector<string > > trans_read_map;
  vector<string> trans_names;
  int i=0;
  string temp1;
  string temp2;
  string temp3;
  while(!cin.eof()){
    cin >> temp1 >> temp2 >> temp3 ;
    trans_read_map[temp3].push_back(temp1);
    if(i%1000000==0)
      cerr << i << " alignments read" << endl;
    i++;
  }
  cerr << "Done reading in bam file" << endl;
  // Loop over all the trans_read_map and fix the names:
  unordered_map<string,vector<string> > trans_read_map_fixed;
  unordered_map<string,vector<string> >::iterator tr_itr=trans_read_map.begin();
  for(;tr_itr!=trans_read_map.end() ; tr_itr++){
    smatch m;
    regex_search(tr_itr->first,m,regex("ENST\\d{11}(.\\d+)?"));
    trans_read_map_fixed[m[0].str()]=tr_itr->second;
  }
  trans_read_map.clear();
  cerr << "Done getting trans ids" << endl;

  /**
   ** now loop over each gene and get all the reads aligning to it
   **/
  unordered_map<string, pair<vector < string> , vector< string > > > gene_reads;
  unordered_map<string,vector<string> >::iterator gt_itr=gene_trans_map.begin();
  //loop over genes and get transcripts id
  for(;gt_itr!=gene_trans_map.end(); gt_itr++){
    vector<string> reads;
    //for all transcripts look up the reads
    for(int t=0; t < gt_itr->second.size(); t++){
      vector<string> reads_to_add = trans_read_map_fixed[gt_itr->second.at(t)];
      reads.insert(reads.end(), reads_to_add.begin(), reads_to_add.end());
    }
    // remove redundancy
    sort( reads.begin(), reads.end() );
    reads.erase( unique( reads.begin(), reads.end() ), reads.end() );   
    //now loop over the reads and separate the read start and ends
    for(int r=0; r<reads.size(); r++){
      //separate the read id and pair end
      string read_id=reads.at(r).substr(0,reads.at(r).find("/"));//reads.at(r).size()-1); 
      char read_end=reads.at(r).back();
      if(read_end=='1') // First of pair (assumes the read IDs ends with 1).
	gene_reads[gt_itr->first].first.push_back(read_id);
      else
	gene_reads[gt_itr->first].second.push_back(read_id); 
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
    int total=0;
    unordered_set<string> temp_set1(g1_r1.begin(),g1_r1.end());
    temp_set1.insert(g2_r2.begin(), g2_r2.end());
    total+=g1_r1.size() + g2_r2.size() - temp_set1.size();

    unordered_set<string> temp_set2(g2_r1.begin(),g2_r1.end());
    temp_set2.insert(g1_r2.begin(), g1_r2.end());
    total+=g1_r2.size() + g2_r1.size() - temp_set2.size();
    spanning_reads[fusion_list.at(f)]=total;
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
