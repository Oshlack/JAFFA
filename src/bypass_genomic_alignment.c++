// Copyright 2019 Nadia Davidson for Peter Mac
// This program is distributed under the GNU
// General Public License. We also ask that you cite this software in
// publications where you made use of it for any part of the data
// analysis.

/** 
 **
 ** Author: Nadia Davidson
 ** Jan. 2020
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
  cerr << "Usage: bypass_genomic_alignment <ref table> <fusion candidates> > <out table>" << endl;
  cerr << endl;
}

//static const int MAX_OVERLAP=15 ;


class Transcript {
public:
  string chrom;
  string exon_starts;
  string exon_ends;
  string strand;
  vector<int> block_starts;
  vector<int> block_ends;

  void set_values(string chr, string es, string ee, string st){
    chrom=chr;
    exon_starts=es;
    exon_ends=ee;
    strand=st;
  }

  void set_starts_and_end(){
    stringstream sses(exon_starts);
    stringstream ssee(exon_ends);
    while( sses.good() && ssee.good()  ){
	string block_start_string;
	string block_end_string;
	getline(sses, block_start_string, ',' );
	getline(ssee, block_end_string, ',' );
	block_starts.push_back(atoi(block_start_string.c_str()));
	block_ends.push_back(atoi(block_end_string.c_str()));
    }
    //remove the last blocks as these are zero (string ends with comma)
    block_starts.pop_back();
    block_ends.pop_back();
  };

  // Get the genomic coordinates of position "pos" along the transcript.
  // first base is at position 1?
  void get_position(int pos, bool start, int &base, int &offset_from_exon){
    if(block_starts.size()==0) set_starts_and_end(); //convert block positions from string to int first time.

    int backbuff=10;
    if(start) pos=pos-backbuff;
    else pos=pos+backbuff;

    //switch position for negative strand transcripts..
    if(strand=="-"){
      int trans_length=0;
      for(int i=0; i< block_starts.size(); i++)
	trans_length+=(block_ends.at(i)-block_starts.at(i));
      pos=trans_length-pos;
    }

    int cumulative_bases=0;
    for(int i=0; i< block_starts.size(); i++){
	int block_size=block_ends.at(i)-block_starts.at(i);
	int current_offset=block_size+cumulative_bases-pos;
	if(current_offset>=0){
	  base=block_ends.at(i)-current_offset;
	  if(start==(strand=="+")){
	    offset_from_exon=current_offset-backbuff;
	    base+=backbuff;
	  }
	  else{
	    offset_from_exon=pos-cumulative_bases-backbuff;
	    base-=backbuff;
	  }
	  return;
	} else {
	  cumulative_bases+=block_size;
	}
      }
  }

};


int main(int argc, char **argv){

  if(argc!=3){
    print_usage();
    exit(1);
  }

  /** 
   ** Now read in the gene to transcript ID mapping 
   **/
  map<string, Transcript > trans_positions;
  ifstream file; 
  file.open(argv[1]);
  if(!(file.good())){
    cerr << "Unable to open file " << argv[3] << endl;
    exit(1);
  }
  string line;
  getline (file,line);
  int n_name=0; int n_strand=0; int n_chrom=0; int n_exonStarts=0; int n_exonEnds=0; 
  istringstream line_stream(line);
  for(int i=0; line_stream ; i++){
    string temp;
    line_stream >> temp;
    if(temp=="name") n_name=i;
    if(temp=="strand") n_strand=i;
    if(temp=="chrom") n_chrom=i;
    if(temp=="exonStarts") n_exonStarts=i;
    if(temp=="exonEnds") n_exonEnds=i;
  }
  while ( getline (file,line) ){
    istringstream line_stream(line);
    string gene; string trans; string chrom;
    string starts; string ends; string strand;
    for(int i=0; line_stream ; i++){
      string temp;
      line_stream >> temp;
      if(i==n_name) trans=temp; 
      if(i==n_strand) strand=temp;
      if(i==n_chrom) chrom=temp;
      if(i==n_exonStarts) starts=temp;
      if(i==n_exonEnds) ends=temp;
    }

    trans_positions[trans].set_values(chrom,starts,ends,strand);
  }
  file.close();
  cerr << "Done reading in transcript IDs" << endl;


  /**
   ** Read in the candidate fusion information
   **/
  file.open(argv[2]);
  if(!(file.good())){
    cerr << "Unable to open file " << argv[1] << endl;
    exit(1);
  }
  while ( getline (file,line) ){
    istringstream line_stream(line);
    string temp; string read; 
    int r_start; int r_end; int r_size;
    line_stream >> read >> r_start >> r_end >> temp >> r_size;
    //beginning...
    string trans_id1 ; int trans_end1; 
    line_stream >> trans_id1 >> trans_end1 ;
    //end
    string trans_id2 ; int trans_start2; 
    line_stream >> trans_id2 >> trans_start2 ;
    int overlap ; string strand;
    line_stream >> r_start >> r_end >> strand ;


    string chrom1, chrom2;
    int base1, base2;
    int offset1, offset2;

    smatch m; //extract the id
    regex_search(trans_id1,m,regex("_(EN[^_]*)__"));
    trans_positions[m[1].str()].get_position(trans_end1+1, true, base1, offset1); //to fix...
    chrom1=trans_positions[m[1].str()].chrom;

    regex_search(trans_id2,m,regex("_(EN[^_]*)__"));
    trans_positions[m[1].str()].get_position(trans_start2-1, false, base2, offset2); //to fix...
    chrom2=trans_positions[m[1].str()].chrom;

    cout << chrom1 << "\t" << base1 << "\t" << offset1 << "\t" 
	 << chrom2 << "\t" << base2 << "\t" << offset2 << "\t" << endl ;

  }
  return(0);
}
