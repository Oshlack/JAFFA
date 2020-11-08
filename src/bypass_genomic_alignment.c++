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


/** this function will print a "fake" .psl format alignment of the
    read to the genome of approx. 100bp. Alignments close
    to exon boundaries will be adjusted. **/
void print_alignment(string read, int read_length, 
		     int r_start, int r_end, string chrom, 
		     int base, int offset, 
		     string read_trans_strand, 
		     string trans_genome_strand, 
		     bool is_start){

  const int FLANK=200;
  const int FAKE_CHROM_LENGTH=1000000000;

  if(offset>=10) offset=0; // only adjust towards exon boundary if the distance to exon boundary is <10bp
  
  //set the start and end positions in the read
  int read_start=r_start-FLANK;
  int read_end=r_start+offset;
  if((read_trans_strand=="+")!=is_start){ //adjust for reverse compliments
    read_start=r_end-offset;
    read_end=r_end+FLANK;
  }
  //set the start and end positions in the genome
  int  genome_start=base-FLANK;
  int  genome_end=base+offset;
  if((trans_genome_strand=="+")!=is_start){ //adjust for reverse compliments
    genome_start=base-offset;
    genome_end=base+FLANK;
  }

  string strand="+";
  if(trans_genome_strand!=read_trans_strand) strand="-";
  
  //output in psl format
  cout << FLANK+offset << "\t0\t0\t0\t0\t0\t0\t0\t"
       << strand << "\t" << read << "\t" << read_length << "\t"
       << read_start << "\t" << read_end << "\t"
       << chrom << "\t" << FAKE_CHROM_LENGTH << "\t" 
       << genome_start << "\t" << genome_end << "\t" 
       << "1" << "\t" << FLANK+offset << ",\t" ;
  if(strand=="+"){ cout << read_start ; }
  else { cout << read_length-(read_start+FLANK+offset) ; }
  cout << ",\t" << genome_start << "," << endl;
}

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
    

    int base, offset;
    string trans;
    smatch m; 

    //extract the id for start of the fusion gene
    regex_search(trans_id1,m,regex("_(EN[^_]*)__"));
    trans=m[1].str();
    trans_positions[trans].get_position(trans_end1, true, base, offset); //to fix...
    print_alignment(read, r_size,
		    r_start, r_end,
		    trans_positions[trans].chrom,
		    base, offset, 
		    strand, 
		    trans_positions[trans].strand, 
		    true);

    //extract the id for end of the fusion gene
    regex_search(trans_id2,m,regex("_(EN[^_]*)__"));
    trans=m[1].str();
    trans_positions[trans].get_position(trans_start2, false, base, offset); //to fix...
    print_alignment(read, r_size,
		    r_start, r_end,
		    trans_positions[trans].chrom,
		    base, offset, 
		    strand, 
		    trans_positions[trans].strand, 
		    false);

  }
  return(0);
}
