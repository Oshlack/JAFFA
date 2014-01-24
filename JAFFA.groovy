
load "JAFFA_stages.groovy"

// The actual pipeline.
Bpipe.run{ fastq_filename_pattern * [ 
       			      prepare_reads + 
			      run_assembly +
			      align_transcripts_to_annotation + 
			      filter_transcripts + 
			      extract_fusion_sequences + 
			      map_reads + 
			      get_spanning_reads +
			      align_transcripts_to_genome + 
			      get_final_list 
			      //clean_up 
			      ] + compile_all_results }
