options(echo=FALSE)

###################################################################################################
# This script takes as input a blat .psl format file from transcripts aligned to the human genome
# (produced in the previous step of the pipeline), read coverage information (produced in the
# second last pipeline stage) and a table of reference annotation information (provided by the pipeline).
# It does the following:
# - works out the position of the break-point in the genome
# - works out the genome gap size
# - works out if there has been a genomic rearrangement or inversion
# - calculates whether the break aligns with an exon-boundary
# - calculates whether the fusion is inframe
# - merges all this information with number of spanning reads and pairs from an earlier stage
# - filters out low coverage / small genomic gap transcripts (?)
# - outputs a table (.summary) of final candidates
#
# written by N. Davidson (nadia.davidson@mcri.edu.au) 9/12/2013
###############################################################################################

# this script can be run like R --args <all the arguments below> < make_final_table.R
args = commandArgs(trailingOnly = TRUE)
blat_table_file=args[1]   # transcripts aligned to the human genome, will be <X>_genome.psl
fusion_info_file=args[2]  # read coverage for the alignments, will be <X>.reads 
trans_table_file=args[3]  # a reference annotation file
gapmin=as.numeric(args[4]) # minimum genomic gap of the transcriptional break-point (in bases). 
exclude=args[5]		  # which "classifications" to remove"
output_file=args[6]       # name of the output file, will be <X>.summary

#maximum number of bases discrepancy between genomic alignment and exons boudary for the break-point to be corrected
OVERHANG_MAX=20
REGGAP=200 #fusions with less than this kb gap and no rearanngments will be flagged as regular

#load all the input files to data.frames
fusion_info<-read.delim(fusion_info_file,stringsAsFactors=F)
transTable=read.table(trans_table_file,header=T,stringsAsFactors=F,comment.char="/")
blat_table<-read.delim(blat_table_file,stringsAsFactors=F,skip=5,header=F)
sgb=split(blat_table,blat_table$V10)

#############  check the contig location in the genome ###########
get_break_pos<-function(n){
   # get the transcript to genome blat results for the nth transcript
   # if there was no genome alignment return nothing
   contig=as.character(fusion_info$transcript[n])
   if(!contig %in% names(sgb)){ return() }
   x=sgb[[contig]]
   break_in_trans1=fusion_info$break_min[n]
   break_in_trans2=fusion_info$break_max[n] #names need to change####<---

   #below is to adjust the break-point based on the genome results.
   get_one<-function(j,brk,is_start){
      start=as.numeric(unlist(strsplit(x[j,]$V20,",")))
      width=as.numeric(unlist(strsplit(x[j,]$V19,",")))
      genome_starts=as.numeric(unlist(strsplit(x[j,]$V21,",")))
      genome_dir=(x[j,]$V9=="+")
      if(!genome_dir){
         start=x[j,]$V11-start-width
	 genome_starts=genome_starts+width
      }
      block=(start <= brk) & ((start + width) >= brk)
      if(sum(block)==0){ return() }
      offset=brk-start[block]
      if(genome_dir){
	genome_pos=genome_starts[block]+offset
      } else { genome_pos=genome_starts[block] - offset + 1}
      genome_chrom=x[j,]$V14

      #now get the direction the alignment is running in
      if(is_start&(brk-x[j,]$V12)<15) return() #magic numbers ... <--------
      if((!is_start)&(x[j,]$V13-brk)<15) return()
      length=x[j,1]
      pid=1-x[j,2]/x[j,1]
      data.frame(genome_chrom,genome_pos,genome_dir,is_start,brk,length,pid)
   }
   res1=lapply(1:(dim(x)[1]),get_one,break_in_trans1,TRUE)
   res2=lapply(1:(dim(x)[1]),get_one,break_in_trans2,FALSE)
   res=rbind(do.call(rbind.data.frame,res1),do.call(rbind.data.frame,res2))
}
show("Getting the location of fusion transcripts in the genome..")
genome_pos=lapply(1:length(fusion_info$transcript),get_break_pos)


#############  get the genomic gap size ###########
# and also filter out transcripts which haven't fully aligned.

#get the best starting match to the genom
#and the best end. Return NA for any transcripts
#which don't have a start and an end.
check_gap<-function(x){
	if(length(x)==0) return()
	chrom=outer(x$genome_chrom,x$genome_chrom,"==")
	start_end=outer(x$is_start,x$is_start,"!=")
	# we are not interested in cases where both sides don't align
	if((any(x$is_start)&any(!x$is_start))==FALSE){ return() }
	#function to check if a rearrangement has occured
	rearr<-function(r){
		r$genome_chrom[1]!=r$genome_chrom[2] ||
                r$genome_dir[1]!=r$genome_dir[2] || #inversion
                ((r$genome_pos[!r$is_start] > r$genome_pos[r$is_start]) != r$genome_dir[1] ) }
		# rearrangement like FGFR3-TACC3 
	#below is how we choose which start and end in the genome to use (in case of
	#multiple starts and ends..	
	get_best_fit<-function(){
		#if there is a start and end close together (less the 200kb default), use them..
                starts=x[x$is_start,]
                ends=x[!x$is_start,]
		close_matrix=(abs(outer(starts$genome_pos,ends$genome_pos,"-"))<REGGAP*1000) &
		      outer(starts$genome_chrom,ends$genome_chrom,"==") &
		      outer(starts$genome_dir,ends$genome_dir,"==")
		close_ind=which(close_matrix,arr.ind=TRUE)
		if(length(close_ind)>0){
		   # use the first instance (not completely correct...)
		   xstart=starts[close_ind[1,1],]
		   xend=ends[close_ind[1,2],]
		} else { 
		   #otherwise choose the two longest
                   xstart=starts[which.max(starts$length*starts$pid),]
                   xend=ends[which.max(ends$length*ends$pid),]
		}
                temp<-rbind(xstart,xend)
                temp$gap<-Inf
		if(length(unique(temp$genome_chrom))==1) temp$gap = abs(diff(temp$genome_pos))/1000
                temp$rearrangement <- rearr(temp)
                return(temp)  }
	if(all(!(start_end&chrom))){ #if there's no start and end on the same chrom. 
		return(get_best_fit())  } # otherwise lets check the gap size...
	distance=abs(outer((x$genome_pos),(x$genome_pos),"-"))
	distance[!(start_end&chrom)]<-NA
	gap=min(distance,na.rm=T)
	min_gap_idx=which(distance==gap,arr.ind=TRUE)[1,]
	#return empty if we don't pass the gap requirement???
	if(gap < gapmin & !rearr(x[min_gap_idx,]) ) return()
	# otherwise lets work out which is the best match
	# (it might not be the one which is closest)
	return( get_best_fit() )
}
message("Calculating gap size in the genome...")
new_genome_pos=lapply(genome_pos,check_gap)

#############  compare against the know annotation ###########
get_frame_info<-function(x){
	if(length(x)==0) return()
	do_one_row<-function(j){
		chrom<-as.character(x[j,]$genome_chrom)
		pos<-x[j,]$genome_pos
		correct_chrom=transTable$chrom==chrom
      		correct_pos=(transTable$txStart<=pos)&(transTable$txEnd>=pos)
		####  if we don't find a match then return..   
		if(sum(correct_pos&correct_chrom)==0) return()
		gene=transTable[correct_pos&correct_chrom,]
		#get the exon positions for each transcripts
		starts=strsplit(gene$exonStart,",")
		ends=strsplit(gene$exonEnd,",")
		dists=lapply(1:length(starts),function(y){ 
		   min(abs(as.integer(c(starts[[y]],ends[[y]]))-pos)) 
                })
		# select the transcript with the closest exon to the break point
		trans=which.min(unlist(dists))
		gene=gene[trans,]
		#which exon in this transc
		dist=sort(c(as.integer(starts[[trans]])+1,as.integer(ends[[trans]]))-pos)
		closest=which(abs(dist)==min(abs(dist)))[1]
		overhang<-dist[closest]   
		closestExon=ceiling(closest[1]/2)
		is_actually_the_start=((gene$strand=="+")==x[j,]$genome_dir)==x[j,]$is_start
		frames=as.integer(unlist(strsplit(as.character(gene$exonFrames),",")))
		frame=frames[closestExon]
		if(is_actually_the_start & frame>=0 ){ #then get the frame of the next exon (if this exon is coding)
		   if(gene$strand=="+") dir=1 else dir=-1
		   nextExon=closestExon + dir
		   if((nextExon < 1) | (nextExon > length(frames))){ 
		      frame=-1 # return non-coding frame if we hit the end of the transcript
		   } else { frame=frames[nextExon] }
		}
		return(data.frame(x[j,],overhang,frame,is_actually_the_start));
	}
   	res=do.call(rbind.data.frame,lapply(1:(dim(x)[1]),do_one_row))   
	if(dim(res)[1]!=2) return()
	### check frame #####
	res$inFrame=NA
	### fix/check overhang... #####
	dir=c(1,1)
        dir[res$genome_dir==F]=-1
        new_break=res$brk+dir*res$overhang
	res$aligns=F
	if(((new_break[!res$is_start]-new_break[res$is_start])==1) &&
		(sum(abs(res$overhang))<OVERHANG_MAX)) { 
	        res$genome_pos=res$genome_pos+res$overhang
		res$brk=new_break
		res$aligns=T
		#now check frame
		if( all(res$frame>=0) & diff(res$frame)==0 & diff(res$is_actually_the_start)!=0){
		   res$inFrame=T
 		} else {
                   res$inFrame=F
		}
	}
	res
}
message("Checking if the fusions are in frame...")
new_new_genome_pos=lapply(new_genome_pos,get_frame_info)

#############  format nicely  ###########
format_positions<-function(x){
	if(length(x)==0){
	   res=data.frame(NA,NA,NA,0,NA,NA,NA,NA,NA)
	   colnames(res)<-c("contig_break","chrom1","base1","chrom2","base2","gap","rearrangement","aligns","inframe")
	   return(res)
	 }
	#there should only be 2 or 0 genomic positions at this point.
	chrom1=x$genome_chrom[1] ; chrom2=x$genome_chrom[2]
	base1=x$genome_pos[1] ;  base2=x$genome_pos[2]
	gap=as.numeric(as.character(x$gap[1])) ; rearrangement=x$rearrangement[1] 
	aligns=all(x$aligns) ; inframe=x$inFrame[1]
	contig_break=min(x$brk)[1] 
	return(data.frame(contig_break,chrom1,base1,chrom2,base2,gap,rearrangement,aligns,inframe))
}
genome_info<-lapply(new_new_genome_pos,format_positions)

message("Merging with read coverage data...")

#############  merge with read coverage and gene name information  ###########
result=cbind(fusion_info[,c(1,4:6)],do.call(rbind.data.frame,genome_info))
fix_names<-function(x){ paste(sort(unlist(strsplit(x,":"))),collapse=":") }
result$fusion_genes<-sapply(result$fusion_genes,fix_names)

#group fusions by break-point
break_string=sapply(paste(result[,"chrom1"],result[,"base1"],":",result[,"chrom2"],result[,"base2"],sep=""),fix_names)
r=split(result,break_string) 
merge_result<-function(x){
   is_short = all(x$spanning_pairs=="-")
   x$spanning_pairs[x$spanning_pairs=="-"]=0
   x$spanning_pairs=as.numeric(x$spanning_pairs)
   
   #choose the representative
   x=x[order(x$spanning_pairs,decreasing=T),]
   x=x[order(x$fusion_genes),]
   x_rep=x[1,]
   if(is_short){ x_rep$spanning_pairs= "-" } else {
      x_rep$spanning_pairs=sum(x$spanning_pairs) }
   x_rep$spanning_reads=sum(x$spanning_reads)
   x_rep
}
result=do.call(rbind.data.frame,lapply(r,merge_result))
cand<-result[!is.na(result$rearrangement),]

########### now classify the candidates #########################

cand=cand[cand$gap>(gapmin/1000),] #remove anything with a gap below 10kb
cand$classification<-"NoSupport"
spanP=cand$spanning_pairs>0
spanR=cand$spanning_reads>0
cand$classification[ spanP & spanR ]<-"LowConfidence"
cand$classification[ cand$aligns & (spanP | spanR ) ]<-"MediumConfidence"
cand$classification[ cand$aligns & spanP & spanR ]<-"HighConfidence"
cand$classification[ (cand$gap<REGGAP) & ( spanP | spanR ) & !cand$rearrangement ]<-"PotentialRegularTranscript"

#remove any group in the exclude list
exclude=unlist(strsplit(exclude," "))
cand=cand[ !cand$classification %in% exclude, ]

write.table(cand,output_file,row.names=F,quote=F,sep="\t")
message("Done producing summary file")