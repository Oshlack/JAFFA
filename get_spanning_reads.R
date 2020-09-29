options(echo=FALSE)
args = commandArgs(trailingOnly = TRUE)
sample=args[1]
candidates_table=args[2]
bam_data_file=args[3]
output_table=args[4]
read_length_list=args[5]  #as.numeric(args[5])
hang=as.numeric(args[6])

show(candidates_table)
show(output_table)

fus_cols_of_interest=c("transcript","break_min","break_max","fusion_genes","contig_length")
fusions=read.delim(candidates_table,header=F,stringsAsFactors=F)
colnames(fusions)<-fus_cols_of_interest
fusions=fusions[,fus_cols_of_interest]

alignments_table=data.frame(read.delim(bam_data_file,header=F,stringsAsFactors=F),
			    read.delim(read_length_list,header=F,stringsAsFactors=F))
sal=split(1:length(alignments_table$V1),alignments_table$V1)

get_spanning_pairs<-function(x){
	contig_length=as.numeric(x[5])
	pos=as.numeric(x[2])
	trans=x[1]
	alignments=alignments_table[sal[[trans]],]
	read_lengths=alignments[,4]
	mean_length=mean(read_lengths)
	if(is.na(mean_length))
		return(0) #case for no aligments
	if((pos<mean_length*2)||((contig_length-pos)<mean_length*2))
		return(0) #flag short contigs
	this_below=alignments$V2<(pos-read_lengths)
	pair_above=alignments$V3>pos
	sum(this_below&pair_above)
}
fusions$spanning_pairs=apply(fusions,1,get_spanning_pairs)

get_spanning_reads<-function(x){
	trans=x[1]
	pos=as.numeric(x[2])
	alignments=alignments_table[sal[[trans]],]
	read_lengths=alignments[,4]
	below_fusion=alignments$V2>(pos+hang-read_lengths)
	above_fusion=alignments$V2<(pos-hang)
        sum(below_fusion&above_fusion)
}
fusions$spanning_reads=apply(fusions,1,get_spanning_reads)

write.table(fusions[,c(1:4,6:7)],output_table,quote=F,row.names=F,sep="\t")
