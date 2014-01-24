args = commandArgs(trailingOnly = TRUE)
sample=args[1]
candidates_table=args[2]
bam_data_file=args[3]
output_table=args[4]
read_length=as.numeric(args[5])
hang=as.numeric(args[6])

show(candidates_table)
show(output_table)

fusions=read.delim(candidates_table,header=F,stringsAsFactors=F)
fusions=unique(fusions[,c(10,22:24,11)])
colnames(fusions)<-c("transcript","break_min","break_max","fusion_genes","contig_length")

alignments_table=read.delim(bam_data_file,header=F,stringsAsFactors=F)
sal=split(1:length(alignments_table$V1),alignments_table$V1)

get_spanning_pairs<-function(x){
	contig_length=as.numeric(x[5])
	pos=as.numeric(x[2])
	if((pos<read_length*2)||((contig_length-pos)<read_length*2))
		return("-") #flag short contigs
	trans=x[1]
	alignments=alignments_table[sal[[trans]],]
	this_below=alignments$V2<(pos-read_length)
	pair_above=alignments$V3>pos
	dim(alignments[this_below&pair_above,])[1]
}
fusions$spanning_pairs=apply(fusions,1,get_spanning_pairs)

get_spanning_reads<-function(x){
	trans=x[1]
	pos=as.numeric(x[2])
	alignments=alignments_table[sal[[trans]],]
	below_fusion=alignments$V2>(pos+hang-read_length)
	above_fusion=alignments$V2<(pos-hang)
        dim(alignments[below_fusion&above_fusion,])[1]
}
fusions$spanning_reads=apply(fusions,1,get_spanning_reads)

write.table(fusions[,c(1:4,6:7)],output_table,quote=F,row.names=F,sep="\t")
