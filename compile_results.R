options(echo=F)

# make a list of the summary files from the pipe-line
args = commandArgs(trailingOnly = TRUE)
out_name=args[1]
sfiles=args[2:length(args)]
summary_files=sfiles
sfiles=sapply(strsplit(sfiles,"/"),function(x){paste(x[(length(x)-1):length(x)],collapse="/")})
dir=sapply(strsplit(sfiles,"/"),function(x){x[1]})
message("Compiling the results from:")
message(paste(dir,collapse=" "))

#check to see if the summary file was even made:
#exists=file.exists(summary_files)
info = file.info(summary_files)
exists = summary_files %in% rownames(info[!is.na(info$size) & info$size>0, ])
if(sum(!exists)>0){
   show("No fusions were found for the following samples:")
   show(dir[!exists])
}
dir=dir[exists]
summary_files=summary_files[exists]

# parse each file:

full_list<-data.frame()

for(i in 1:length(dir)){
      single_list=read.delim(summary_files[i],stringsAsFactors=F)
      single_list$sample<-rep(dir[i],length(single_list$transcript))
      #rearrange
      n=dim(single_list)[2]
      full_list<-rbind(full_list,single_list)
}
colnames(full_list)<-gsub("_"," ",colnames(full_list))
colnames(full_list)[colnames(full_list)=="transcript"]<-"contig"
colnames(full_list)[colnames(full_list)=="gap"]<-"gap (kb)"

#reorder
full_list=full_list[,c("sample","fusion genes","chrom1","base1","strand1","chrom2","base2","strand2",
	     "gap (kb)","spanning pairs","spanning reads",
	     "inframe","aligns","rearrangement",
	     "contig","contig break","classification","known")]
#	     "contig","contig break","classification","known","geneCounts1","geneCounts2","exon1","exon2")]

#first order on the gap size (useful for direct mode where the spanning reads are mostly just 1)
full_list=full_list[order(as.numeric(full_list$`gap (kb)`),decreasing=F),] 

#now order on the number of supporting reads
supporting_reads=full_list$`spanning pairs`
supporting_reads[supporting_reads=="-"]<-0
supporting_reads=as.numeric(supporting_reads)+full_list$`spanning reads`

#full_list=full_list[order(as.numeric(full_list$`spanning reads`),decreasing=T),]
full_list=full_list[order(supporting_reads,decreasing=T),]
v=split(full_list,full_list$classification)

#then order on classification
full_list=rbind(v[["HighConfidence"]],v[["MediumConfidence"]],v[["LowConfidence"]],v[["PotentialTransSplicing"]],
	v[["PotentialRunThrough"]])
write.csv(full_list,paste(out_name,".csv",sep=""),row.names=F)

message(paste("Done writing output ",out_name,".csv",sep=""))

