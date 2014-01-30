
# make a list of the summary files from the pipe-line
args = commandArgs(trailingOnly = TRUE)
out_name=args[1]
sfiles=args[2:length(args)]
dir=sapply(strsplit(sfiles,"/"),function(x){x[1]})
show(dir)

summary_files=paste(getwd(),sfiles,sep="/")
show(summary_files)
#check to see if the summary file was even made:
exists=file.exists(summary_files)
show(exists)
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
show(colnames(full_list))
full_list=full_list[,c("sample","fusion genes","chrom1","base1","chrom2","base2",
	     "gap (kb)","spanning pairs","spanning reads",
	     "inframe","aligns","rearrangement",
	     "contig","contig break","classification")]
full_list=full_list[order(as.numeric(full_list$`spanning reads`),decreasing=T),]
v=split(full_list,full_list$classification)
full_list=rbind(v[["HighConfidence"]],v[["MediumConfidence"]],v[["LowConfidence"]],v[["PotentialRegularTranscript"]])
write.csv(full_list,paste(out_name,".csv",sep=""),row.names=F)



