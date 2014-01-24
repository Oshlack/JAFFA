
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
full_list=full_list[order(full_list$`gap (kb)`,decreasing=T),]
v=split(full_list,full_list$classification)
full_list=rbind(v[["HighConfidence"]],v[["MediumConfidence"]],v[["LowConfidence"]],v[["PotentialRegularTranscript"]])
write.csv(full_list,paste(out_name,".csv",sep=""),row.names=F)


##########  get the sequences and igv snap-shots ##################

system(paste("rm ",out_name,".fa",sep=""))
system(paste("rm ",out_name,"_igv_command ; mkdir igv_snapshots",sep=""))
system("echo \"@HD\tVN:1.0\tSO:unsorted\" > header.temp ")
system("rm body.temp ")
for(i in 1:dim(full_list)[1]){

      samp=full_list[i,]$sample
      trans=full_list[i,]$contig
      pos=full_list[i,]$`contig break min`
      fus=full_list[i,]$`fusion genes`

      ### code to make the fasta file
      fusions_file=paste(samp,"/",samp,".fusions.fa",sep="")
      new_id=paste(samp,fus,trans,sep="---")
    
      system(paste("echo \">",new_id,"\" >> ",out_name,".fa",sep=""))
      comm=paste("grep -A1 \"^>",trans," \" ",fusions_file," | grep -v \"^>\" >> ",out_name,".fa",sep="")
      system(comm)

      ### code to make the IGV batch script
      add_line_to_file<-function(line,file=paste(out_name,"igv_command",sep="_")){
         system(paste("echo \"",line,"\" >> ",file,sep=""))
      }

      add_line_to_file("new")
      add_line_to_file(paste("snapshotDirectory ",getwd(),"/igv_snapshots",sep=""))
      add_line_to_file(paste("genome ",getwd(),"/",fusions_file,sep=""))
      add_line_to_file(paste("load ",getwd(),"/",samp,"/",samp,".sam.sorted.bam",sep=""))
      add_line_to_file(paste("goto ",trans,sep=""))
      add_line_to_file(paste("region ",trans," ",pos," ",pos,sep=""))
      add_line_to_file(paste("snapshot ",new_id,"-FULL.png",sep=""))
      add_line_to_file(paste("goto ",trans,":",pos-50,"-",pos+50,sep=""))
      add_line_to_file(paste("snapshot ",new_id,".png",sep=""))

      ### code to concatenate the bam files
      bam=paste(samp,"/",samp,".sorted.bam",sep="")

      ## make a new header
      cmd=paste("samtools view -H ",bam," | grep ",trans,sep="")
      cmd=paste(cmd," | sed \'s/",trans,"/",new_id,"/g\' >> header.temp",sep="")
      show(cmd)
      system(cmd)

      ## get the alignments
      cmd=paste("samtools view ",bam," | grep ",trans,sep="")
      cmd=paste(cmd," | sed \'s/",trans,"/",new_id,"/g\' >> body.temp",sep="")
      show(cmd)
      system(cmd)
}

system(paste("cat header.temp > ",out_name,".sam",sep=""))
system(paste("cat body.temp >> ",out_name,".sam",sep=""))
system(paste("sam2bam ",out_name,".sam",sep=""))

print("finished")

