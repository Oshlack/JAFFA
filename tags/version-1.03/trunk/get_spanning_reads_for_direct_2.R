#options(echo=FALSE)

options(stringsAsFactors=F)
args = commandArgs(trailingOnly = TRUE)

## load the table of candidates
candidates=read.delim(args[1],header=F)

## fix the ids to be alphabetically ordered
candidates$V4=sapply(strsplit(candidates$V4,":"),function(x){ paste(sort(x),collapse=":") })
colnames(candidates)<-c("transcript","break_min","break_max","fusion_genes","length")
candidates$length<-NULL
candidates$spanning_pairs<-0
candidates$spanning_reads<-1

## load the count data
sr=read.delim(args[2],header=F)

for( n in 1:length(sr[,1])){
   entries=candidates$fusion_genes==sr[n,1]
   spanning_rs=sr[n,2]
   if(spanning_rs>0){ #distribute evelyn amongst instances of this fusion in the list
      scounts=table(rep_len(1:sum(entries),spanning_rs))
      if(length(scounts)<sum(entries))
         scounts=c(scounts,rep(0,sum(entries)-length(scounts)))
      candidates$spanning_pairs[entries]<-scounts
   }
}

write.table(candidates,args[3],sep="\t",row.names=F,col.names=T,quote=F)
