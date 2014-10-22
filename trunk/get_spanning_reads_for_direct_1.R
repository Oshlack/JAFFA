#options(echo=FALSE)

options(stringsAsFactors=F)

args = commandArgs(trailingOnly = TRUE)

## load the table of candidates
candidates=read.delim(args[1],header=F)

## load the reference annotation which maps gene IDs to 
ref=read.delim(args[2],stringsAsFactors=F)[,c("name","name2")]

## load the actually gene ids as they are stored in the bam file
ref_contigs=read.delim(args[3],header=F,stringsAsFactors=F)[,1]

## replace the contig ids in the ref table
names(ref_contigs)<-sapply(strsplit(ref_contigs,"_"),function(x){x[3]})
ref$name<-ref_contigs[ref$name]
sref=split(ref$name,ref$name2)
sref=sapply(sref,function(x){paste(x,collapse=" ")})


#now make the table which need to be run on the .bam file
genes=unique(candidates$V4)
sgenes=lapply(strsplit(genes,":"),function(x){sort(x)})
sgenes=unique(sgenes)
g1=sapply(sgenes,function(x){ sref[[x[1]]]  })
g2=sapply(sgenes,function(x){ sref[[x[2]]]  })
genes_new=sapply(sgenes,function(x){ paste(x[1],x[2],sep=":") })
write.table(data.frame(genes_new,g1,g2),args[4],sep="?",row.names=F,col.names=F,quote=F)
