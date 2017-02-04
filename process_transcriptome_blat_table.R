
library(IRanges)

args = commandArgs(trailingOnly = TRUE)
blat_table=args[1]
output_file=args[2]
gap_size=as.numeric(args[3])

MAX_OVERLAP=15 #maximum+1 number of bases that two genes can share at the break-point (<15)

#get the gene names
ref_to_symbols=read.table(args[4],header=T,comment.char="/",stringsAsFactors=F) #we want to read the first comment line
symbols<-as.character(ref_to_symbols$name2)
names(symbols)<-as.character(ref_to_symbols$name)

#get the gene positions
gene_positions<-data.frame(ref_to_symbols$chrom,ref_to_symbols$txStart,ref_to_symbols$txEnd,ref_to_symbols$name,stringsAsFactors=F)
unique_genes=match(ref_to_symbols$name,ref_to_symbols$name)
gene_positions<-unique(gene_positions[unique_genes,])
rownames(gene_positions)<-as.factor(gene_positions$ref_to_symbols.name)
gene_positions$ref_to_symbols.name<-NULL
colnames(gene_positions)<-c("chrom","start","end")

#load the .psl file and split based on de novo transcript ID
rows<-read.table(file=blat_table,skip=5,stringsAsFactors=F, nrows=5)
classes<-sapply(rows, class)
blat_results=read.table(file=blat_table,skip=5,stringsAsFactors=F,colClasses = classes,comment.char="")

split_results=split(1:length(blat_results$V10),blat_results$V10)

#filter out all transcripts which are only covered by one reference transcript
j<<-1
multi_gene<-function(x){
   if(j %% 10 == 0  ) { show(j) } ; j<<-j+1 

   #if one reference transcript covered the whole range then return.
   if(length(x)==1) return()
   starts=blat_results$V12[x]
   ends=blat_results$V13[x]
   maxs=which(ends==max(ends))
   mins=which(starts==min(starts))
   if(any(mins %in% maxs)) return()

   #now get just the non-redundant set of transcripts 
   ranges=IRanges(starts,ends)
   overs=findOverlaps(ranges,type="within",drop.self=TRUE,drop.redundant=TRUE,select="arbitrary")
   reduced=is.na(overs)
   regions=ranges[reduced]

   #where are these regions in the genome?
   genes=blat_results$V14[x[reduced]] 
   # here we are expecting the fasta ids to be in the format: ">hg19_annotation_geneName__otherStuff"
   temp_names=sapply(strsplit(genes,"__"),function(y){y[1]})
   gene_names=sapply(strsplit(temp_names,"_"),function(y){ paste(y[3:length(y)],collapse="_")  })
   gp=gene_positions[gene_names,]
   if(!any(!is.na(gp))){ #throw an error if we can't get the genomic range
        print(paste("Stopping because the gene,",gene_names,"can't be found in",args[4]))
        stop() }

    #split by chrom 
    split_chroms=split(1:length(gene_names),gp$chrom)
    #split by region on chromosome
    new_ranges<-IRanges()
    for(sc in split_chroms){
	if(length(sc)==1){ 
	   new_ranges=c(new_ranges,regions[sc])
	} else {
	   ir=IRanges(gp$start[sc],gp$end[sc]+gap_size)
      	   iru=reduce(ir)
      	   overlaps=findOverlaps(ir,iru,select="arbitrary")
	   for(w in split(sc,overlaps)){ new_ranges=c(new_ranges,reduce(regions[w])) }
        }
     }
     cov=coverage(new_ranges)

     # now look for the hallmark of a fusion. We should see a coverage pattern
     # of 1 2 1 (1 gene , 2 genes overlap a few bases, then the other gene)
     covValue=runValue(cov)
     pattern_match<-function(n){ isTRUE(all.equal(covValue[n:(n+2)],c(1,2,1))) }

     if(length(covValue)>=3){ 
         temp_pos=which(sapply( 1:(length(runValue(cov))-2),pattern_match))+1
	 temp_pos=suppressWarnings(temp_pos[which(runLength(cov)[temp_pos]==min(runLength(cov)[temp_pos]))[1]])
         ov=runLength(cov)[temp_pos] #how big is the overlap between genes
      	 if(!is.na(ov) & ov<MAX_OVERLAP & ov>0){ #when two gene overlap by more than 14 bases it probably an assembler mistake
	       # set the fusion point to the first base which has
	       # a coverage of 2 ??
	       base_before=sum(runLength(cov)[1:(temp_pos-1)])
	       base_after=sum(runLength(cov)[1:temp_pos])+1

	       rbr=blat_results[x[reduced],]

	       get_gene_symbol_at_base<-function(b){
	           g=which( (rbr$V12 < b) & (b < rbr$V13 ))
		   g=g[which(rbr[g,]$V1==max(rbr[g,]$V1))[1]]
		   symbols[gene_names[g]] 
	       }
	       fusion_genes<-paste(get_gene_symbol_at_base(base_before),
				   get_gene_symbol_at_base(base_after),sep=":")

	       show(rbr)
	       show(cov)
	       return(data.frame(rbr$V10[1],base_before,base_after,fusion_genes,rbr$V11[1]))
	}
     }
}
candidates=do.call("rbind",lapply(split_results,multi_gene))

show(paste(dim(candidates)[1],"fusion preliminary candidates identified"))
write.table(candidates,file=output_file,sep="\t",row.names=F,quote=F,col.names=F)

