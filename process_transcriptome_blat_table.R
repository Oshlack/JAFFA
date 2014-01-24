
library(GenomicRanges)
library(IRanges)

args = commandArgs(trailingOnly = TRUE)
blat_table=args[1]
output_file=args[2]
gap_size=as.numeric(args[3])

#get the gene names
ref_to_symbols=read.table(args[4],header=T,comment.char="/") #we want to read the first comment line
symbols<-as.character(ref_to_symbols$name2)
names(symbols)<-as.character(ref_to_symbols$name)

#get the gene positions
gene_positions<-data.frame(ref_to_symbols$chrom,ref_to_symbols$txStart,ref_to_symbols$txEnd,ref_to_symbols$name)
unique_genes=match(ref_to_symbols$name,ref_to_symbols$name)
gene_positions<-unique(gene_positions[unique_genes,])
rownames(gene_positions)<-gene_positions$ref_to_symbols.name
gene_positions$ref_to_symbols.name<-NULL
colnames(gene_positions)<-c("chrom","start","end")

#load the .psl file and split based on de novo transcript ID
rows<-read.table(file=blat_table,skip=5,stringsAsFactors=F, nrows=5)
classes<-sapply(rows, class)
blat_results=read.table(file=blat_table,skip=5,stringsAsFactors=F,colClasses = classes)

split_results=split(1:length(blat_results$V10),blat_results$V10)

#filter out all transcripts which are only covered by one reference transcript
interesting<-c()
j=1
for(x in split_results){
   if(j %% 500 == 0  ) { show(j) }
   if(length(x)!=1){ 
      # if there is more than one refSeq transcript which 
      # matches the de novo transcript
      # then check whether the two (or more) refSeq
      # transcripts cover different parts of the de novo
      # transcript
      maxs=which(blat_results$V13[x]==max(blat_results$V13[x]))
      mins=which(blat_results$V12[x]==min(blat_results$V12[x]))
      if(!any(mins %in% maxs)){ 
          ranges=unique(IRanges(blat_results$V12[x],blat_results$V13[x]))
	  u=reduce(ranges)
	  regions_accounted_for=sum(ranges %in% u)
#	  regions_accounted_for=length(subjectHits(findOverlaps(ranges,u,type="equal")))
	  regions=length(u)
	  if(regions_accounted_for < regions) { interesting<-c(interesting,j) }
      }
   }
   j=j+1
}

show(paste(length(interesting),"fusion candidates"))
#now start filtering based on overlap and distance between genes.
interesting2<-c()
candidates<-data.frame()
j=1
for(spx in split_results[interesting]){
      show(j) ; j=j+1
      x=blat_results[spx,]
      # here we are expecting the fasta ids to be in the format: ">hg19_annotation_geneName__otherStuff"
      temp_names=sapply(strsplit(x$V14,"__"),function(y){y[1]})
      gene_names=sapply(strsplit(temp_names,"_"),function(y){ paste(y[3:length(y)],collapse="_")  })
      gp=gene_positions[gene_names,]
      if(!any(!is.na(gp))){ #throw an error if we can't get the genomic range
	print(paste("Stopping because the gene,",gene_names,"can't be found in",args[4])) 
	stop() }

      #now we check for overlap between the genes
      granges=GRanges(gp$chrom,IRanges(gp$start,gp$end+gap_size)) #<-add gap size to here
      gunion=as.data.frame(reduce(granges))
      granges=as.data.frame(granges)
      #The GRanges findOverlaps is really slow, so we'll use the IRanges version:
      overlaps=findOverlaps(IRanges(granges$start,granges$end),IRanges(gunion$start,gunion$end))
      #now check the chrom. names
      overlaps=overlaps[gunion$seqnames[subjectHits(overlaps)]==granges$seqnames[queryHits(overlaps)]]
      ranges=IRanges(x$V12,x$V13)

      new_ranges<-IRanges()
      for(i in 1:length(unique(unique(subjectHits(overlaps))))){
      	    which_to_combine=queryHits(overlaps)[subjectHits(overlaps)==i]
            new_ranges=c(new_ranges,reduce(ranges[which_to_combine]))
      }
      cov=coverage(new_ranges)
      # now look for the hallmark of a fusion. We should see a coverage pattern
      # of 1 2 1 (1 gene , 2 genes overlap a few bases, then the other gene)
      covValue=runValue(cov)
      pattern_match<-function(n){ isTRUE(all.equal(covValue[n:(n+2)],c(1,2,1))) }

      if(length(covValue)>=3){ 
         temp_pos=which(sapply( 1:(length(runValue(cov))-2),pattern_match))+1
	 temp_pos=temp_pos[which(runLength(cov)[temp_pos]==min(runLength(cov)[temp_pos]))[1]]
         ov=runLength(cov)[temp_pos] #how big is the overlap between genes
      	 if(!is.na(ov) & ov<15 & ov>0){ #when two gene overlap by more than 15 bases it probably an assembler mistake
	       # set the fusion point to the first base which has
	       # a coverage of 2 ??
	       base_before=sum(runLength(cov)[1:(temp_pos-1)])
	       base_after=sum(runLength(cov)[1:temp_pos])+1
	       x$fusionPointMin<-rep(base_before,dim(x)[1])
	       x$fusionPointMax<-rep(base_after,dim(x)[1])
	       get_gene_symbol_at_base<-function(b){
	           g=which( (x$V12 < b) & (b < x$V13 ))
		   g=g[which(x[g,]$V1==max(x[g,]$V1))[1]]
		   symbols[gene_names[g]] 
	       }
	       x$fusion_genes<-paste(get_gene_symbol_at_base(base_before),
				     get_gene_symbol_at_base(base_after),sep=":")

	       show(x)
	       show(cov)
	       candidates<-rbind(candidates,x)
	       interesting2<-c(interesting2,interesting[j])
	}
     }
}

show(paste(length(interesting2),"fusions which look real"))
write.table(candidates,file=output_file,sep="\t",row.names=F,quote=F,col.names=F)

