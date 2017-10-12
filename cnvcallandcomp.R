library(copynumber)
library(sequenza)
library(ASCAT)

######################MY FUNCTIONS##############

printsegmentation = function (datai,outdir='segmentation_chromosomes',chromosomes) {
	dir.create(outdir)
    for (chr in 1:length(chromosomes)) {
        chromosome=chromosomes[chr]
        pdf(paste0(outdir,"/chr",chromosomes[chr],"_view.pdf"))
        try(chromosome.view(mut.tab = datai$mutations[[chromosome]], baf.windows = datai$BAF[[chromosome]], ratio.windows = datai$ratio[[chromosome]], min.N.ratio = 1, segments = datai$segments[[chromosome]], main = datai$chromosomes[chromosome]))
        dev.off()
    }
}

##Function to get the data for a chromosome
prepDataChr = function (gc.stats,file,chr,chr.vect,gc.vect,min.reads.baf,window,overlap,gz= TRUE) {

    file.lines <- gc.stats$file.metrics[which(chr.vect == 
                                                chr), ] ##get the lines from the gc stats
    seqz.data <- read.seqz(file, gz = gz, n.lines = c(file.lines$start, 
                                                      file.lines$end)) #read them
    seqz.data$adjusted.ratio <- round(seqz.data$depth.ratio/gc.vect[as.character(seqz.data$GC.percent)], 
                                      3)#adjusts the depth with the GC content
    seqz.hom <- seqz.data$zygosity.normal == "hom" #this is analyzed in the python script apparently
    seqz.het <- seqz.data[!seqz.hom, ] #this is analyzed in the python script apparently
    het.filt <- seqz.het$good.reads >= min.reads.baf #this filters the heterozigous that have less than min.reads.baf out. What about the normal??
    seqz.r.win <- windowValues(x = seqz.data$adjusted.ratio, 
                               positions = seqz.data$position, chromosomes = seqz.data$chromosome, 
                               window = window, overlap = overlap, weight = seqz.data$depth.normal) ##Binned adjusted depth data every 1M positions
    
    if (nrow(seqz.het) > 0) { ##We have heterozygous positions
        seqz.b.win <- windowValues(x = seqz.het$Bf[het.filt], 
                             positions = seqz.het$position[het.filt], chromosomes = seqz.het$chromosome[het.filt], 
                             window = window, overlap = overlap, weight = seqz.het$good.reads[het.filt]) ##Binned BAF very 1M positions
    } else { #just one segment              
        seqz.b.win <- list()
        seqz.b.win[[1]] <- data.frame(start = min(seqz.data$position, na.rm = TRUE), end = max(seqz.data$position, 
                                na.rm = TRUE), mean = 0.5, q0 = 0.5, q1 = 0.5, N = 1)
    }
    return(list(file.lines=file.lines,seqz.data=seqz.data,seqz.hom=seqz.hom,seqz.het=seqz.het,het.filt=het.filt,seqz.r.win=seqz.r.win,seqz.b.win=seqz.b.win))
}

##Function to segment multiple samples at once

find.breaks.multi= function (seqz.baf1, seqz.baf2, gamma= 80, penalty=25, verbose = FALSE, seg.algo = "asmultipcf", 
                            assembly= "hg19", gender=c("XX","XX"), sexchromosomes=c("X","Y"), refine=FALSE, ...) 
{
    require(copynumber)
    require(sequenza)
    require(ASCAT)
	
	dummyarm=""
	
    chromosome1 <- gsub(x = seqz.baf1$chromosome, pattern = "chr", replacement = "")
    chromosome2 <- gsub(x = seqz.baf2$chromosome, pattern = "chr", replacement = "")

    logR1 = data.frame(chrom = chromosome1, pos = seqz.baf1$position, 
        s1 = log2(seqz.baf1$adjusted.ratio))
    logR2 = data.frame(chrom = chromosome2, pos = seqz.baf2$position,
        s2= log2(seqz.baf2$adjusted.ratio))

    if (seg.algo == "asmultipcf") {
    
        #Merge the two samples and remove NAs in the common column
        logR=merge(logR1,logR2[,-1],by = "pos",all=TRUE)
        ichrom=unique(logR$chrom)
        ichrom=ichrom[!is.na(ichrom)]
        logR$chrom=ichrom
        logR=logR[,c(2,1,3,4)] #correct order
       
        ##Generate BAF dataframes 
        BAF1 = data.frame(chrom = chromosome1, pos = seqz.baf1$position, 
            s1 = seqz.baf1$Bf)
        BAF2 = data.frame(chrom = chromosome2, pos = seqz.baf2$position, 
            s2 = seqz.baf2$Bf)

        #Merge the two samples and remove NAs in the common column
        BAF=merge(BAF1,BAF2[,-1],by = "pos",all=TRUE)
        BAF$chrom=ichrom
        BAF=BAF[,c(2,1,3,4)]

        #Get arms
        arms=copynumber:::getArms(BAF$chrom,BAF$pos,cyto.data=getAnywhere(assembly)$objs[[1]]) 
        pints=which(arms=="p")
        qints=which(arms=="q")
        
        ####ATTENTION: runAscat does not work properly with only one chromosome arm, to solve this, I am making a dummy arm when we have data only for one
        if(length(pints)==0 || length(qints)==0)
        {
        	dummyarm=if(length(pints)==0) "p" else "q"
        	maxpos=max(logR$pos)
        	templogR=logR
        	templogR$pos=templogR$pos+maxpos
        	logR=rbind(logR,templogR)
        	tempBAF=BAF
        	tempBAF$pos=tempBAF$pos+maxpos
        	BAF=rbind(BAF,tempBAF)
        	arms=c(arms,rep(dummyarm,length(arms)))
        	pints=which(arms=="p")
        	qints=which(arms=="q")
        }
        
        chrlist=list(pints,qints)
        chrlist=Filter(length,chrlist)

        #ASCAT data format
        thisdataASCAT=list(chr=chrlist,samples=c("s1","s2"),Tumor_LogR=logR[,c(3,4)],Tumor_BAF=BAF[,c(3,4)],sexchromosomes=sexchromosomes,gender=gender,SNPpos=data.frame(chrs=arms,pos=BAF$pos))

        #We are already using only heterozygous loci, so none of them are germline
        gg=list(germlinegenotypes=rep(FALSE,nrow(BAF)))

        ##We use a full run of ascat to get the combined segmentation from asmultipcf (this may not be ideal)
        thisdataASCAT <- ascat.asmultipcf(thisdataASCAT, ascat.gg = gg, penalty = penalty, wsample=NULL,
                       selectAlg="exact",refine=refine)
        
        thisdataASCAT$chrs=unique(arms)

        resASCAT=ascat.runAscat(thisdataASCAT)
        allele.seg=resASCAT$segments
        colnames(allele.seg)=c("sample","arm", "start.pos", "end.pos", "nMajor", "nMinor")
        allele.seg$chrom=ichrom
        allele.seg1=allele.seg[allele.seg$sample=="s1" & allele.seg$arm!=dummyarm,] 
        allele.seg2=allele.seg[allele.seg$sample=="s2"& allele.seg$arm!=dummyarm,]
        
        ##DEBUG
        write.table(x=allele.seg1,file=paste0("s1ASCAT_",ichrom,".txt"))
        write.table(x=allele.seg2,file=paste0("s2ASCAT_",ichrom,".txt"))


    }
    else if (seg.algo == "multipcf") {
        logR1.wins <- copynumber::winsorize(logR1, verbose = verbose)
        logR2.wins <- copynumber::winsorize(logR2, verbose = verbose)
        logR.wins=merge(logR1.wins,logR2.wins[,-1],by = "pos",all=FALSE) ##multipcf does not allow missing data. Accordingly, I am getting the intersection
        logR.wins=logR.wins[,c(2,1,3,4)]
        
        if (nrow(logR.wins) > 400) { ##With a small number of probes we can use the exact algorithm
 
            allele.seg <- copynumber::multipcf(data = logR.wins, verbose = verbose, 
                                        gamma = gamma, ...)
        } else {
            allele.seg <- copynumber::multipcf(data = logR.wins, verbose = verbose,
                                        gamma = gamma, fast=FALSE, ...)
        }

        ##multipcf allways generates the same breakpoints for the two
        allele.seg1=allele.seg
        allele.seg2=allele.seg
    }
    else {
        stop("Supported segmentation algorithms are only 'asmultipcf' or 'multipcf' from the ASCAT and copynumber package respectively.")
    }

    if (length(grep("chr", seqz.baf1$chromosome)) > 0) {
        allele.seg1$chrom <- paste("chr", allele.seg1$chrom, sep = "")
    }
    breaks1 <- allele.seg1[, c("chrom", "start.pos", "end.pos", 
        "arm")]
    not.uniq <- which(breaks1$end.pos == c(breaks1$start.pos[-1], 
        0))
    breaks1$end.pos[not.uniq] <- breaks1$end.pos[not.uniq] - 1

    if (length(grep("chr", seqz.baf2$chromosome)) > 0) {
        allele.seg2$chrom <- paste("chr", allele.seg2$chrom, sep = "")
    }
    breaks2 <- allele.seg2[, c("chrom", "start.pos", "end.pos", 
        "arm")]
    not.uniq <- which(breaks2$end.pos == c(breaks2$start.pos[-1], 
        0))
    breaks2$end.pos[not.uniq] <- breaks2$end.pos[not.uniq] - 1
    
    return(list(s1=breaks1,s2=breaks2))
}

##Function separated from sequenza-extract for clarity ##There is some problems here
merge.breaks <- function(breaks, breaks.het) {
    merged.breaks <- unique(sort(c(breaks$start.pos, 
                         breaks$end.pos, breaks.het$start.pos, 
                         breaks.het$end.pos)))
    merged.breaks <- merged.breaks[diff(merged.breaks) > 1]
    merged.start <- merged.breaks
    merged.start[-1] <- merged.start[-1] + 1
    breaks <- data.frame(chrom = unique(breaks$chrom), 
               start.pos = merged.start[-(length(merged.start))], 
               end.pos = merged.breaks[-1])
}

merge.breaks.byarm = function(breaks, breaks.het) {
    chr.p <- merge.breaks(breaks[breaks$arm == "p", ], breaks.het[breaks.het$arm == "p", ])
    chr.q <- merge.breaks(breaks[breaks$arm == "q", ], breaks.het[breaks.het$arm == "q", ])
    breaks <- rbind(chr.p, chr.q)
}

merge.multibreaks = function(breaks1, breaks2)
{
    return(list(s1=merge.breaks.byarm(breaks1$s1,breaks2$s1),s2=merge.breaks.byarm(breaks1$s2,breaks2$s2)))
}

##Modified from sequenza.extract, to use asmultipcf and multipcf. I may add them back in the future. It requires ASCAT >=2.5
my.sequenza.extract.paired=function (file1, file2, gz = TRUE, window = 1e+06, overlap = 1, gamma = 80, penalty=25,
          kmin = 10, gamma.pcf = 140, kmin.pcf = 40, mufreq.treshold = 0.1, 
          min.reads = 40, min.reads.normal = 10, min.reads.baf = 1, 
          max.mut.types = 1, min.type.freq = 0.9, min.fw.freq = 0, 
          verbose = TRUE, chromosome.list = NULL, breaks.method = "het", breaks= NULL,
          assembly = "hg19", weighted.mean = TRUE, normalization.method = "mean", 
          gc.stats1 = NULL, gc.stats2 = NULL, gender=c("XX","XX"), sexchromosomes=c("X","Y")) 
{

    ####TODO:
    ###### I should take into consideration the fact that we may not have data for a specific DNA stretch for one sample while we do for the other. So far I am not considering this. Doing it at a nucletoid level makes no sense either
     
    require(sequenza)
    require(ASCAT)
 
    ##Get GCstats
    if (is.null(gc.stats1)) {
        gc.stats1 <- gc.sample.stats(file1, gz = gz) #get precalculated gc data
    }
    chr.vect1 <- as.character(gc.stats1$file.metrics$chr)

    if (is.null(gc.stats2)) {
        gc.stats2 <- gc.sample.stats(file2, gz = gz) #get precalculated gc data
    }
    chr.vect2 <- as.character(gc.stats2$file.metrics$chr)
  
    chr.vect=intersect(chr.vect1,chr.vect2) ##Intersection, we only want chromosomes in common
 
    ##Use given breaks if given 
    if (is.null(dim(breaks))) {
        breaks.all <- NULL
    } else {
        breaks.all <- breaks
    }
 
    ##Get normalization values   
    if (normalization.method != "mean") {
        gc.vect1 <- setNames(gc.stats1$raw.median, gc.stats1$gc.values) ##gc normalization values
        gc.vect2 <- setNames(gc.stats2$raw.median, gc.stats2$gc.values) ##gc normalization values
    } else {
        gc.vect1 <- setNames(gc.stats1$raw.mean, gc.stats1$gc.values)
        gc.vect2 <- setNames(gc.stats2$raw.mean, gc.stats2$gc.values)
    }
   
    ##Initiate variables 
    windows.baf1 <- list()
    windows.ratio1 <- list()
    mutation.list1 <- list()
    segments.list1 <- list()
    coverage.list1 <- list()
    windows.baf2 <- list()
    windows.ratio2 <- list()
    mutation.list2 <- list()
    segments.list2 <- list()
    coverage.list2 <- list()

    ##This is common
    if (is.null(chromosome.list)) { ##Final list of chromosomes
        chromosome.list <- chr.vect
    } else {
        chromosome.list <- chromosome.list[chromosome.list %in% chr.vect]
    }
  
  
    for (chr in chromosome.list) { ##Analysis per chr
        if (verbose) {
            message("Processing ", chr, ": ", appendLF = FALSE)
        }
    
        data1=prepDataChr(gc.stats1,file1,chr,chr.vect1,gc.vect1,min.reads.baf,window,overlap,gz)
        data2=prepDataChr(gc.stats2,file2,chr,chr.vect2,gc.vect2,min.reads.baf,window,overla,gz)
  
        if (nrow(data1$seqz.het) > 0 && nrow(data2$seqz.het) > 0) ##Het post for the two, we can use asmultipcf in het and/or multipcf in hom/all
        {
            breaks.method.i <- breaks.method
            
            if (is.null(breaks.all)) {
##TODO: Indicate that here we are not applyting kmin for anything, gamma only for multipcf, and are applying a penalty for asmultipcf
                if (breaks.method.i == "full") { #Uses both homozygous and heterozygous positions. The first using pcf on the depth, the second aspcf.
                    multibreaks <- find.breaks.multi(data1$seqz.data, data2$seqz.data, penalty = penalty, assembly = assembly, seg.algo = "multipcf") #multipcf
                    multibreaks.het <- try(find.breaks.multi(data1$seqz.het, data2$seqz.het, penalty = penalty, gamma = gamma, assembly = assembly ,seg.algo = "asmultipcf"), silent = FALSE) #asmultipcf
                
                    if (!is.null(breaks.het)) {# merge homozygous and heterozygous breaks if the second worked
                        multibreaks=merge.multibreaks(multibreaks,multibreaks.het)                    
                    }
                }   
                else if (breaks.method.i == "het") { #only uses BAF
                    multibreaks <- try(find.breaks.multi(data1$seqz.het, data2$seqz.het, gamma = gamma, penalty = penalty, assembly = assembly, seg.algo = "asmultipcf"), silent = FALSE)
                }
##Commented out since I am in a rush right now. I should reimplement this in the future
#            else if (breaks.method.i == "fast") { #method fast uses the binned windows for BAF
#                BAF <- data.frame(chrom = chr, pos = c(seqz.b.win[[1]]$start, 
#                                                     tail(seqz.b.win[[1]]$end, n = 1)), s1 = c(seqz.b.win[[1]]$mean, 
#                                                     tail(seqz.b.win[[1]]$mean, n = 1)))
#                BAF$s1[is.na(BAF$s1)] <- 0
#                logR <- data.frame(chrom = chr, pos = c(seqz.r.win[[1]]$start, 
#                                                      tail(seqz.r.win[[1]]$end, n = 1)), s1 = c(log2(seqz.r.win[[1]]$mean), 
#                                                      log2(tail(seqz.r.win[[1]]$mean, n = 1))))
#                not.cover <- is.na(logR$s1)
#                BAF <- BAF[!not.cover, ]
#                logR <- logR[!not.cover, ]
#                logR.wins <- copynumber::winsorize(logR, verbose = FALSE)
#                allele.seg <- copynumber::aspcf(logR = logR.wins, 
#                                                BAF = BAF, baf.thres = c(0, 0.5), verbose = FALSE, 
#                                                gamma = gamma, kmin = kmin)
#                if (length(grep("chr", chr)) > 0) {
#                    allele.seg$chrom <- paste("chr", allele.seg$chrom, 
#                                                sep = "")
#                }
#              
#                breaks <- allele.seg[, c("chrom", "start.pos", 
#                                       "end.pos")]
#                not.uniq <- which(breaks$end.pos == c(breaks$start.pos[-1], 
#                                                    0))
#                breaks$end.pos[not.uniq] <- breaks$end.pos[not.uniq] - 
#                1
#            }
                else {
                    #stop("The implemented segmentation methods are 'full', 'het' and 'fast'.")
                    stop("The implemented segmentation methods are 'full' and 'het'.")
                }
           
##working with multibreaks now
 
                ##Adds information to the new segments, this is for the sequenza model
                if (!is.null(multibreaks$s1) & nrow(multibreaks$s1) > 0) {
                    seg.s1 <- segment.breaks(seqz.tab = data1$seqz.data, 
                                     breaks = multibreaks$s1, min.reads.baf = min.reads.baf, 
                                     weighted.mean = weighted.mean)
                } else { #No breaks, segments=full chromosomes
                    seg.s1 <- segment.breaks(data1$seqz.data, breaks = data.frame(chrom = chr, 
                                                                    start.pos = min(data1$seqz.data$position, na.rm = TRUE), 
                                                                    end.pos = max(data1$seqz.data$position, na.rm = TRUE)), 
                                     weighted.mean = weighted.mean)
                }
                
                ##Adds information to the new segments, this is for the sequenza model
                if (!is.null(multibreaks$s2) & nrow(multibreaks$s2) > 0) {
                    seg.s2 <- segment.breaks(seqz.tab = data2$seqz.data, 
                                     breaks = multibreaks$s2, min.reads.baf = min.reads.baf, 
                                     weighted.mean = weighted.mean)
                } else { #No breaks, segments=full chromosomes
                    seg.s2 <- segment.breaks(data2$seqz.data, breaks = data.frame(chrom = chr, 
                                                                    start.pos = min(data2$seqz.data$position, na.rm = TRUE), 
                                                                    end.pos = max(data2$seqz.data$position, na.rm = TRUE)), 
                                     weighted.mean = weighted.mean)
                }
            }
    
        } else { #There is no Het information for at least one of them, so we can only do multipcf
    
            if (breaks.method == "full") { ##Then full = pcf only
                multibreaks <- find.breaks.multi(data1$seqz.data, data2$seqz.data, penalty = penalty, assembly = assembly, seg.algo = "multipcf") #multipcf
            }
            else { ## and the rest, whole chromosomes
                breaks1 = data.frame(chrom = chr, start.pos = min(data1$seqz.data$position, 
                                                             na.rm = TRUE), end.pos = max(data1$seqz.data$position, 
                                                             na.rm = TRUE))
                breaks2 = data.frame(chrom = chr, start.pos = min(data2$seqz.data$position, 
                                                             na.rm = TRUE), end.pos = max(data2$seqz.data$position, 
                                                             na.rm = TRUE))
                multibreaks=list(s1=breaks1,s2=breaks2)
            }

            seg.s1 <- segment.breaks(data1$seqz.data, breaks = multibreaks$s1, 
                                   weighted.mean = weighted.mean)
            seg.s2 <- segment.breaks(data2$seqz.data, breaks = multibreaks$s2, 
                                   weighted.mean = weighted.mean)
        }
        
        mut.tab1 <- mutation.table(data1$seqz.data, mufreq.treshold = mufreq.treshold, 
                                  min.reads = min.reads, min.reads.normal = min.reads.normal, 
                                  max.mut.types = max.mut.types, min.type.freq = min.type.freq, 
                                  min.fw.freq = min.fw.freq, segments = seg.s1)
        windows.baf1[[which(chromosome.list == chr)]] <- data1$seqz.b.win[[1]]
        windows.ratio1[[which(chromosome.list == chr)]] <- data1$seqz.r.win[[1]]
        mutation.list1[[which(chromosome.list == chr)]] <- mut.tab1
        segments.list1[[which(chromosome.list == chr)]] <- seg.s1
        coverage.list1[[which(chromosome.list == chr)]] <- data.frame(sum = sum(as.numeric(data1$seqz.data$depth.tumor), 
                                                                               na.rm = TRUE), N = length(data1$seqz.data$depth.tumor))
        
        mut.tab2 <- mutation.table(data2$seqz.data, mufreq.treshold = mufreq.treshold, 
                                  min.reads = min.reads, min.reads.normal = min.reads.normal, 
                                  max.mut.types = max.mut.types, min.type.freq = min.type.freq, 
                                  min.fw.freq = min.fw.freq, segments = seg.s2)
        windows.baf2[[which(chromosome.list == chr)]] <- data2$seqz.b.win[[1]]
        windows.ratio2[[which(chromosome.list == chr)]] <- data2$seqz.r.win[[1]]
        mutation.list2[[which(chromosome.list == chr)]] <- mut.tab2
        segments.list2[[which(chromosome.list == chr)]] <- seg.s2
        coverage.list2[[which(chromosome.list == chr)]] <- data.frame(sum = sum(as.numeric(data2$seqz.data$depth.tumor), 
                                                                               na.rm = TRUE), N = length(data2$seqz.data$depth.tumor))
        if (verbose) {
          message("A:", nrow(mut.tab1), " variant calls; ", nrow(data1$seqz.het), 
                  " heterozygous positions; ", sum(data1$seqz.hom), " homozygous positions.")
          message("B:", nrow(mut.tab2), " variant calls; ", nrow(data2$seqz.het), 
                  " heterozygous positions; ", sum(data2$seqz.hom), " homozygous positions.")
        }
    }
    
    names(windows.baf1) <- chromosome.list
    names(windows.ratio1) <- chromosome.list
    names(mutation.list1) <- chromosome.list
    names(segments.list1) <- chromosome.list
    
    coverage.list1 <- do.call(rbind, coverage.list1)
    coverage1 <- sum(coverage.list1$sum)/sum(coverage.list1$N)
    
    names(windows.baf2) <- chromosome.list
    names(windows.ratio2) <- chromosome.list
    names(mutation.list2) <- chromosome.list
    names(segments.list2) <- chromosome.list
    
    coverage.list2 <- do.call(rbind, coverage.list2)
    coverage2 <- sum(coverage.list2$sum)/sum(coverage.list2$N)

    return(list(A=list(BAF = windows.baf1, ratio = windows.ratio1, mutations = mutation.list1, 
              segments = segments.list1, chromosomes = chromosome.list, 
              gc = gc.stats1, avg.depth1 = round(coverage1, 0)),B=list(BAF = windows.baf2, ratio = windows.ratio2, mutations = mutation.list2, 
              segments = segments.list2, chromosomes = chromosome.list, 
              gc = gc.stats2, avg.depth = round(coverage2, 0))))
}

##WARNING: STATISTICS ARE JUST COPIED WHEN WE SPLIT FRAGMENTS, ONLY GENOTYPES AND BREAKPOINTS SHOULD BE USED
syncsegments= function(segs1,segs2) {

	chrs=unique(c(as.character(segs1$chromosome),as.character(segs2$chromosome)))
	breakpoints=sort(unique(c(segs1$start.pos,segs1$end.pos,segs2$start.pos,segs2$end.pos)))
	maxnsegments=length(breakpoints)-1-(length(chrs)-1)
	outsegments1=data.frame(chromosome=character(maxnsegments),start.pos=numeric(maxnsegments),end.pos=numeric(maxnsegments),Bf=numeric(maxnsegments),N.BAF=numeric(maxnsegments),sd.BAF=numeric(maxnsegments),depth.ratio=numeric(maxnsegments),N.ratio=numeric(maxnsegments),sd.ratio=numeric(maxnsegments),CNt=numeric(maxnsegments),A=numeric(maxnsegments),B=numeric(maxnsegments),LPP=numeric(maxnsegments),stringsAsFactors=FALSE)
	outsegments2=data.frame(chromosome=character(maxnsegments),start.pos=numeric(maxnsegments),end.pos=numeric(maxnsegments),Bf=numeric(maxnsegments),N.BAF=numeric(maxnsegments),sd.BAF=numeric(maxnsegments),depth.ratio=numeric(maxnsegments),N.ratio=numeric(maxnsegments),sd.ratio=numeric(maxnsegments),CNt=numeric(maxnsegments),A=numeric(maxnsegments),B=numeric(maxnsegments),LPP=numeric(maxnsegments),stringsAsFactors=FALSE)
	
	iseg=1
	
	for (nchr in 1:length(chrs)) {
		chr=chrs[nchr]
		segs1i=segs1[segs1$chromosome==chr,]
		segs2i=segs2[segs2$chromosome==chr,]
		breakpoints=sort(unique(c(segs1i$start.pos,segs1i$end.pos,segs2i$start.pos,segs2i$end.pos)))
				
		for (nbr in 2:length(breakpoints)) {	
			b1=breakpoints[nbr-1]
			b2=breakpoints[nbr]
			
			#print(paste0("Breakpoints ",b1," and",b2))
			nseg1=which(segs1i$start.pos<=b1 & segs1i$end.pos>=b2)
			nseg2=which(segs2i$start.pos<=b1 & segs2i$end.pos>=b2)
			
			if (length(nseg1)==1 & length(nseg2)==1) {
				#chromosome start.pos   end.pos        Bf N.BAF     sd.BAF depth.ratio N.ratio  sd.ratio CNt A B       LPP
				outsegments1[iseg,]=c(chr,b1,b2,segs1i[nseg1,4:ncol(outsegments1)])
				outsegments2[iseg,]=c(chr,b1,b2,segs2i[nseg2,4:ncol(outsegments2)])
				iseg=iseg+1
			} else if (length(nseg1)>1 | length(nseg2)>1) {
				stop(paste0("More than one compatible breakpoint for ",b1," and ",b2))
			}
		}
	}
	return(list(segmentsA=outsegments1[1:iseg-1,],segmentsB=outsegments2[1:iseg-1,]))
}

getgenotypesimilarity= function(segments) {
	segsA=segments$segmentsA[,c("A","B")]
	segsB=segments$segmentsB[,c("A","B")]
	##The comparison retuns a allele by allele comparison, the apply allows me to see if both are true, then divide by the number of rows
	sum(apply(MARGIN=1,X=(segsA==segsB),FUN=function(x){sum(x)==length(x)}))/nrow(segsA)
}

################################################

#Main
#####
options = commandArgs(trailingOnly = TRUE)

name1=options[1]
name2=options[2]
outdir=options[3]

id1=sub(pattern=".seqz.gz",replacement="",x=name1)
id1=basename(id1)
id2=sub(pattern=".seqz.gz",replacement="",x=name2)
id2=basename(id2)

dir.create(outdir)
ncores=as.numeric(options[4])
patient=options[5]

chromosomes=c(seq(1,22),"X")

#min.reads=20,min.reads.normal=10,min.reads.baf=20,mufreq.treshold=0.25,max.mut.types=3
 
##TODO: ADD filtering parameters. I gotta think better about them

datacomp=my.sequenza.extract.paired(file1=name1, file2=name2, chromosome.list=chromosomes,mufreq.treshold=0.10, min.reads=20,min.reads.normal=10, max.mut.types =3) #gamma and kmin not being used, penalty not being changed.
##mufreq.treshold default
##minreads lowered to 20, we don't have so much coverage to put 40
##minreads normal default
##minreads baf default (1, all reads are used for the segment average, this is not the BAF per site)
##method default (baf)
##max mut types relaxed to 3
##min.type.freq default, I do not fully understand it
##min.fw.freq default. We are not considering differences between strands


##Pending, new funtion to print them using a par to compare chr by chr the two samples
printsegmentation(datacomp$A,outdir=paste0(outdir,"/A_chr"),chromosomes=chromosomes)
printsegmentation(datacomp$B,outdir=paste0(outdir,"/B_chr"),chromosomes=chromosomes)

cpA <- sequenza.fit(datacomp$A,mc.cores=ncores,chromosome.list=chromosomes,female=TRUE,N.BAF.filter=5)
cpB <- sequenza.fit(datacomp$B,mc.cores=ncores,chromosome.list=chromosomes,female=TRUE,N.BAF.filter=5)
##CNt max, 20 copies in a segment
##N.ratio.filter default, a segment needs 10 positions to be included in the ratio
##N.BAF.filter  positions in the segment to be used for BAF. Since just one read is enough, I will use 5 
##segment.filter 1e6bp, default
##ratio.priority false. The ratio does not have priority over BAF (with priority, CNt is determined using ratio, then with BAF the different alleles)

sequenza.results(datacomp$A,cp.table=cpA,sample.id=id1,out.dir=outdir,chromosome.list=chromosomes,female=TRUE)
sequenza.results(datacomp$B,cp.table=cpB,sample.id=id2,out.dir=outdir,chromosome.list=chromosomes,female=TRUE)

segmentsA=read.table(paste0(outdir,"/",id1,"_segments.txt"),header=TRUE)
segmentsB=read.table(paste0(outdir,"/",id2,"_segments.txt"),header=TRUE)

finalsegments=syncsegments(segmentsA,segmentsB)

propgenotypes=getgenotypesimilarity(finalsegments)

results=data.frame(patient=patient,prop=propgenotypes,nseg=nrow(finalsegments$segmentsA),cellA=get.ci(cpA)$max.cellularity,cellB=get.ci(cpB)$max.cellularity,ploidyA=get.ci(cpA)$max.ploidy,ploidyB=get.ci(cpB)$max.ploidy)

write.csv(x=results,file="results.csv",quote = FALSE,row.names = FALSE)
