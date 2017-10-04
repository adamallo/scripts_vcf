library(sequenza)
options = commandArgs(trailingOnly = TRUE)

name=options[1]
outdir=options[2]
id=sub(pattern=".seqz.gz",replacement="",x=name)
id=basename(id)

#dir.create(outdir)
ncores=as.numeric(options[3])

mydataplotgc=read.seqz(name) ##This takes like 3 min
gcstats=gc.norm(x = mydataplotgc$depth.ratio, gc = mydataplotgc$GC.percent)
gc.vect <- setNames(gcstats$raw.mean, gcstats$gc.values)
mydataplotgc$adjusted.ratio <- mydataplotgc$depth.ratio / gc.vect[as.character(mydataplotgc$GC.percent)]


pdf(paste0(outdir,"/GCnormalization.pdf"))
par(mfrow = c(1,2), cex = 1, las = 1, bty = 'l')
matplot(gcstats$gc.values, gcstats$raw, type = 'b', col = 1, pch = c(1, 19, 1), lty = c(2, 1, 2), xlab = 'GC content (%)', ylab = 'Uncorrected depth ratio')
legend('topright', legend = colnames(gcstats$raw), pch = c(1, 19, 1))
hist2(mydataplotgc$depth.ratio, mydataplotgc$adjusted.ratio, breaks = prettyLog, key = vkey, panel.first = abline(0, 1, lty = 2), xlab = 'Uncorrected depth ratio', ylab = 'GC-adjusted depth ratio')
dev.off()

chromosomes=c(seq(1,22),"X")
mydata=sequenza.extract(file=name,chromosome.list=chromosomes) ##I am excluding the GL "chromosomes" and the mitocondria since they generate an error in sequenza, and the Y, since these samples are all from females

for (chr in 1:length(chromosomes)) {
print(chromosomes[chr])
chromosome=chromosomes[chr]
pdf(paste0(outdir,"/chromosome",chromosomes[chr],"_view.pdf"))
try(chromosome.view(mut.tab = mydata$mutations[[chromosome]], baf.windows = mydata$BAF[[chromosome]], ratio.windows = mydata$ratio[[chromosome]], min.N.ratio = 1, segments = mydata$segments[[chromosome]], main = mydata$chromosomes[chromosome]))
dev.off()
}

cp <- sequenza.fit(mydata,mc.cores=ncores,chromosome.list=chromosomes,female=TRUE,XY=c("X","Y"))
sequenza.results(mydata,cp.table=cp,sample.id=id,out.dir=outdir,chromosome.list=chromosomes,female=TRUE,XY=c("X","Y"))
