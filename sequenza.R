library(sequenza)
options = commandArgs(trailingOnly = TRUE)

name=options[1]
outdir=options[2]
id=sub(pattern=".seqz.*$",replacement="",x=name)
id=basename(id)

dir.create(outdir,showWarnings=FALSE)
ncores=as.numeric(options[3])

chromosomes=c(seq(1,22),"X")
mydata=sequenza.extract(file=name,chromosome.list=chromosomes,min.reads=20,min.reads.normal=10,min.reads.baf=20,mufreq.treshold=0.25,max.mut.types=3) ##I am excluding the GL "chromosomes" and the mitocondria since they generate an error in sequenza, and the Y, since these samples are all from females

for (chr in 1:length(chromosomes)) {
print(chromosomes[chr])
chromosome=chromosomes[chr]
pdf(paste0(outdir,"/chromosome",chromosomes[chr],"_view.pdf"))
try(chromosome.view(mut.tab = mydata$mutations[[chromosome]], baf.windows = mydata$BAF[[chromosome]], ratio.windows = mydata$ratio[[chromosome]], min.N.ratio = 1, segments = mydata$segments[[chromosome]], main = mydata$chromosomes[chromosome]))
dev.off()
}

cp <- sequenza.fit(mydata,mc.cores=ncores,chromosome.list=chromosomes,female=TRUE)
sequenza.results(mydata,cp.table=cp,sample.id=id,out.dir=outdir,chromosome.list=chromosomes,female=TRUE)
