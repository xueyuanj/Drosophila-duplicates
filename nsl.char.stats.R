current= setwd("/data/nSL/result")

labels <- read.table("melDuplicates.eucDist.classified")
neocid <- labels[labels$V9=="neochild",1]
concid <- labels[labels$V9=="cons",1]
specid <- labels[labels$V9=="spec",1]

charstats <- read.table("nsl.child.005.r549.anno.norm.charstats")
singles <- read.table("nsl.single.005.r549.anno.norm.charstats")

neoc <- charstats[charstats$V1%in%neocid,]
con <- charstats[charstats$V1%in%concid,]
spec <- charstats[charstats$V1%in%specid,]

boxplot(neoc$V2, con$V2, spec$V2, singles$V2, outline = FALSE, notch = TRUE, names = c("neoc", "con", "spec", "single"), main="nSL score(highest)", col = c("mintcream", "lightcyan", "lavender", "lightgray"))

boxplot(neoc$V3, con$V3, spec$V3, singles$V3, outline = FALSE, notch = TRUE, names = c("neoc", "con", "spec", "single"), main="nSL score(lowest)", col = c("mintcream", "lightcyan", "lavender", "lightgray"))

boxplot(neoc$V4, con$V4, spec$V4, singles$V4, outline = FALSE, notch = TRUE, names = c("neoc", "con", "spec", "single"), main="nSL score(median)", col = c("mintcream", "lightcyan", "lavender", "lightgray"))



############################################
######## All of them
############################################
allchild<- read.table("child.all.005.r549.anno.nsl.norm")
alls<-read.table("single.all.005.r549.anno.nsl.norm")

allneo <- allchild[allchild$V1%in%neocid,3]
allcon <- allchild[allchild$V1%in%concid,3]
allspe <- allchild[allchild$V1%in%specid,3]

boxplot( abs(allneo),abs(allcon),abs(allspe), abs(alls$V3), outline = FALSE, notch = TRUE, names = c("neoc", "con", "spec", "single"), main="nSL score(all, abs)", col = c("mintcream", "lightcyan", "lavender", "lightgray"))

wilcox.test(allneo,allcon)
#P-value:                         0.004514279
wilcox.test(allneo,allspe)
#P-value:                         0.7646048
wilcox.test(allneo,alls$V3)
#P-value:                         0.0002136397
wilcox.test(allcon, allspe)
#P-value:                         0.01353568
wilcox.test(allcon, alls$V3)
#P-value:                         0.1395031
wilcox.test(allspe, alls$V3)
#P-value:                         0.05646592


charstats <- read.table("nsl.child.005.r549.anno.norm.charstats.abs")
singles <- read.table("nsl.single.005.r549.anno.norm.charstats.abs")

neoc <- charstats[charstats$V1%in%neocid,]
con <- charstats[charstats$V1%in%concid,]
spec <- charstats[charstats$V1%in%specid,]

boxplot(neoc$V2, con$V2, spec$V2, singles$V2, outline = FALSE, notch = TRUE, names = c("neoc", "con", "spec", "single"), main="nSL score(highest,abs)", col = c("mintcream", "lightcyan", "lavender", "lightgray"))

boxplot(neoc$V3, con$V3, spec$V3, singles$V3, outline = FALSE, notch = TRUE, names = c("neoc", "con", "spec", "single"), main="nSL score(lowest, abs)", col = c("mintcream", "lightcyan", "lavender", "lightgray"))

boxplot(neoc$V4, con$V4, spec$V4, singles$V4, outline = FALSE, notch = TRUE, names = c("neoc", "con", "spec", "single"), main="nSL score(median, abs)", col = c("mintcream", "lightcyan", "lavender", "lightgray"))


##########################################################
################ top 5% and top 1%########################
higest.c <- charstats[,1:2]
lowest.c <- charstats[,c(1,3)]
med.c <- charstats[,c(1,4)]

h.top5 <- head(higest.c[order(higest.c$V2, decreasing = T), 1], 100)
h.top1 <- head(higest.c[order(higest.c$V2, decreasing = T), 1], 1)


l.top5 <- head(lowest.c[order(lowest.c$V3, decreasing = T), 1], 0.05*nrow(lowest.c))
l.top1 <- head(lowest.c[order(lowest.c$V3, decreasing = T), 1], 1)

m.top5 <- head(med.c[order(med.c$V4, decreasing = T), 1], 0.05*nrow(med.c))
m.top1 <- head(med.c[order(med.c$V4, decreasing = T), 1], 1)

labels[labels$V1%in%h.top5, c(1,9)]
labels[labels$V1%in%l.top5, c(1,9)]
labels[labels$V1%in%m.top5, c(1,9)]

labels[labels$V1%in%h.top1, c(1,9)]
labels[labels$V1%in%l.top1, c(1,9)]
labels[labels$V1%in%m.top1, c(1,9)]





current.1= setwd("/data/nSL/standnsl")

allgenes <- read.table("nsl.allgene.005.r549.anno.norm.charstats.abs")
all.h <- allgenes[,1:2]
all.l <- allgenes[,c(1,3)]
all.m <- allgenes[,c(1,4)]

all.h.top5 <- head(all.h[order(all.h$V2, decreasing = T), ], 0.05*nrow(all.h))
all.h.top1 <- head(all.h[order(all.h$V2, decreasing = T), ], 0.01*nrow(all.h))

all.l.top5 <- head(all.l[order(all.l$V3, decreasing = T), ], 0.05*nrow(all.l))
all.l.top1 <- head(all.l[order(all.l$V3, decreasing = T), ], 0.01*nrow(all.l))

all.m.top5 <- head(all.m[order(all.m$V4, decreasing = T), ], 0.05*nrow(all.m))
all.m.top1 <- head(all.m[order(all.m$V4, decreasing = T), ], 0.01*nrow(all.m))

labels[labels$V1%in%all.h.top5$V1, c(1,9)]
labels[labels$V1%in%all.h.top1$V1, c(1,9)]

labels[labels$V1%in%all.l.top5$V1, c(1,9)]
labels[labels$V1%in%all.l.top1$V1, c(1,9)]

labels[labels$V1%in%all.m.top5$V1, c(1,9)]
labels[labels$V1%in%all.m.top1$V1, c(1,9)]


all.h.top5[all.h.top5$V1%in%labels$V1,]



################################################################
# individual gene on different chromosome arms
################################################################

all.2l <- read.table("nsl.2L.005.r549.anno.norm.charstats.abs")
all.2r <- read.table("nsl.2R.005.r549.anno.norm.charstats.abs")
all.3l <- read.table("nsl.3L.005.r549.anno.norm.charstats.abs")
all.3r <- read.table("nsl.3R.005.r549.anno.norm.charstats.abs")
all.x <- read.table("nsl.X.005.r549.anno.norm.charstats.abs")


all.2l.h <- all.2l[,1:2]
all.2l.l <- all.2l[, c(1,3)]
all.2l.m <- all.2l[, c(1,4)]
all.2r.h <- all.2r[,1:2]
all.2r.l <- all.2r[, c(1,3)]
all.2r.m <- all.2r[, c(1,4)]
all.3l.h <- all.3l[,1:2]
all.3l.l <- all.3l[, c(1,3)]
all.3l.m <- all.3l[,1:4]
all.3r.h <- all.3r[, c(1,2)]
all.3r.l <- all.3r[, c(1,3)]
all.3r.m <- all.3r[,c(1,4)]
all.x.h <- all.x[, c(1,2)]
all.x.l <- all.x[, c(1,3)]
all.x.m <- all.x[, c(1,4)]

sink("nsl.chromarm.top5&top1")
all.2l.h.top5 <- head(all.2l.h[order(all.2l.h$V2, decreasing = T), ], 0.05*nrow(all.2l.h))
all.2l.h.top1 <- head(all.2l.h[order(all.2l.h$V2, decreasing = T), ], 0.01*nrow(all.2l.h))
labels[labels$V1%in%all.2l.h.top5$V1, c(1,9)]
labels[labels$V1%in%all.2l.h.top1$V1, c(1,9)]

all.2l.l.top5 <- head(all.2l.l[order(all.2l.l$V3, decreasing = T), ], 0.05*nrow(all.2l.l))
all.2l.l.top1 <- head(all.2l.l[order(all.2l.l$V3, decreasing = T), ], 0.01*nrow(all.2l.l))
labels[labels$V1%in%all.2l.l.top5$V1, c(1,9)]
labels[labels$V1%in%all.2l.l.top1$V1, c(1,9)]

all.2l.m.top5 <- head(all.2l.m[order(all.2l.m$V4, decreasing = T), ], 0.05*nrow(all.2l.m))
all.2l.m.top1 <- head(all.2l.m[order(all.2l.m$V4, decreasing = T), ], 0.01*nrow(all.2l.m))
labels[labels$V1%in%all.2l.m.top5$V1, c(1,9)]
labels[labels$V1%in%all.2l.m.top1$V1, c(1,9)]




all.2r.h.top5 <- head(all.2r.h[order(all.2r.h$V2, decreasing = T), ], 0.05*nrow(all.2r.h))
all.2r.h.top1 <- head(all.2r.h[order(all.2r.h$V2, decreasing = T), ], 0.01*nrow(all.2r.h))
labels[labels$V1%in%all.2r.h.top5$V1, c(1,9)]
labels[labels$V1%in%all.2r.h.top1$V1, c(1,9)]

all.2r.l.top5 <- head(all.2r.l[order(all.2r.l$V3, decreasing = T), ], 0.05*nrow(all.2r.l))
all.2r.l.top1 <- head(all.2r.l[order(all.2r.l$V3, decreasing = T), ], 0.01*nrow(all.2r.l))
labels[labels$V1%in%all.2r.l.top5$V1, c(1,9)]
labels[labels$V1%in%all.2r.l.top1$V1, c(1,9)]

all.2r.m.top5 <- head(all.2r.m[order(all.2r.m$V4, decreasing = T), ], 0.05*nrow(all.2r.m))
all.2r.m.top1 <- head(all.2r.m[order(all.2r.m$V4, decreasing = T), ], 0.01*nrow(all.2r.m))
labels[labels$V1%in%all.2r.m.top5$V1, c(1,9)]
labels[labels$V1%in%all.2r.m.top1$V1, c(1,9)]




all.3l.h.top5 <- head(all.3l.h[order(all.3l.h$V2, decreasing = T), ], 0.05*nrow(all.3l.h))
all.3l.h.top1 <- head(all.3l.h[order(all.3l.h$V2, decreasing = T), ], 0.01*nrow(all.3l.h))
labels[labels$V1%in%all.3l.h.top5$V1, c(1,9)]
labels[labels$V1%in%all.3l.h.top1$V1, c(1,9)]

all.3l.l.top5 <- head(all.3l.l[order(all.3l.l$V3, decreasing = T), ], 0.05*nrow(all.3l.l))
all.3l.l.top1 <- head(all.3l.l[order(all.3l.l$V3, decreasing = T), ], 0.01*nrow(all.3l.l))
labels[labels$V1%in%all.3l.l.top5$V1, c(1,9)]
labels[labels$V1%in%all.3l.l.top1$V1, c(1,9)]

all.3l.m.top5 <- head(all.3l.m[order(all.3l.m$V4, decreasing = T), ], 0.05*nrow(all.3l.m))
all.3l.m.top1 <- head(all.3l.m[order(all.3l.m$V4, decreasing = T), ], 0.01*nrow(all.3l.m))
labels[labels$V1%in%all.3l.m.top5$V1, c(1,9)]
labels[labels$V1%in%all.3l.m.top1$V1, c(1,9)]



all.3r.h.top5 <- head(all.3r.h[order(all.3r.h$V2, decreasing = T), ], 0.05*nrow(all.3r.h))
all.3r.h.top1 <- head(all.3r.h[order(all.3r.h$V2, decreasing = T), ], 0.01*nrow(all.3r.h))
labels[labels$V1%in%all.3r.h.top5$V1, c(1,9)]
labels[labels$V1%in%all.3r.h.top1$V1, c(1,9)]

all.3r.l.top5 <- head(all.3r.l[order(all.3r.l$V3, decreasing = T), ], 0.05*nrow(all.3r.l))
all.3r.l.top1 <- head(all.3r.l[order(all.3r.l$V3, decreasing = T), ], 0.01*nrow(all.3r.l))
labels[labels$V1%in%all.3r.l.top5$V1, c(1,9)]
labels[labels$V1%in%all.3r.l.top1$V1, c(1,9)]

all.3r.m.top5 <- head(all.3r.m[order(all.3r.m$V4, decreasing = T), ], 0.05*nrow(all.3r.m))
all.3r.m.top1 <- head(all.3r.m[order(all.3r.m$V4, decreasing = T), ], 0.01*nrow(all.3r.m))
labels[labels$V1%in%all.3r.m.top5$V1, c(1,9)]
labels[labels$V1%in%all.3r.m.top1$V1, c(1,9)]



all.x.h.top5 <- head(all.x.h[order(all.x.h$V2, decreasing = T), ], 0.05*nrow(all.x.h))
all.x.h.top1 <- head(all.x.h[order(all.x.h$V2, decreasing = T), ], 0.01*nrow(all.x.h))
labels[labels$V1%in%all.x.h.top5$V1, c(1,9)]
labels[labels$V1%in%all.x.h.top1$V1, c(1,9)]

all.x.l.top5 <- head(all.x.l[order(all.x.l$V3, decreasing = T), ], 0.05*nrow(all.x.l))
all.x.l.top1 <- head(all.x.l[order(all.x.l$V3, decreasing = T), ], 0.01*nrow(all.x.l))
labels[labels$V1%in%all.x.l.top5$V1, c(1,9)]
labels[labels$V1%in%all.x.l.top1$V1, c(1,9)]

all.x.m.top5 <- head(all.x.m[order(all.x.m$V4, decreasing = T), ], 0.05*nrow(all.x.m))
all.x.m.top1 <- head(all.x.m[order(all.x.m$V4, decreasing = T), ], 0.01*nrow(all.x.m))
labels[labels$V1%in%all.x.m.top5$V1, c(1,9)]
labels[labels$V1%in%all.x.m.top1$V1, c(1,9)]
sink()


nsl.2l.all <-read.table("dmel.005.2L.r549.anno.norm.tab")
nsl.2r.all <-read.table("dmel.005.2R.r549.anno.norm.tab")
nsl.3l.all <-read.table("dmel.005.3L.r549.anno.norm.tab")
nsl.3r.all <-read.table("dmel.005.3R.r549.anno.norm.tab")
nsl.x.all <-read.table("dmel.005.X.r549.anno.norm.tab")

top5.nsl.2l<-nsl.2l.all[abs(nsl.2l.all$V4)>1.93932, 3]
top5.nsl.2r<-nsl.2r.all[abs(nsl.2r.all$V4)>1.9327, 3]
top5.nsl.3l<-nsl.3l.all[abs(nsl.3l.all$V4)>1.948978, 3]
top5.nsl.3r<-nsl.3r.all[abs(nsl.3r.all$V4)>1.926, 3]
top5.nsl.x<-nsl.x.all[abs(nsl.x.all$V4)>1.96644, 3]

all.x.h.top5 <- head(nsl.x.all[order(abs(nsl.x.all$V4), decreasing = T), ], 0.05*nrow(nsl.x.all))
all.2l.h.top5 <- head(nsl.2l.all[order(abs(nsl.2l.all$V4), decreasing = T), ], 0.05*nrow(nsl.2l.all))
all.2r.h.top5 <- head(nsl.2r.all[order(abs(nsl.2r.all$V4), decreasing = T), ], 0.05*nrow(nsl.2r.all))
all.3l.h.top5 <- head(nsl.3l.all[order(abs(nsl.3l.all$V4), decreasing = T), ], 0.05*nrow(nsl.3l.all))
all.3r.h.top5 <- head(nsl.3r.all[order(abs(nsl.3r.all$V4), decreasing = T), ], 0.05*nrow(nsl.3r.all))

top5.nsl.2l<-nsl.2l.all[abs(nsl.2l.all$V4)>1.96386, 3]
top5.nsl.2r<-nsl.2r.all[abs(nsl.2r.all$V4)>1.936453, 3]
top5.nsl.3l<-nsl.3l.all[abs(nsl.3l.all$V4)>1.948978, 3]
top5.nsl.3r<-nsl.3r.all[abs(nsl.3r.all$V4)>1.939785, 3]
top5.nsl.x<-nsl.x.all[abs(nsl.x.all$V4)>1.979804, 3]


length(unique(top5.nsl.2l))
length(unique(top5.nsl.2r))
length(unique(top5.nsl.3l))
length(unique(top5.nsl.3r))
length(unique(top5.nsl.x))

#################################################################################
######################## The same for iHS  ######################################
current.1 <- setwd("/data/iHS/ihs_table")
ihs.2l <- read.table("dmel.005.2L.gene.norm.tab")
ihs.2r <- read.table("dmel.005.2R.gene.norm.tab")
ihs.3l <- read.table("dmel.005.3L.gene.norm.tab")
ihs.3r <- read.table("dmel.005.3R.gene.norm.tab")
ihs.x <- read.table("dmel.005.X.gene.norm.tab")

ihs.cut.2l <- head(ihs.2l[order(abs(ihs.2l$V4), decreasing = T), ], 0.05*nrow(ihs.2l))
tail(ihs.cut.2l)
#1.962162
ihs.cut.2r <- head(ihs.2r[order(abs(ihs.2r$V4), decreasing = T), ], 0.05*nrow(ihs.2r))
tail(ihs.cut.2r)
#1.979288
ihs.cut.3l <- head(ihs.3l[order(abs(ihs.3l$V4), decreasing = T), ], 0.05*nrow(ihs.3l))
tail(ihs.cut.3l)
#2.053539
ihs.cut.3r <- head(ihs.3r[order(abs(ihs.3r$V4), decreasing = T), ], 0.05*nrow(ihs.3r))
tail(ihs.cut.3r)
#1.948990
ihs.cut.x <- head(ihs.x[order(abs(ihs.x$V4), decreasing = T), ], 0.05*nrow(ihs.x))
tail(ihs.cut.x)
#1.972646

length(unique(ihs.cut.2l$V3))
#176
length(unique(ihs.cut.2r$V3))
#136
length(unique(ihs.cut.3l$V3))
#163
length(unique(ihs.cut.3r$V3))
#202
length(unique(ihs.cut.x$V3))
#141

labels[labels$V1%in%ihs.cut.2l$V3, c(1,9)]
labels[labels$V1%in%ihs.cut.2r$V3, c(1,9)]
labels[labels$V1%in%ihs.cut.3l$V3, c(1,9)]
labels[labels$V1%in%ihs.cut.3r$V3, c(1,9)]
labels[labels$V1%in%ihs.cut.x$V3, c(1,9)]

currunt.3 <- setwd("/data/HKA/M_script/introns")
nest.id <- read.table("mel.nest.host.id")
my.nest.list <- read.table("dmel.all.intron.to.be.removed.within")
labels[labels$V1%in%nest.id$V1 , c(1,9)]
labels[labels$V1%in%my.nest.list$V1 , c(1,9)]
