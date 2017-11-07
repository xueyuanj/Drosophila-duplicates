library(EnvStats)
library(BSDA)

display.brewer.all(NULL)
display.brewer.pal(4, "Reds")
brewer.pal(4, "Pastel1")
brewer.pal(4, "Reds")

display.brewer.all(n=NULL)
####################### CDS #############################
setwd("~/")
path2 <- setwd("/data/HKA/M_script/neo_con")
# sink("hka.cds_intron.g4.txt")

cds.s <- read.table("single.w10000s1.r549.cds.hka")
cds.con <- read.table("hka.con.cds.w10000s1.r549.child")
cds.neo <- read.table("hka.neo.cds.w10000s1.r549.child")
cds.spec <- read.table("dmel.r549.spec.w10000s1.cds.hka")


# print("CDS")

# median(cds.neo$V3)
# median(cds.spec$V3)
# median(cds.con$V3)
# median(cds.s$V3)

# twoSamplePermutationTestLocation(cds.neo$V3, cds.spec$V3, fcn="median", alternative="two.sided", n.permutations=1000)

# twoSamplePermutationTestLocation(cds.neo$V3, cds.con$V3, fcn="median", alternative="two.sided", n.permutations=1000)

# twoSamplePermutationTestLocation(cds.neo$V3, cds.s$V3, fcn="median", alternative="two.sided", n.permutations=1000)

# twoSamplePermutationTestLocation(cds.spec$V3, cds.s$V3, fcn="median", alternative="two.sided", n.permutations=1000)

# twoSamplePermutationTestLocation(cds.spec$V3, cds.con$V3, fcn="median", alternative="two.sided", n.permutations=1000)


SIGN.test(cds.neo$V3, md=0, alternative = "two.sided", conf.level = 0.95)
SIGN.test(cds.spec$V3, md=0, alternative = "two.sided", conf.level = 0.95)
SIGN.test(cds.con$V3, md=0, alternative = "two.sided", conf.level = 0.95)
SIGN.test(cds.s$V3, md=0, alternative = "two.sided", conf.level = 0.95)




boxplot(cds.neo$V3, cds.spec$V3, cds.con$V3, cds.s$V3, lwd=2, ylim=range(-60,47),cex.main=3, cex.axis=1.5,  notch = TRUE, outline = FALSE, main="Coding Region", names=c("neofunctionalized","specialized", "conserved", "single"), col=c("#CB181D", "#FB6A4A","#FCAE91", "#FEE5D9"))
boxplot(cds.neo$V3, cds.spec$V3, cds.con$V3, cds.s$V3, lwd=2, ylim=range(-60,47),cex.main=3, cex.axis=1.5,  notch = TRUE, outline = FALSE, main="Coding Region", names=c("neofunctionalized","specialized", "conserved", "single"), col=c("#FBB4AE", "#B3CDE3","#CCEBC5", "#DECBE4"))

par(mgp=c(2.2,1,0))
boxplot(cds.s$V3, cds.neo$V3, cds.spec$V3, cds.con$V3,  boxwex=0.5, lwd=5, ylim=range(-60,47),cex.main=3, cex.axis=1.5,   notch = TRUE, outline = FALSE, main="Coding Regions", names=c("Single-copy", "Neofunctionalized","Specialized", "Conserved" ),border = c("lightgray", "#F03B20", "#FFEDA0", "lightblue"))

boxplot(cds.s$V3, cds.neo$V3, cds.spec$V3, cds.con$V3,  boxwex=0.5, lwd=5, ylim=range(-60,47),cex.main=3, cex.axis=1.5,   notch = TRUE, outline = FALSE, main="Coding Regions", names=c("Single-copy", "Neofunctionalized","Specialized", "Conserved" ),border= c("gray57", "#EF3B2C","#FED976",  "#1D91C0"))

boxplot(cds.s$V3, cds.con$V3,  cds.spec$V3, cds.neo$V3, boxwex=0.5,  ylim=range(-60,47),cex.main=3, cex.axis=1.5,   notch = TRUE, outline = FALSE, main="Coding Regions", lwd=3,border= c("black", "#225EA8","#FED976",  "#B10026" ), col = c("gray", "lightblue", "#FFFFB2","#EF3B2C"), names=c("Single-copy","Conserved","Specialized" , "Neofunctionalized"))

abline(h = 0, lty = 2, lwd=2)
title( ylab=expression("Signed    "^2) , cex.lab=1.5)

title( ylab="Signed  " , cex.lab=1.5)
text(2.5,43,expression(paste(italic("P < 0.001")," for all pairwise comparisons")), cex=1.5)

text(3,21, expression(italic("\u03B2")~plain(":")~italic("\u03B1")^2) )


# ########### Introns ############


# sink("hka.intron.txt")

current= setwd("/data/HKA/M_script/introns")

rm.con=read.table("dmel.rm.all.nest.within.r549.conchild.w10000s1.intron.hka")
rm.neo=read.table("dmel.rm.all.nest.within.r549.neochild.w10000s1.intron.hka")
rm.spe=read.table("dmel.rm.all.nest.within.r549.spechild.w10000s1.intron.hka")
rm.s=read.table("dmel.rm.all.nest.within.r549.singles.w10000s1.intron.hka")

median(rm.neo$V3)
median(rm.con$V3)
median(rm.spe$V3)
median(rm.s$V3)
# print("Intron")


# median(intron.neo$V2)
# median(intron.spec$V3)
# median(intron.con$V2)
# median(intron.s$V2)

# twoSamplePermutationTestLocation(intron.neo$V2, intron.spec$V3, fcn="median", alternative="two.sided", n.permutations=1000)

# twoSamplePermutationTestLocation(intron.neo$V2, intron.con$V2, fcn="median", alternative="two.sided", n.permutations=1000)

# twoSamplePermutationTestLocation(intron.neo$V2, intron.s$V2, fcn="median", alternative="two.sided", n.permutations=1000)

# twoSamplePermutationTestLocation(intron.spec$V3, intron.s$V2, fcn="median", alternative="two.sided", n.permutations=1000)

# twoSamplePermutationTestLocation(intron.spec$V3, intron.con$V2, fcn="median", alternative="two.sided", n.permutations=1000)

# sink()


SIGN.test(intron.neo$V2, md=0, alternative = "two.sided", conf.level = 0.95)
SIGN.test(intron.spec$V3, md=0, alternative = "two.sided", conf.level = 0.95)
SIGN.test(intron.con$V2, md=0, alternative = "two.sided", conf.level = 0.95)
SIGN.test(intron.s$V2, md=0, alternative = "two.sided", conf.level = 0.95)


binom.test(length(rm.neo[rm.neo$V3>0,]$V3),length(rm.neo$V3), p=0.5, alternative = "two.sided", conf.level = 0.95 )
binom.test(length(rm.con[rm.con$V3>0,]$V3),length(rm.con$V3), p=0.5, alternative = "two.sided", conf.level = 0.95 )
binom.test(length(rm.spec[rm.spec$V3>0,]$V3),length(rm.spec$V3), p=0.5, alternative = "two.sided", conf.level = 0.95 )
binom.test(length(rm.s[rm.s$V3>0,]$V3),length(rm.s$V3), p=0.5, alternative = "two.sided", conf.level = 0.95 )



boxplot(intron.neo$V2, intron.spec$V3, intron.con$V2, intron.s$V2,ylim=range(-300,260), notch = TRUE, outline = FALSE, main="HKA on Introns", names=c("neo","spec", "conchild", "single"), col=c("lavender", "lightcyan","mintcream","lightgray"))
boxplot(intron.neo$V2, intron.spec$V3, intron.con$V2, intron.s$V2, lwd=2, ylim=range(-300,260),cex.main=3, cex.axis=1.5,  notch = TRUE, outline = FALSE, main="Intron Region", names=c("neofunctionalized","specialized", "conserved", "single"), col=c("#CB181D", "#FB6A4A","#FCAE91", "#FEE5D9"))
boxplot(intron.neo$V2, intron.spec$V3, intron.con$V2, intron.s$V2, lwd=2, ylim=range(-300,260),cex.main=3, cex.axis=1.5,  notch = TRUE, outline = FALSE, main="Intron Region", names=c("neofunctionalized","specialized", "conserved", "single"), col=c("#FBB4AE", "#B3CDE3","#CCEBC5", "#DECBE4"))

boxplot(intron.neo$V2, intron.spec$V3, intron.con$V2, intron.s$V2,ylim=range(-300,260), lwd=2, ylab=expression(paste("Signed"," ", chi^2)) ,cex.main=3,cex.lab=0.95, cex.axis=1.5,  notch = TRUE, outline = FALSE, main="Introns", names=c("Neofunctionalized","Specialized", "Conserved", "Single-copy"), col=c("#FBB4AE", "#B3CDE3","#CCEBC5", "#DECBE4"))
text(3.7,250,"P<0.001 for all comparisons", cex=1.25)

boxplot(rm.s$V3, rm.neo$V3, rm.spe$V3, rm.con$V3,  boxwex=0.5, ylim=range(-90,80),lwd=5, cex.main=3, cex.axis=1.5,   notch = TRUE, outline = FALSE, main="Introns", names=c("Single-copy", "Neofunctionalized","Specialized", "Conserved" ),border =  c("gray57", "#EF3B2C","#FED976",  "#1D91C0"))

boxplot(rm.s$V3, rm.con$V3,  rm.spe$V3, rm.neo$V3, boxwex=0.5, ylim=range(-95,80), cex.main=3, cex.axis=1.5,   notch = TRUE, outline = FALSE, main="Introns", lwd=3,border= c("black", "#225EA8","#FED976",  "#B10026" ), col = c("gray", "lightblue", "#FFFFB2","#EF3B2C"), names=c("Single-copy","Conserved","Specialized" , "Neofunctionalized"))


abline(h = 0, lty = 2, lwd=2)
title( ylab="Signed  " , cex.lab=1.5)
title( ylab=expression("Signed    "^2) , cex.lab=1.5)

#title( ylab=expression(paste("Signed"," ", Chi^2)) , cex.lab=1.5)
text(2.5,72,expression(paste(italic("P < 0.001")," for all pairwise comparisons")), cex=1.5)


########### 5 UTR   ############


 # sink("hka.3utr.txt")


utr5.s <- read.table("dmel.r549.single.w10000s1.5utr.hka")

utr5.con <- read.table("hka.con.fivep.w10000s1.r549.child")
utr5.neo <- read.table("hka.neo.fivep.w10000s1.r549.child")
utr5.spec <- read.table("dmel.r549.spec.w10000s1.5utr.hka")

# print("utr5")


# median(utr5.neo$V3)
# median(utr5.spec$V2)
# median(utr5.con$V3)
# median(utr5.s$V2)

# twoSamplePermutationTestLocation(utr5.neo$V3, utr5.spec$V2, fcn="median", alternative="two.sided", n.permutations=1000)

# twoSamplePermutationTestLocation(utr5.neo$V3, utr5.con$V3, fcn="median", alternative="two.sided", n.permutations=1000)

# twoSamplePermutationTestLocation(utr5.neo$V3, utr5.s$V2, fcn="median", alternative="two.sided", n.permutations=1000)

# twoSamplePermutationTestLocation(utr5.spec$V2, utr5.s$V2, fcn="median", alternative="two.sided", n.permutations=1000)

# twoSamplePermutationTestLocation(utr5.spec$V2, utr5.con$V3, fcn="median", alternative="two.sided", n.permutations=1000)

boxplot(utr5.neo$V3, utr5.spec$V2, utr5.con$V3, utr5.s$V2,ylim=range(-80,85), notch = TRUE, outline = FALSE, main="HKA on 5'UTR", names=c("neo","spec", "conchild", "single"), col=c("lavender", "lightcyan","mintcream","lightgray"))
boxplot(utr5.neo$V3, utr5.spec$V2, utr5.con$V3, utr5.s$V2, lwd=2, ylim=range(-80,85),cex.main=3, cex.axis=1.5,  notch = TRUE, outline = FALSE, main="5' UTR Region", names=c("neofunctionalized","specialized", "conserved", "single"), col=c("#CB181D", "#FB6A4A","#FCAE91", "#FEE5D9"))
boxplot(utr5.neo$V3, utr5.spec$V2, utr5.con$V3, utr5.s$V2, lwd=2, ylim=range(-80,85),cex.main=3, cex.axis=1.5,  notch = TRUE, outline = FALSE, main="5' UTR Region", names=c("neofunctionalized","specialized", "conserved", "single"), col=c("#FBB4AE", "#B3CDE3","#CCEBC5", "#DECBE4"))

boxplot(utr5.neo$V3, utr5.spec$V2, utr5.con$V3, utr5.s$V2,ylim=range(-80,85), lwd=2, ylab=expression(paste("Signed"," ", chi^2)) ,cex.main=3,cex.lab=1.25, cex.axis=1.25,  notch = TRUE, outline = FALSE, main="5' UTRs", names=c("Neofunctionalized","Specialized", "Conserved", "Single-copy"), col=c("#FBB4AE", "#B3CDE3","#CCEBC5", "#DECBE4"))



SIGN.test(utr5.neo$V3, md=0, alternative = "two.sided", conf.level = 0.95)
SIGN.test(utr5.spec$V2, md=0, alternative = "two.sided", conf.level = 0.95)
SIGN.test(utr5.con$V3, md=0, alternative = "two.sided", conf.level = 0.95)
SIGN.test(utr5.s$V2, md=0, alternative = "two.sided", conf.level = 0.95)


boxplot(utr5.s$V2, utr5.neo$V3, utr5.spec$V2, utr5.con$V3,  boxwex=0.5, lwd=5, ylim=range(-70,80),cex.main=3, cex.axis=1.5,   notch = TRUE, outline = FALSE, main="5' UTRs", names=c("Single-copy", "Neofunctionalized","Specialized", "Conserved" ),border =  c("gray57", "#EF3B2C","#FED976",  "#1D91C0"))

boxplot(utr5.s$V2, utr5.con$V3, utr5.spec$V2, utr5.neo$V3,   boxwex=0.5, ylim=range(-80,85),cex.main=3, cex.axis=1.5,   notch = TRUE, outline = FALSE, main="5' UTRs", lwd=3,border= c("black", "#225EA8","#FED976",  "#B10026" ), col = c("gray", "lightblue", "#FFFFB2","#EF3B2C"), names=c("Single-copy","Conserved","Specialized" , "Neofunctionalized"))

abline(h = 0, lty = 2, lwd=2)
title( ylab="Signed  " , cex.lab=1.5)
title( ylab=expression("Signed    "^2) , cex.lab=1.5)

#title( ylab=expression(paste("Signed"," ", Chi^2)) , cex.lab=1.5)


#lines(c(1,2), c(70,70))
#text(1.5,29,"***", cex=2)
lines(c(1,4), c(82,82))
text(2.5,84,"*", cex=2)
lines(c(1,1), c(82,80.5))
lines(c(4,4), c(82,80.5))

lines(c(2,3), c(68,68))
lines(c(3,3), c(68,66.5))
lines(c(2,2), c(68,66.5))

text(2.5,70,"*", cex=2)
#lines(c(2,4), c(80,80))
#text(3,44,"***", cex=2)


lines(c(2,4), c(75,75))
lines(c(2,2), c(75,73.5))
lines(c(4,4), c(75,73.5))
text(3,77,"***", cex=2)
#lines(c(3,4), c(65,65))
#text(3.5,48,"***", cex=2)

########### 3 UTR   ############


utr3.s <- read.table("dmel.r549.single.w10000s1.3utr.hka")

utr3.con <- read.table("hka.con.threep.w10000s1.r549.child")
utr3.neo <- read.table("hka.neo.threep.w10000s1.r549.child")
utr3.spec <- read.table("dmel.r549.spec.w10000s1.3utr.hka")


# print("utr3")


# median(utr3.neo$V3)
# median(utr3.spec$V2)
# median(utr3.con$V3)
# median(utr3.s$V2)

# twoSamplePermutationTestLocation(utr3.neo$V3, utr3.spec$V2, fcn="median", alternative="two.sided", n.permutations=1000)

# twoSamplePermutationTestLocation(utr3.neo$V3, utr3.con$V3, fcn="median", alternative="two.sided", n.permutations=1000)

# twoSamplePermutationTestLocation(utr3.neo$V3, utr3.s$V2, fcn="median", alternative="two.sided", n.permutations=1000)

# twoSamplePermutationTestLocation(utr3.spec$V2, utr3.s$V2, fcn="median", alternative="two.sided", n.permutations=1000)

# twoSamplePermutationTestLocation(utr3.spec$V2, utr3.con$V3, fcn="median", alternative="two.sided", n.permutations=1000)

# sink()

SIGN.test(utr3.neo$V3, md=0, alternative = "two.sided", conf.level = 0.95)
SIGN.test(utr3.spec$V2, md=0, alternative = "two.sided", conf.level = 0.95)
SIGN.test(utr3.con$V3, md=0, alternative = "two.sided", conf.level = 0.95)
SIGN.test(utr3.s$V2, md=0, alternative = "two.sided", conf.level = 0.95)



boxplot(utr3.neo$V3, utr3.spec$V2, utr3.con$V3, utr3.s$V2, ylim=range(-60,50),notch = TRUE, outline = FALSE, main="HKA on 3'UTR", names=c("neo","spec", "conchild", "single"), col=c("lavender", "lightcyan","mintcream","lightgray"))
boxplot(utr3.neo$V3, utr3.spec$V2, utr3.con$V3, utr3.s$V2, ylim=range(-60,50), lwd=2, cex.main=3, cex.axis=1.5,  notch = TRUE, outline = FALSE, main="3' UTR Region", names=c("neofunctionalized","specialized", "conserved", "single"), col=c("#CB181D", "#FB6A4A","#FCAE91", "#FEE5D9"))
boxplot(utr3.neo$V3, utr3.spec$V2, utr3.con$V3, utr3.s$V2, ylim=range(-60,50),lwd=2, cex.main=3, cex.axis=1.5,  notch = TRUE, outline = FALSE, main="3' UTR Region", names=c("neofunctionalized","specialized", "conserved", "single"), col=c("#FBB4AE", "#B3CDE3","#CCEBC5", "#DECBE4"))

boxplot(utr3.neo$V3, utr3.spec$V2, utr3.con$V3, utr3.s$V2, ylim=range(-60,50), lwd=2, ylab=expression(paste("Signed"," ", chi^2)) ,cex.main=3,cex.lab=1.25, cex.axis=1.25,  notch = TRUE, outline = FALSE, main="3' UTRs", names=c("Neofunctionalized","Specialized", "Conserved", "Single-copy"), col=c("#FBB4AE", "#B3CDE3","#CCEBC5", "#DECBE4"))
text(3.7,45,"P<0.001 for all comparisons", cex=1.25)

boxplot(utr3.s$V2, utr3.neo$V3, utr3.spec$V2, utr3.con$V3,  boxwex=0.5, lwd=5, ylim=range(-60,50),cex.main=3, cex.axis=1.5,   notch = TRUE, outline = FALSE, main="3' UTRs", names=c("Single-copy", "Neofunctionalized","Specialized", "Conserved" ),border =  c("gray57", "#EF3B2C","#FED976",  "#1D91C0"))

boxplot(utr3.s$V2, utr3.con$V3, utr3.spec$V2,  utr3.neo$V3,  boxwex=0.5,  ylim=range(-85,50),cex.main=3, cex.axis=1.5,   notch = TRUE, outline = FALSE, main="3' UTRs",lwd=3,border= c("black", "#225EA8","#FED976",  "#B10026" ), col = c("gray", "lightblue", "#FFFFB2","#EF3B2C"), names=c("Single-copy","Conserved","Specialized" , "Neofunctionalized"))


abline(h = 0, lty = 2, lwd=2)
title( ylab="Signed  " , cex.lab=1.5)
title( ylab=expression("Signed    "^2) , cex.lab=1.5)
#title( ylab=expression(paste("Signed"," ", Chi^2)) , cex.lab=1.5)
text(2.5,43,expression(paste(italic("P < 0.001")," for all pairwise comparisons")), cex=1.5)


############################## whole gene region ########################
# sink("hka.wholegene.txt")

whole.neo <- c(utr3.neo$V3, utr5.neo$V3, intron.neo$V2, cds.neo$V3)
whole.spec <- c(utr3.spec$V2, utr5.spec$V2, intron.spec$V2, cds.spec$V3)
whole.con <- c(utr3.con$V3, utr5.con$V3, intron.con$V2, cds.con$V3)
whole.s <- c(utr3.s$V2, utr5.s$V2, intron.s$V2, cds.s$V3)

# print("Whole")

# median(whole.neo)
# median(whole.spec)
# median(whole.con)
# median(whole.s)

# twoSamplePermutationTestLocation(whole.neo, whole.spec, fcn="median", alternative="two.sided", n.permutations=1000)

# twoSamplePermutationTestLocation(whole.neo, whole.con, fcn="median", alternative="two.sided", n.permutations=1000)

# twoSamplePermutationTestLocation(whole.neo, whole.s, fcn="median", alternative="two.sided", n.permutations=1000)

# twoSamplePermutationTestLocation(whole.spec, whole.s, fcn="median", alternative="two.sided", n.permutations=1000)

# twoSamplePermutationTestLocation(whole.spec, whole.con, fcn="median", alternative="two.sided", n.permutations=1000)

# twoSamplePermutationTestLocation(whole.s, whole.con, fcn="median", alternative="two.sided", n.permutations=1000)


# sink()

boxplot(whole.neo, whole.spec, whole.con, whole.s, notch = TRUE, outline = FALSE,ylim=range(-230,200),  main="HKA on whole gene region", names=c("neo","spec", "conchild", "single"), col=c("lavender", "lightcyan","mintcream","lightgray"))


lines(c(1,2), c(150,150))
text(1.5,155,"***", cex=2)
lines(c(1,3), c(170,170))
text(2,175,"***", cex=2)
lines(c(2,3), c(120,120))
text(2.5,125,"***", cex=2)
lines(c(2,4), c(140,140))
text(3,145,"***", cex=2)
lines(c(1,4), c(185,185))
text(2.5,190,"***", cex=2)
lines(c(3,4), c(130,130))
text(3.5,5,"***", cex=2)
