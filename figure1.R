current1=setwd("/data/Progress_Report/2017-3-15")

# Ka/Ks
k.neo=read.table("kaks.neo.2.child.txt")
k.spec=read.table("child.spec.kaks.txt")
k.con=read.table("kaks.con.child.txt")
#k.p=read.table("parent.sim.paml.neocon.txt")
k.s=read.table("single_sim_paml.txt")

#k.pneo=read.table("parent.sim.paml.neo.txt")
#k.pcon=read.table("parent.sim.paml.con.txt")

median(k.neo$V2)
#[1] 0.3299
median(k.con$V3)
#[1] 0.2314
#median(k.p$V5)
median(k.s$V5)
#[1] 0.0795

#median(k.pneo$V5)
#median(k.pcon$V5)
twoSamplePermutationTestLocation(k.neo$V2, k.spec$V2, fcn="median", alternative="two.sided", n.permutations=1000)

twoSamplePermutationTestLocation(k.neo$V2, k.con$V3, fcn="median", alternative="two.sided", n.permutations=1000)

twoSamplePermutationTestLocation(k.neo$V2, k.s$V5, fcn="median", alternative="two.sided", n.permutations=1000)

twoSamplePermutationTestLocation(k.spec$V2, k.s$V5, fcn="median", alternative="two.sided", n.permutations=1000)

twoSamplePermutationTestLocation(k.spec$V2, k.con$V3, fcn="median", alternative="two.sided", n.permutations=1000)

twoSamplePermutationTestLocation(k.s$V5, utr3.con$V3, fcn="median", alternative="two.sided", n.permutations=1000)

wilcox.test(k.neoc$V2, k.con$V3)
#wilcox.test(k.neoc$V3, k.p$V5)
#wilcox.test(k.con$V3, k.p$V5)
wilcox.test(k.neoc$V2, k.s$V5)
wilcox.test(k.con$V3, k.s$V5)
#wilcox.test(k.p$V5, k.s$V5)

#wilcox.test(k.pneo$V5, k.pcon$V5)

boxplot(k.neo$V2,k.spec$V2, k.con$V3, k.s$V5,ylim=range(0,1.2), lwd=2, cex.main=3, cex.axis=1.5,notch = TRUE, outline = FALSE, main="Ka/Ks ratio", names=c("neo","spec", "conchild", "single"), col=c("lavender", "lightcyan","mintcream", "light gray"))
boxplot(k.neo$V2,k.spec$V2, k.con$V3, k.s$V5,ylim=range(0,1.2), lwd=2, cex.main=3, cex.axis=1.5,notch = TRUE, outline = FALSE, main="Ka/Ks Ratio", names=c("neofunctionalized","specialized", "conserved", "single"), col=c("#CB181D", "#FB6A4A","#FCAE91", "#FEE5D9"))
boxplot(k.neo$V2,k.spec$V2, k.con$V3, k.s$V5,ylim=range(0,1.2), lwd=2, cex.main=3, cex.axis=1.5,notch = TRUE, outline = FALSE, main="Ka/Ks Ratio", names=c("neofunctionalized","specialized", "conserved", "single"), col=c("#FBB4AE", "#B3CDE3","#CCEBC5", "#DECBE4"))

boxplot(k.neo$V2,k.spec$V2, k.con$V3, k.s$V5,ylim=range(0,1.1), ylab="Ka/Ks", lwd=2, cex.lab=1.5, cex.axis=1.5,notch = TRUE, outline = FALSE,  names=c("Neofunctionalized","Specialized", "Conserved", "Single-copy"), col=c("#FBB4AE", "#B3CDE3","#CCEBC5", "#DECBE4"))

boxplot(k.s$V5,k.neo$V2,k.spec$V2, k.con$V3, ylim=range(0,1.1), ylab=expression(paste(italic("K"['a']), italic(" / "), italic("K"['s']))), boxwex=0.8, lwd=4, cex.lab=1.2, cex.axis=1.2,notch = TRUE, outline = FALSE, border= c("lightgray", "#F03B20", "#FFEDA0", "lightblue"),  names=c("Single-copy", "Neofunctionalized","Specialized", "Conserved" ))

#V1
boxplot(res=5000, k.s$V5, k.con$V3, k.spec$V2,k.neo$V2, ylim=range(0,1.2), ylab=expression(paste(italic("K"['a']), italic(" / "), italic("K"['s']))), boxwex=0.5, lwd=2, cex.lab=1.2, cex.axis=1.2,notch = TRUE, outline = FALSE, border= c("black", "#225EA8","#FED976",  "#B10026" ), col = c("gray", "lightblue", "#FFFFB2","#EF3B2C"), names=c("Single-copy","Conserved","Specialized" , "Neofunctionalized"))

lines(c(1,2), c(1.05,1.05))
lines(c(1,1), c(1.04,1.05) )
lines(c(2,2), c(1.04,1.05) )
text(1.5,1.07,"***", cex=2)
#lines(c(1,3), c(1.04,1.04))

lines(c(1,3), c(1.1,1.1))
lines(c(1,1), c(1.1,1.09) )
lines(c(3,3), c(1.1,1.09) )
text(2,1.12,"**", cex=2)

lines(c(1,4), c(1.15,1.15))
lines(c(1,1), c(1.14,1.15) )
lines(c(4,4), c(1.14,1.15) )
text(2.5,1.17,"***", cex=2)
#lines(c(2,3), c(1.09,1.09))
#text(2.5,1,"***", cex=2)

lines(c(3,4), c(1,1))
lines(c(3,3), c(0.99,1) )
lines(c(4,4), c(0.99,1) )
text(3.5,1.02,"*", cex=2)

#V2
boxplot(k.s$V5,k.neo$V2,k.spec$V2, k.con$V3, ylim=range(0,1.1), ylab=expression(paste(italic("K"['a']), italic(" / "), italic("K"['s']))), boxwex=0.5, lwd=3, cex.lab=1.2, cex.axis=1.2,notch = TRUE, outline = FALSE, border= c("gray57", "#EF3B2C","#FED976",  "#1D91C0"), names=c("Single-copy", "Neofunctionalized","Specialized", "Conserved" ))
#V3
boxplot(k.s$V5,k.neo$V2,k.spec$V2, k.con$V3, ylim=range(0,1.1), ylab=expression(paste(italic("K"['a']), italic(" / "), italic("K"['s']))), boxwex=0.8, lwd=1, cex.lab=1.2, cex.axis=1.2,notch = TRUE, outline = FALSE, col= c("lightgray", "#EF3B2C","#FFEDA0",  "lightblue"), names=c("Single-copy", "Neofunctionalized","Specialized", "Conserved" ))


display.brewer.all(type = "all")

brewer.pal(7,"YlOrRd")
brewer.pal(7, "YlGnBu")
brewer.pal(7, "Reds")
display.brewer.pal(9,"YlOrRd")
display.brewer.pal(7,"YlGnBu")
display.brewer.pal(7,"Reds")


lines(c(1,2), c(0.8,0.8))
lines(c(1,1), c(0.79,0.8) )
lines(c(2,2), c(0.79,0.8) )
text(1.5,0.82,"***", cex=2)
#lines(c(1,3), c(1.04,1.04))

lines(c(1,4), c(1.06,1.06))
lines(c(1,1), c(1.06,1.05) )
lines(c(4,4), c(1.06,1.05) )
text(2.5,1.08,"***", cex=2)
#lines(c(2,3), c(1.09,1.09))
#text(2.5,1,"***", cex=2)

lines(c(1,3), c(0.97,0.97))
lines(c(1,1), c(0.97,0.96) )
lines(c(3,3), c(0.97,0.96) )
text(2,0.99,"**", cex=2)


lines(c(2,3), c(0.89,0.89))
lines(c(3,3), c(0.88,0.89) )
lines(c(2,2), c(0.88,0.89) )
text(2.5,0.91,"*", cex=2)

# Neutrality index
mk.neoc=read.table("child.neo.mk.txt")
mk.spec=read.table("child.spec.mk.txt")
mk.con=read.table("conc.mk.alpha.raw.txt")
#mk.p=read.table("parent.sim.mk.neocon.txt")
mk.s=read.table("single.mk.alpha.raw.txt")

#mk.pneo=read.table("parent.sim.mk.neo.txt")
#mk.pcon=read.table("parent.sim.mk.con.txt")

neutral <- function(x){
  return(-(x-1))
}

ni.neoc=neutral(mk.neoc$V2)
ni.spec=neutral(mk.spec$V2)
ni.con=neutral(mk.con$V6)
#ni.p=neutral(mk.p$V6)
ni.s=neutral(mk.s$V6)

#ni.pneo=neutral(mk.pneo$V6)
#ni.pcon=neutral(mk.pcon$V6)

median(ni.neoc)
#[1] 0.9489796
median(ni.spec)
#[1] 1.219697
median(ni.con)
#[1] 0.6607143
#median(ni.p)
median(ni.s)
#[1] 1.615385

#median(ni.pneo)
#median(ni.pcon)

#wilcox.test(ni.neoc, ni.con)
#0.1737588
#wilcox.test(ni.neoc, ni.p)
#wilcox.test(ni.con, ni.p)
#wilcox.test(ni.neoc, ni.s)
#7.446311e-06

#wilcox.test(ni.con, ni.s)
#0.000425603
#wilcox.test(ni.p, ni.s)

#wilcox.test(ni.pneo, ni.pcon)

twoSamplePermutationTestLocation(ni.neoc, ni.spec, fcn="median", alternative="two.sided", n.permutations=1000)
# 0.251
twoSamplePermutationTestLocation(ni.neoc, ni.con, fcn="median", alternative="two.sided", n.permutations=1000)
# 0.321
twoSamplePermutationTestLocation(ni.neoc, ni.s, fcn="median", alternative="two.sided", n.permutations=1000)
#0.015
twoSamplePermutationTestLocation(ni.spec, ni.s, fcn="median", alternative="two.sided", n.permutations=1000)
# 0.371
twoSamplePermutationTestLocation(ni.spec, ni.con, fcn="median", alternative="two.sided", n.permutations=1000)
#0.07
twoSamplePermutationTestLocation(ni.s, ni.con, fcn="median", alternative="two.sided", n.permutations=1000)
#0.04


boxplot(ni.neoc, ni.spec, ni.con, ni.s, ylim=range(0,7),notch = TRUE, outline = FALSE, main="Neutrality index", names=c("neo", "spec", "conchild","single"), col=c("lavender", "lightcyan","mintcream", "lightgray"))
lines(c(1,2), c(5,5))
lines(c(1,3), c(5.5,5.5))
lines(c(1,4), c(6.7,6.7))
text(2.5,6.8,"*", cex=2)
lines(c(2,3), c(5.25,5.25))
lines(c(3,4), c(6,6))
text(3.5,6.2,"*", cex=2)
lines(c(2,4), c(6.5, 6.5))

# HKA
h.neoc=read.table("hka.specneo.child.txt")
h.con=read.table("hka.con.cds.w10000s1.r549.child")
#h.p=read.table("parent.sim.hka.neocon.txt")
h.s=read.table("single.w10000s1.r549.cds.hka")

#h.pneo=read.table("parent.sim.hka.neo.txt")
#h.pcon=read.table("parent.sim.hka.con.txt")

median(h.neoc$V2)
#[1] 0.527613
median(h.con$V3)
#[1] -13.55966
#median(h.p$V3)
median(h.s$V3)
#[1] -1.261694
#median(h.pneo$V3)
#median(h.pcon$V3)

wilcox.test(h.neoc$V2, h.con$V3)
#0
wilcox.test(h.neoc$V2, h.s$V3)
#0
wilcox.test(h.con$V3, h.s$V3)
#0
#wilcox.test(h.p$V3, h.s$V3)
#wilcox.test(h.neoc$V3, h.p$V3)
#wilcox.test(h.con$V3, h.p$V3)

#wilcox.test(h.pneo$V3, h.pcon$V3)

boxplot(h.neoc$V2, h.con$V3,h.s$V3,ylim=range(-60,45), notch = TRUE, outline = FALSE, main="HKA chi-square", names=c("neospec", "conchild", "single"), col=c("lavender", "lightcyan","mintcream"))
lines(c(1,2), c(26,26))
text(1.5,29,"***", cex=2)
lines(c(2,3), c(42,42))
text(2.5,44,"***", cex=2)
lines(c(1,3), c(36,36))
text(2,39,"***", cex=2)
