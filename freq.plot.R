current=setwd("/data/HKA/M_script/FreqSpectrum/results")



### CDS CPS ###
cds.child=read.table("hka.r549.child.cds.freqspec")
#hist(cds.child$V2, main="SNP frequency spectrum of child copies in CDS", breaks = 100, xlab="Allele frequency", col="red")
cds.child$group=rep("child", nrow(cds.child))
cds.child[["group"]]=rep("child", nrow(cds.child))
cds.child[,"group"]=rep("child", nrow(cds.child))
cds.child$group="child"


cds.parent=read.table("hka.r549.parent.cds.freqspec")
#hist(cds.parent$V2, main="SNP frequency spectrum of parent copies in CDS", breaks = 100, xlab="Allele frequency", col="red")
cds.parent$group=rep("child", nrow(cds.parent))
cds.parent[["group"]]=rep("child", nrow(cds.parent))
cds.parent[,"group"]=rep("child", nrow(cds.parent))
cds.parent$group="parent"
View(cds.parent)

cds.single=read.table("hka.r549.single.cds.freqspec")
cds.single$group=rep("child", nrow(cds.single))
cds.single[["group"]]=rep("child", nrow(cds.single))
cds.single[,"group"]=rep("child", nrow(cds.single))
cds.single$group="single"
View(cds.single)



cds.cps=rbind(cds.child, cds.parent, cds.single)
cds.cps=na.omit(cds.cps)

ggplot(cds.cps, aes(x=V2, fill=group))+
  geom_histogram(aes(y=c(..count..[..group..==1]/sum(..count..[..group..==1]),
                         ..count..[..group..==2]/sum(..count..[..group..==2]),
                         ..count..[..group..==3]/sum(..count..[..group..==3]))),
                 bins = 50, position = "dodge",alpha=0.8 )+
  ggtitle("Frequency spectrum for child, parent and single copy genes of CDS region")+
  xlab("Allele frequency")+
  ylab("Percent")+
  scale_fill_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A") )



### CDS neocon ###
cds.child.neo=read.table("hka.r549.neochild.cds.freqspec")
#hist(cds.child.neo$V2, main="SNP frequency spectrum of neochild copies in CDS", breaks = 100, xlab="Allele frequency", col="red")
cds.child.neo$group=rep("child", nrow(cds.child.neo))
cds.child.neo[["group"]]=rep("child", nrow(cds.child.neo))
cds.child.neo[,"group"]=rep("child", nrow(cds.child.neo))
cds.child.neo$group="child.neo"


cds.child.con=read.table("hka.r549.conchild.cds.freqspec")
#hist(cds.child.con$V2, main="SNP frequency spectrum of conchild copies in CDS", breaks = 100, xlab="Allele frequency", col="red")
cds.child.con$group=rep("child", nrow(cds.child.con))
cds.child.con[["group"]]=rep("child", nrow(cds.child.con))
cds.child.con[,"group"]=rep("child", nrow(cds.child.con))
cds.child.con$group="child.con"

cds.neocon=rbind(cds.child.neo, cds.child.con, cds.single)


ggplot(cds.neocon, aes(x=V2, fill=group))+
  geom_histogram(aes(y=c(..count..[..group..==1]/sum(..count..[..group..==1]),
                         ..count..[..group..==2]/sum(..count..[..group..==2]),
                         ..count..[..group..==3]/sum(..count..[..group..==3]))),
                 bins = 50, position = "dodge",alpha=0.8 )+
  ggtitle("Frequency spectrum for neofunctionalized, conserved child copy and single copy genes of CDS region")+
  xlab("Allele frequency")+
  ylab("Percent")+
  scale_fill_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A") )



### 5utr CPS ###

d.5utr.child=read.table("hka.r549.child.5utr.freqspec")
#hist(d.5utr.child$V2, main="SNP frequency spectrum of child copies in d.5utr", breaks = 100, xlab="Allele frequency", col="red")
d.5utr.child$group=rep("child", nrow(d.5utr.child))
d.5utr.child[["group"]]=rep("child", nrow(d.5utr.child))
d.5utr.child[,"group"]=rep("child", nrow(d.5utr.child))
d.5utr.child$group="child"


d.5utr.parent=read.table("hka.r549.parent.5utr.freqspec")
#hist(d.5utr.parent$V2, main="SNP frequency spectrum of parent copies in d.5utr", breaks = 100, xlab="Allele frequency", col="red")
d.5utr.parent$group=rep("child", nrow(d.5utr.parent))
d.5utr.parent[["group"]]=rep("child", nrow(d.5utr.parent))
d.5utr.parent[,"group"]=rep("child", nrow(d.5utr.parent))
d.5utr.parent$group="parent"

d.5utr.single=read.table("hka.r549.single.5utr.freqspec")
d.5utr.single$group=rep("child", nrow(d.5utr.single))
d.5utr.single[["group"]]=rep("child", nrow(d.5utr.single))
d.5utr.single[,"group"]=rep("child", nrow(d.5utr.single))
d.5utr.single$group="single"


ggplot(cds.cps, aes(x=V2, fill=group))+
  geom_histogram(aes(y=c(..count..[..group..==1]/sum(..count..[..group..==1]),
                         ..count..[..group..==2]/sum(..count..[..group..==2]),
                         ..count..[..group..==3]/sum(..count..[..group..==3]))),
                 bins = 50, position = "dodge",alpha=0.8 )+
  ggtitle("Frequency spectrum for child, parent and single copy genes of 5' UTR region")+
  xlab("Allele frequency")+
  ylab("Percent")+
  scale_fill_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A") )





### 5utr neocon ###

d.5utr.child.neo=read.table("hka.r549.neochild.5utr.freqspec")
#hist(d.5utr.child.neo$V2, main="SNP frequency spectrum of neochild copies in d.5utr", breaks = 100, xlab="Allele frequency", col="red")
d.5utr.child.neo$group=rep("child", nrow(d.5utr.child.neo))
d.5utr.child.neo[["group"]]=rep("child", nrow(d.5utr.child.neo))
d.5utr.child.neo[,"group"]=rep("child", nrow(d.5utr.child.neo))
d.5utr.child.neo$group="child.neo"


d.5utr.child.con=read.table("hka.r549.conchild.5utr.freqspec")
#hist(d.5utr.child.con$V2, main="SNP frequency spectrum of conchild copies in d.5utr", breaks = 100, xlab="Allele frequency", col="red")
d.5utr.child.con$group=rep("child", nrow(d.5utr.child.con))
d.5utr.child.con[["group"]]=rep("child", nrow(d.5utr.child.con))
d.5utr.child.con[,"group"]=rep("child", nrow(d.5utr.child.con))
d.5utr.child.con$group="child.con"

d.5utr.neocon=rbind(d.5utr.child.neo, d.5utr.child.con, d.5utr.single)


ggplot(d.5utr.neocon, aes(x=V2, fill=group))+
  geom_histogram(aes(y=c(..count..[..group..==1]/sum(..count..[..group..==1]),
                         ..count..[..group..==2]/sum(..count..[..group..==2]),
                         ..count..[..group..==3]/sum(..count..[..group..==3]))),
                 bins = 50, position = "dodge",alpha=0.8 )+
  ggtitle("Frequency spectrum for neofunctionalized, conserved child copy and single copy genes of 5' UTR region")+
  xlab("Allele frequency")+
  ylab("Percent")+
  scale_fill_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A") )




### 3utr CPS ###

d.3utr.child=read.table("hka.r549.child.3utr.freqspec")
#hist(d.3utr.child$V2, main="SNP frequency spectrum of child copies in d.3utr", breaks = 100, xlab="Allele frequency", col="red")
d.3utr.child$group=rep("child", nrow(d.3utr.child))
d.3utr.child[["group"]]=rep("child", nrow(d.3utr.child))
d.3utr.child[,"group"]=rep("child", nrow(d.3utr.child))
d.3utr.child$group="child"


d.3utr.parent=read.table("hka.r549.parent.3utr.freqspec")
#hist(d.3utr.parent$V2, main="SNP frequency spectrum of parent copies in d.3utr", breaks = 100, xlab="Allele frequency", col="red")
d.3utr.parent$group=rep("child", nrow(d.3utr.parent))
d.3utr.parent[["group"]]=rep("child", nrow(d.3utr.parent))
d.3utr.parent[,"group"]=rep("child", nrow(d.3utr.parent))
d.3utr.parent$group="parent"


d.3utr.single=read.table("hka.r549.single.3utr.freqspec")
d.3utr.single$group=rep("child", nrow(d.3utr.single))
d.3utr.single[["group"]]=rep("child", nrow(d.3utr.single))
d.3utr.single[,"group"]=rep("child", nrow(d.3utr.single))
d.3utr.single$group="single"


ggplot(cds.cps, aes(x=V2, fill=group))+
  geom_histogram(aes(y=c(..count..[..group..==1]/sum(..count..[..group..==1]),
                         ..count..[..group..==2]/sum(..count..[..group..==2]),
                         ..count..[..group..==3]/sum(..count..[..group..==3]))),
                 bins = 50, position = "dodge",alpha=0.8 )+
  ggtitle("Frequency spectrum for child, parent and single copy genes of 3' UTR region")+
  xlab("Allele frequency")+
  ylab("Percent")+
  scale_fill_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A") )





### 3utr neocon ###

d.3utr.child.neo=read.table("hka.r549.neochild.3utr.freqspec")
#hist(d.3utr.child.neo$V2, main="SNP frequency spectrum of neochild copies in d.3utr", breaks = 100, xlab="Allele frequency", col="red")
d.3utr.child.neo$group=rep("child", nrow(d.3utr.child.neo))
d.3utr.child.neo[["group"]]=rep("child", nrow(d.3utr.child.neo))
d.3utr.child.neo[,"group"]=rep("child", nrow(d.3utr.child.neo))
d.3utr.child.neo$group="child.neo"


d.3utr.child.con=read.table("hka.r549.conchild.3utr.freqspec")
#hist(d.3utr.child.con$V2, main="SNP frequency spectrum of conchild copies in d.3utr", breaks = 100, xlab="Allele frequency", col="red")
d.3utr.child.con$group=rep("child", nrow(d.3utr.child.con))
d.3utr.child.con[["group"]]=rep("child", nrow(d.3utr.child.con))
d.3utr.child.con[,"group"]=rep("child", nrow(d.3utr.child.con))
d.3utr.child.con$group="child.con"

d.3utr.neocon=rbind(d.3utr.child.neo, d.3utr.child.con, d.3utr.single)


ggplot(d.3utr.neocon, aes(x=V2, fill=group))+
  geom_histogram(aes(y=c(..count..[..group..==1]/sum(..count..[..group..==1]),
                         ..count..[..group..==2]/sum(..count..[..group..==2]),
                         ..count..[..group..==3]/sum(..count..[..group..==3]))),
                 bins = 50, position = "dodge",alpha=0.8 )+
  ggtitle("Frequency spectrum for neofunctionalized, conserved child copy and single copy genes of 3' UTR region")+
  xlab("Allele frequency")+
  ylab("Percent")+
  scale_fill_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A") )










