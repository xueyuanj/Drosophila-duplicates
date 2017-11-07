################# Old strigent Rule ########################
######################  D.ana  ##############################

#1) plot k values on the y-axis
#2) plot four mechanisms separately
#3) plot only newly arosed genes at each time point


# the K on y-axis
relax.ana.f=read.table("~/Documents/Progress_Report/2017-1-18/ana.all.strigent.txt")
head(relax.ana.f)

relax.ana.f$V5 <- as.character(relax.ana.f$V5)
relax.ana.f$V2 <- as.character(relax.ana.f$V2)
relax.ana.f$V1 <- as.character(relax.ana.f$V1)

relax.ana.f$V5[relax.ana.f$V5=="ana.n_a.strigent.txt"] <- 1
relax.ana.f$V5[relax.ana.f$V5=="ana.n_ms.strigent.txt"] <- 2
relax.ana.f$V5[relax.ana.f$V5=="ana.n_mel.strigent.txt"] <- 3
relax.ana.f$V5[relax.ana.f$V5=="ana.n_a.strigent.miss.txt"] <- 1
relax.ana.f$V5[relax.ana.f$V5=="ana.n_ms.strigent.miss.txt"] <- 2
relax.ana.f$V5[relax.ana.f$V5=="ana.n_mel.strigent.miss.txt"] <- 3

ggplot(data = relax.ana.f, aes(x=V5, y=log(V3),group=V1))+
  geom_line(size=0.3, colour="lightblue")+
  geom_point(shape=1,size=3, colour="lightblue")+
  geom_point( data = relax.ana.f[relax.ana.f$V4<0.05,], size=2, colour="red" )+
  geom_hline(aes(yintercept=log(1)), colour="black", linetype="dashed")+
  labs(title="Intesified selection indicated by RELAX (D.ana)", x="Branch label", y="log(K)")

ggplot(data = relax.ana.f[relax.ana.f$V2=="neochild",], aes(x=V5, y=log(V3),group=V1))+
  geom_line(size=0.3, colour="lightblue")+
  geom_point(shape=1,size=3, colour="lightblue")+
  geom_point( data = relax.ana.f[relax.ana.f$V4<0.05&relax.ana.f$V2=="neochild",], size=2, colour="red" )+
  geom_hline(aes(yintercept=log(1)), colour="black", linetype="dashed")+
  labs(title="Intesified selection indicated by RELAX (D.ana, neochild)", x="Branch label", y="log(K)")

ggplot(data = relax.ana.f[relax.ana.f$V2=="cons",], aes(x=V5, y=log(V3),group=V1))+
  geom_line(size=0.3, colour="lightblue")+
  geom_point(shape=1,size=3, colour="lightblue")+
  geom_point( data = relax.ana.f[relax.ana.f$V4<0.05&relax.ana.f$V2=="cons",], size=2, colour="red" )+
  geom_hline(aes(yintercept=log(1)), colour="black", linetype="dashed")+
  labs(title="Intesified selection indicated by RELAX (D.ana, conchild)", x="Branch label", y="log(K)")

ggplot(data = relax.ana.f[relax.ana.f$V2=="spec",], aes(x=V5, y=log(V3),group=V1))+
  geom_line(size=0.3, colour="lightblue")+
  geom_point(shape=1,size=3, colour="lightblue")+
  geom_point( data = relax.ana.f[relax.ana.f$V4<0.05&relax.ana.f$V2=="spec",], size=2, colour="red" )+
  geom_hline(aes(yintercept=log(1)), colour="black", linetype="dashed")+
  labs(title="Intesified selection indicated by RELAX (D.ana, specialized)", x="Branch label", y="log(K)")

ggplot(data = relax.ana.f[relax.ana.f$V2=="neoparent",], aes(x=V5, y=log(V3),group=V1))+
  geom_line(size=0.3, colour="lightblue")+
  geom_point(shape=1,size=3, colour="lightblue")+
  geom_point( data = relax.ana.f[relax.ana.f$V4<0.05&relax.ana.f$V2=="neoparent",], size=2, colour="red" )+
  geom_hline(aes(yintercept=log(1)), colour="black", linetype="dashed")+
  labs(title="Intesified selection indicated by RELAX (D.ana, neoparent)", x="Branch label", y="log(K)")





##########################################################################################################################


relax.pse.f=read.table("/data/allspecies/miss.sister/pk_values/pse.all.strigent.txt")
head(relax.pse.f)

relax.pse.f$V5 <- as.character(relax.pse.f$V5)
relax.pse.f$V2 <- as.character(relax.pse.f$V2)
relax.pse.f$V1 <- as.character(relax.pse.f$V1)

relax.pse.f$V5[relax.pse.f$V5=="pse.n_pse.k_pvalue.txt"] <- 1
relax.pse.f$V5[relax.pse.f$V5=="pse.n_ana.k_pvalue.txt"] <- 2
relax.pse.f$V5[relax.pse.f$V5=="pse.n_ms.k_pvalue.txt"] <- 3
relax.pse.f$V5[relax.pse.f$V5=="pse.n_mel.k_pvalue.txt"] <- 4
relax.pse.f$V5[relax.pse.f$V5=="pse.n_pse.strigent.txt"] <- 1
relax.pse.f$V5[relax.pse.f$V5=="pse.n_ana.strigent.txt"] <- 2
relax.pse.f$V5[relax.pse.f$V5=="pse.n_ms.strigent.txt"] <- 3
relax.pse.f$V5[relax.pse.f$V5=="pse.n_mel.strigent.txt"] <- 4

target<- relax.pse.f$V1[relax.pse.f$V4<0.05]
target.u <- unique(target)

sub.pse <- relax.pse.f[relax.pse.f$V1%in%target.u,]

#################################################################
##########     FINAL VERSION   ##################################
#################################################################
View(relax.pse.f)
ggplot(data = relax.pse.f, aes(x=V5, y=log10(V3),group=V1))+
  geom_line(size=0.5, aes(group=V1), col="lightblue")+
  #geom_line(size=0.3, colour=ifelse(relax.pse.f$V4<0.05, "red", "lightblue") )+
  #geom_line(data = sub.pse, size=0.2, colour="red" )+
  geom_point(shape=1,size=3, colour="lightblue")+
  geom_point( data = relax.pse.f[relax.pse.f$V4<0.05,], size=2, colour="red" )+
  geom_hline(aes(yintercept=log10(1)), colour="black", linetype="dashed")+
  labs( x="Test Branch", y="")+#expression(paste(italic(log)[10],"(", italic("K"), ")"  )) )+
  theme_classic()+
  theme(text = element_text(size = 20),panel.border = element_rect(colour = "black", fill=NA, size=0.5))





ggplot(data = relax.pse.f[relax.pse.f$V2=="neochild",], aes(x=V5, y=log(V3),group=V1))+
  geom_line(size=0.3, colour="lightblue")+
  geom_point(shape=1,size=3, colour="lightblue")+
  geom_point( data = relax.pse.f[relax.pse.f$V4<0.05&relax.pse.f$V2=="neochild",], size=2, colour="indianred" )+
  geom_hline(aes(yintercept=log(1)), colour="black", linetype="dashed")+
  labs(title="Intesified selection indicated by RELAX (D.pse, neochild)", x="Branch label", y="log(K)")

ggplot(data = relax.pse.f[relax.pse.f$V2=="cons",], aes(x=V5, y=log(V3),group=V1))+
  geom_line(size=0.3, colour="lightblue")+
  geom_point(shape=1,size=3, colour="lightblue")+
  geom_point( data = relax.pse.f[relax.pse.f$V4<0.05&relax.pse.f$V2=="cons",], size=2, colour="red" )+
  geom_hline(aes(yintercept=log(1)), colour="black", linetype="dashed")+
  labs(title="Intesified selection indicated by RELAX (D.pse, conchild)", x="Branch label", y="log(K)")

ggplot(data = relax.pse.f[relax.pse.f$V2=="spec",], aes(x=V5, y=log(V3),group=V1))+
  geom_line(size=0.3, colour="lightblue")+
  geom_point(shape=1,size=3, colour="lightblue")+
  geom_point( data = relax.pse.f[relax.pse.f$V4<0.05&relax.pse.f$V2=="spec",], size=2, colour="red" )+
  geom_hline(aes(yintercept=log(1)), colour="black", linetype="dashed")+
  labs(title="Intesified selection indicated by RELAX (D.pse, specialized)", x="Branch label", y="log(K)")

ggplot(data = relax.pse.f[relax.pse.f$V2=="neoparent",], aes(x=V5, y=log(V3),group=V1))+
  geom_line(size=0.3, colour="lightblue")+
  geom_point(shape=1,size=3, colour="lightblue")+
  geom_point( data = relax.pse.f[relax.pse.f$V4<0.05&relax.pse.f$V2=="neoparent",], size=2, colour="red" )+
  geom_hline(aes(yintercept=log(1)), colour="black", linetype="dashed")+
  labs(title="Intesified selection indicated by RELAX (D.pse, neoparent)", x="Branch label", y="log(K)")

