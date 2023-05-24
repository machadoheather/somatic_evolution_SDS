## April 2022 HEM

# Plotting metadata

## setup
library(ggplot2)
library(reshape2)
library(plyr)

## Plotting SDS metadata of patients for figures 1b and 1c

## checking colonies analyzed per sample
burden = read.table("../data/meta_analysis_10samples_allSDS8_buccal_Dec2021_SDS/burden_by_sample.txt", header=T, stringsAsFactors = F, sep="\t")
ncolonies = data.frame(tapply(burden$colony, INDEX=burden$internal_id, FUN=function(X) length(unique(X))))
colnames(ncolonies) = "n.colonies"
ncolonies$ID = rownames(ncolonies)
# SDS1 SDS10  SDS2  SDS3  SDS4  SDS5  SDS6  SDS7  SDS8  SDS9 
# 18    27    38    23    18   100    27    35    10    27 

## genotype data
genotypes = read.csv("../data/germlineSBDS_manualcuration_10samples.csv", header=T, stringsAsFactors = F)
genotypes2 = merge(genotypes, ncolonies, by.x="Pt.ID", by.y="ID")
genotypes2$Pt.ID = factor(genotypes2$Pt.ID, levels=genotypes2$Pt.ID[c(1,3:10,2)])
genotypes_melt = melt(genotypes2, id.vars=colnames(genotypes2)[c(1:6,9:10)] )


# splice code a pink shades, coding as blue shades
ggplot(genotypes_melt, aes(x=Pt.ID, y=Age.y.1, group=variable, fill=value))+
  geom_col(position=position_dodge2(padding=0), width=0.5)+
  theme_bw()+
  scale_fill_manual(values=c( "pink", "light blue4", "light blue1", "light blue3") )+
  xlab("Individual")+
  ylab("Age (years)")+
  #geom_label(aes(x=Pt.ID, y=Age.y.1+1, label=n.colonies), colour = "black",fill=NA)+
  geom_label(aes(x=Pt.ID, y=-1, label=n.colonies), colour = "black",fill=NA, size=3)+
  theme()
ggsave("cohort_ages_genotypes_barplot.pdf", width=7, height=3)



## hsc phenotyping
phen = read.csv("../data/Frozen_donor_StemCellTech_w.SDSpatients_18-02-2022.csv", stringsAsFactors = F, header=T)
phen$Patient.ID = revalue(phen$Patient.ID, replace=c("D1"="N1","D2"="N2","D3"="N3"))
phen$Patient.ID = factor(phen$Patient.ID, levels=sort(unique(phen$Patient.ID))[c(1:4,6:11,5,12)])
phen$Patient_Donor = revalue(phen$Patient_Donor, replace=c("HD"="normal"))
phen$Patient_Donor = factor(phen$Patient_Donor, levels=c("SDS","normal"))
mymean = data.frame(tapply(phen$CD34_pct_of_live, INDEX=phen$Patient_Donor, FUN=mean))
colnames(mymean) = "mean"
mymean$group = c(1,2)


ggplot(subset(phen, Patient.ID!="unknown"), aes(x=Patient_Donor, y=CD34_pct_of_live, col=Patient.ID, shape=Sample.type, size=Sample.type) )+
  #geom_jitter(aes(x=Patient_Donor, y=CD34_pct_of_live, col=Patient.ID, shape=Sample.type, size=Sample.type), width=0.2)+
  geom_point(aes(x=Patient_Donor, y=CD34_pct_of_live, col=Patient.ID, shape=Sample.type, size=Sample.type), position=position_dodge2(width=0.2))+
  scale_color_manual(values=c("grey1","grey2","grey3","red","blue","green","pink","purple","cyan","light blue","dark red"))+
  scale_y_continuous(trans='log10')+
  scale_shape_manual(values=c(19,2))+
  scale_size_manual(values=c(2,4))+
  geom_errorbar(aes(xmin=group-0.3, xmax=group+0.3,  y=mean), data=mymean,
               inherit.aes=FALSE, width=0, linetype=1)+
  xlab("")+
  ylab("CD34+ (% MNCs)")+
  theme_bw()+
  annotation_logticks(sides = "l")+
  theme(legend.text = element_text(size=7),
#        legend.key.size = element_rect(6)
        )
ggsave("cd34percent_scatter.pdf", width=3, height=3)



## bone marrow samples only (peripheral blood lower)
ggplot(subset(phen, Patient.ID!="unknown" & Sample.type=="BM"), aes(x=Patient_Donor, y=CD34_pct_of_live, fill=Patient.ID, size=Patient_Donor) )+
  geom_point(position=position_dodge2(width=0.4), shape=23)+
  scale_fill_manual(values=c("grey1","grey2","grey3","red","blue","green","pink","purple","cyan","light blue","dark red"))+
  geom_errorbar(aes(xmin=group-0.4, xmax=group+0.4,  y=mean), data=mymean,
                inherit.aes=FALSE, width=0, linetype=1)+
  xlab("")+
  ylab("CD34+ (% MNCs)")+
  scale_size_manual(values=c(3,2))+
  theme_bw()+
  theme(legend.text = element_text(size=7),
        #        legend.key.size = element_rect(6)
  )
ggsave("cd34percent_BMonly_scatter.pdf", width=3, height=3)
