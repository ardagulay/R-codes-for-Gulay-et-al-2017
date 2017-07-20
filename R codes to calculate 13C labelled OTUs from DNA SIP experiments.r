## R codes to calculate 13C labelled OTUs from DNA SIP experiments
## Written by Dr. Arda Gülay (2017), DTU @ Technical University of Denmark.

library(reshape2)

## OTUtable_normalized is the file acquired after applying R codes for sample normalization "R codes to normalize sequence libraries from SIP experiments"

OTUtable_normalized->b_clean
melt(b_clean)->b_clean

##Then, modify the column names as follows : "OTU" "variable" "value"
colnames(b_clean) <- c("OTU", "variable","value")

coefficent<-0
coefficent_list<-c()

for(j in 1:length(unique(b_clean$OTU)))
{
    coefficent<-0
    b_clean[j,1]->OTUname
    paste("^",OTUname,"$", sep="")->OTUname
    b_clean[grep(OTUname, b_clean$OTU, ignore.case=F),]->OTU_grep
    for(k in 1:length(row.names(samples)))
    { 
        samples[k,]->samplename
        OTU_grep[sort(grep(samplename, OTU_grep$variable, ignore.case=F)),]->OTU_grep_selected
        coefficent<-0
		total_density<-0
			for(f in 1:length(row.names(OTU_grep_selected)))
			{ 
            OTU_grep_selected$variable[f]->fraction_name
            paste("^",fraction_name,"$", sep="")->fraction_name
            OTU_grep_selected[grep(fraction_name, OTU_grep_selected$variable, ignore.case=F),]->fraction_line
            DNA_density[grep(fraction_name, DNA_density$Sample, ignore.case=F),]->density_dna_line
			coefficent<-coefficent+(as.double(as.character(density_dna_line$Density))*as.double(as.character(fraction_line$value)))
			total_density<-coefficent+(as.double(as.character(density_dna_line$Density)))
			}
        coefficent_list<-rbind(coefficent_list,cbind(as.character(fraction_line$OTU), as.character(samplename), total_density, as.double(as.character(coefficent)), as.double(as.character(((as.double(as.character(coefficent)))/as.double(as.character(total_density)))))))
  	}
}
as.data.frame(coefficent_list)->coefficent_list
colnames(coefficent_list) <- c("OTU", "sample","Yij","Xjk*Yijk", "Wij")
write.table(coefficent_list,"coefficent_list_DNA_manfield_99_phase1.txt", sep="\t")

criteria1_list<-c()
for(j in 1:length(unique(coefficent_list$V1)))
{
   unique(coefficent_list$V1)[j]->OTUname
   paste("^",OTUname,"$", sep="")->OTUname
   coefficent_list[grep(OTUname, coefficent_list$V1, ignore.case=F),]->OTU_grep
   OTU_grep[grep("DNA.R1.C1", OTU_grep$V2, ignore.case=F),]->DNA.R1.C1
   OTU_grep[grep("DNA.R1.C2", OTU_grep$V2, ignore.case=F),]->DNA.R1.C2
   OTU_grep[grep("DNA.R1.C3", OTU_grep$V2, ignore.case=F),]->DNA.R1.C3
   OTU_grep[grep("DNA.R1.C4", OTU_grep$V2, ignore.case=F),]->DNA.R1.C4
   OTU_grep[grep("DNA.R2.C1", OTU_grep$V2, ignore.case=F),]->DNA.R2.C1
   OTU_grep[grep("DNA.R2.C2", OTU_grep$V2, ignore.case=F),]->DNA.R2.C2
   OTU_grep[grep("DNA.R2.C3", OTU_grep$V2, ignore.case=F),]->DNA.R2.C3
   OTU_grep[grep("DNA.R2.C4", OTU_grep$V2, ignore.case=F),]->DNA.R2.C4

   
   suppressWarnings(as.double(as.character(DNA.R1.C1$V3))/(as.double(as.character(DNA.R1.C2$V3))+1))->DNA_C1vsC2
   suppressWarnings(as.double(as.character(DNA.R1.C4$V3))/(as.double(as.character(DNA.R1.C3$V3))+1))->DNA_C4vsC3
   suppressWarnings(as.double(as.character(DNA.R2.C1$V3))/(as.double(as.character(DNA.R2.C2$V3))+1))->DNA_C5vsC6
   suppressWarnings(as.double(as.character(DNA.R2.C4$V3))/(as.double(as.character(DNA.R2.C3$V3))+1))->DNA_C8vsC7
   if (identical(DNA_C1vsC2, numeric(0))==TRUE){DNA_C1vsC2<-0}
   if (identical(DNA_C4vsC3, numeric(0))==TRUE){DNA_C4vsC3<-0}
   if (identical(DNA_C5vsC6, numeric(0))==TRUE){DNA_C5vsC6<-0}
   if (identical(DNA_C8vsC7, numeric(0))==TRUE){DNA_C8vsC7<-0}
   criteria1_list<-rbind(criteria1_list,cbind(as.character(unique(coefficent_list$V1)[j]),"DN A_C1vsC2",as.double(as.character(DNA_C1vsC2))))
   criteria1_list<-rbind(criteria1_list,cbind(as.character(unique(coefficent_list$V1)[j]),"DNA_C4vsC3",as.double(as.character(DNA_C4vsC3))))
   criteria1_list<-rbind(criteria1_list,cbind(as.character(unique(coefficent_list$V1)[j]),"DNA_C5vsC6",as.double(as.character(DNA_C5vsC6))))
   criteria1_list<-rbind(criteria1_list,cbind(as.character(unique(coefficent_list$V1)[j]),"DNA_C8vsC7",as.double(as.character(DNA_C8vsC7))))
}
as.data.frame(criteria1_list)->DNA_SIP_13C_OTUs
colnames(DNA_SIP_13C_OTUs) <- c("OTU", "Sample","BD_shift")
write.table(DNA_SIP_13C_OTUs,"DNA_SIP_13C_OTUs.txt", sep="\t")




