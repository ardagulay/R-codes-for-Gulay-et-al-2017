library(data.table)
library(scatterplot3d)

##CsCl 
##Light range :  until 1.73 DNA de if amount in heavy col labelled / amount in light col non-labelled
##Heavy range: after 1.73
##CsTFA 
##light range: until 1.80
##Heavy range after 1.80

##If the data is coming from the data normalization, it uses the otutable() function of the phyloseq:
##Therefore melt using (reshape2):
##OTUtable_normalized is the file acquired after applying R codes for sample normalization "R codes to normalize sequence libraries from SIP experiments"

library(reshape2)
OTUtable_normalized->b_clean
melt(b_clean)->b_clean

##Then, modify the column names as follows : OTU variable value
colnames(b_clean) <- c("OTU", "variable","value")

#YHL=DNA/RNA concentration in the HL fractions
#XHL= weights of representative heavy fractions of the labeled pulse
# YLL=DNA/RNA concentration in the LL fractions
#XLL=weights of representative light fractions of the labeled pulse
#BD= the mean buoyant density of the ith phylo-type in the labeled incubation

#sigma= is the standard deviation of the RNA distribution derived from the RNA concentrations across the gradient
#Sigma has to be calculated beforehand (see Zemb et al., 2012) and the values should be used as a input here.

sigma <- c(0.000758664, 0.000758664, 0.001939269,0.001939269,0.001684653,0.001684653,0.000738739,0.000738739)
names(sigma) <- c("R1.C1.F", "R1.C2.F", "R1.C3.F","R1.C4.F","R2.C1.F","R2.C2.F","R2.C3.F","R2.C4.F")

##Determine the labelled and unlabelled runs

labelled<-c("DNA.R1.C1","DNA.R1.C4","DNA.R2.C1","DNA.R2.C4","R1.C1.F","R1.C4.F","R2.C1.F","R2.C4.F")

##non_labelled<-c("DNA.R1.C2","DNA.R1.C3","DNA.R2.C2","DNA.R2.C3","R1.C2.F","R1.C3.F","R2.C2.F","R2.C3.F")

coefficent<-0
coefficent_list_manfield<-c()
for(j in 1:length(unique(b_clean$OTU)))
{
    coefficent<-0
    b_clean[j,1]->OTUname
    paste("^",OTUname,"$", sep="")->OTUname
    b_clean[grep(OTUname, b_clean$OTU, ignore.case=F),]->OTU_grep
    print(OTUname)
    for(k in 1:length(row.names(samples)))
    { 
        samples[k,]->samplename
        OTU_grep[sort(grep(samplename, OTU_grep$variable, ignore.case=F)),]->OTU_grep_selected
        if (samplename%in%labelled)
        {
            if(grepl("DNA", samplename)){
             next
            } else{
                DNA_density[grep(samplename, DNA_density$Sample, ignore.case=F),]->samples_selected 
                samples_selected[which(samples_selected$Density>=1.80),]->samples_selected_high
                samples_selected[which(samples_selected$Density<1.80),]->samples_selected_low
                OTU_grep_selected[grepl(paste(paste("^",samples_selected_high$Sample,"$", sep=""),collapse="|"),OTU_grep_selected$variable, ignore.case=F),]->heavy_samples_RNA
                OTU_grep_selected[grepl(paste(paste("^",samples_selected_low$Sample,"$", sep=""),collapse="|"),OTU_grep_selected$variable, ignore.case=F),]->light_samples_RNA
                XHL<-samples_selected[which(samples_selected$Sample==(as.character((heavy_samples_RNA[which(heavy_samples_RNA$value==max(heavy_samples_RNA$value)),])$variable)[1])),]$Density
                XLL<-samples_selected[which(samples_selected$Sample==(as.character((light_samples_RNA[which(light_samples_RNA$value==max(light_samples_RNA$value)),])$variable)[1])),]$Density
                YHL<-as.double(as.character((heavy_samples_RNA[which(heavy_samples_RNA$value==max(heavy_samples_RNA$value)),])$value)[1])
                YLL<-as.double(as.character((light_samples_RNA[which(light_samples_RNA$value==max(light_samples_RNA$value)),])$value)[1])
					if (YLL==0&YHL==0){next} else if (YLL==0){YLL<-1E-20} else if (YHL==0){YHL<-1E-20}
						if (samplename%in%labelled){Control<-"TRUE"} else {Control<-"FALSE"}
                BD<-abs(((2*((sigma[[as.character(samplename)]])^2)*(log((YLL)/(YHL))))-((XHL)^2)+((XLL)^2))/(2*(XHL-XLL)))
				coefficent_list_manfield<-rbind(coefficent_list_manfield,cbind(as.character(OTU_grep$OTU[1]), as.character(samplename), YHL,XHL,YLL,XLL,BD,Control))
                print(rbind(coefficent_list_manfield,cbind(as.character(OTU_grep$OTU[1]), as.character(samplename), YHL,XHL,YLL,XLL,BD,Control)))
            }} else {
                if(grepl("DNA", samplename)){
                next
                } else {
                DNA_density[grep(samplename, DNA_density$Sample, ignore.case=F),]->samples_selected 
                samples_selected[which(samples_selected$Density>=1.80),]->samples_selected_high
                samples_selected[which(samples_selected$Density<1.80),]->samples_selected_low
                OTU_grep_selected[grepl(paste(paste("^",samples_selected_high$Sample,"$", sep=""),collapse="|"),OTU_grep_selected$variable, ignore.case=F),]->heavy_samples_RNA
                OTU_grep_selected[grepl(paste(paste("^",samples_selected_low$Sample,"$", sep=""),collapse="|"),OTU_grep_selected$variable, ignore.case=F),]->light_samples_RNA
                XHL<-samples_selected[which(samples_selected$Sample==(as.character((heavy_samples_RNA[which(heavy_samples_RNA$value==max(heavy_samples_RNA$value)),])$variable)[1])),]$Density
                XLL<-samples_selected[which(samples_selected$Sample==(as.character((light_samples_RNA[which(light_samples_RNA$value==max(light_samples_RNA$value)),])$variable)[1])),]$Density
                YHL<-as.double(as.character((heavy_samples_RNA[which(heavy_samples_RNA$value==max(heavy_samples_RNA$value)),])$value)[1])
                YLL<-as.double(as.character((light_samples_RNA[which(light_samples_RNA$value==max(light_samples_RNA$value)),])$value)[1])
					if (YLL==0&YHL==0){next} else if (YLL==0){YLL<-1E-20} else if (YHL==0){YHL<-1E-20}
                BD<-abs(((2*((sigma[[as.character(samplename)]])^2)*(log((YLL)/(YHL))))-((XHL)^2)+((XLL)^2))/(2*(XHL-XLL)))
						if (samplename%in%labelled){Control<-"TRUE"} else {Control<-"FALSE"}
				coefficent_list_manfield<-rbind(coefficent_list_manfield,cbind(as.character(OTU_grep$OTU[1]), as.character(samplename), YHL,XHL,YLL,XLL,BD,Control))
                print(rbind(coefficent_list_manfield,cbind(as.character(OTU_grep$OTU[1]), as.character(samplename), YHL,XHL,YLL,XLL,BD,Control)))
                }}}}
                as.data.frame(coefficent_list_manfield)-> coefficent_list_manfield
                colnames(coefficent_list_manfield) <- c("OTU", "sample","YHL","XHL", "YLL","XLL","BD","IS_CONTROL")
                write.table(coefficent_list_manfield,"coefficent_list_manfield_RNA_24.05.17_3mil.99.n10.json.txt", sep="\t")


criteria_manfield_list<-c()
for(j in unique(coefficent_list_manfield$OTU))
{
   paste("^",j,"$", sep="")->OTUname
   coefficent_list_manfield[grep(OTUname, coefficent_list_manfield$OTU, ignore.case=F),]->OTU_grep

	if("R1.C1.F"%in% OTU_grep$sample&"R1.C2.F"%in% OTU_grep$sample=="TRUE"){
   R1.C1vsC2<-suppressWarnings(as.double(as.character((OTU_grep[grep("R1.C1.F", OTU_grep$sample, ignore.case=F),])$BD))-as.double(as.character((OTU_grep[grep("R1.C2.F", OTU_grep$sample, ignore.case=F),])$BD)))
   criteria_manfield_list<-rbind(criteria_manfield_list,cbind(j,"R1.C1vsC2",as.double(as.character((OTU_grep[grep("R1.C2.F", OTU_grep$sample, ignore.case=F),])$BD)),as.double(as.character((OTU_grep[grep("R1.C1.F", OTU_grep$sample, ignore.case=F),])$BD)),as.double(as.character(R1.C1vsC2))))
																			}
		if("R1.C3.F"%in% OTU_grep$sample&"R1.C4.F"%in% OTU_grep$sample=="TRUE"){
   R1.C4vsC3<-suppressWarnings(as.double(as.character((OTU_grep[grep("R1.C4.F", OTU_grep$sample, ignore.case=F),])$BD))-as.double(as.character((OTU_grep[grep("R1.C3.F", OTU_grep$sample, ignore.case=F),])$BD))) 
   criteria_manfield_list<-rbind(criteria_manfield_list,cbind(j,"R1.C4vsC3",as.double(as.character((OTU_grep[grep("R1.C4.F", OTU_grep$sample, ignore.case=F),])$BD)),as.double(as.character((OTU_grep[grep("R1.C3.F", OTU_grep$sample, ignore.case=F),])$BD)),as.double(as.character(R1.C4vsC3))))
																				}

			if("R2.C1.F"%in% OTU_grep$sample&"R2.C2.F"%in% OTU_grep$sample=="TRUE"){
   R2.C1vsC2<-suppressWarnings(as.double(as.character((OTU_grep[grep("R2.C1.F", OTU_grep$sample, ignore.case=F),])$BD))-as.double(as.character((OTU_grep[grep("R2.C2.F", OTU_grep$sample, ignore.case=F),])$BD))) 
   criteria_manfield_list<-rbind(criteria_manfield_list,cbind(j,"R2.C1vsC2",as.double(as.character((OTU_grep[grep("R2.C1.F", OTU_grep$sample, ignore.case=F),])$BD)),as.double(as.character((OTU_grep[grep("R2.C2.F", OTU_grep$sample, ignore.case=F),])$BD)),as.double(as.character(R2.C1vsC2))))
																					}

				if("R2.C4.F"%in% OTU_grep$sample&"R2.C3.F"%in% OTU_grep$sample=="TRUE"){
   R2.C4vsC3<-suppressWarnings(as.double(as.character((OTU_grep[grep("R2.C4.F", OTU_grep$sample, ignore.case=F),])$BD))-as.double(as.character((OTU_grep[grep("R2.C3.F", OTU_grep$sample, ignore.case=F),])$BD))) 
   criteria_manfield_list<-rbind(criteria_manfield_list,cbind(j,"R2.C4vsC3",as.double(as.character((OTU_grep[grep("R2.C4.F", OTU_grep$sample, ignore.case=F),])$BD)),as.double(as.character((OTU_grep[grep("R2.C3.F", OTU_grep$sample, ignore.case=F),])$BD)),as.double(as.character(R2.C4vsC3))))
																						}
} 

as.data.frame(criteria_manfield_list)->criteria_manfield_list
colnames(criteria_manfield_list) <- c("OTU", "Sample","WLAB","WLİGHT","BD_shift")
criteria_manfield_list[ which(as.double(as.character(criteria_manfield_list$BD_shift))>0),]->RNA_SIP_13C_OTUs
write.table(RNA_SIP_13C_OTUs,"RNA_SIP_13C_OTUs.txt", sep="\t")



## If you would like to plot BD shifts as a funciton of OTUs in different samples, you can use the below codes.
##-------------------------------------------------------------------------------------------------
##par(mfrow=c(1, 4))
##for(k in sort(as.character(unique(RNA_SIP_13C_OTUs$Sample))))
 ##{
    ##RNA_SIP_13C_OTUs[grep(k, RNA_SIP_13C_OTUs$Sample, ignore.case=F),]->sample_grep
    ##setorder(RNA_SIP_13C_OTUs, -BD_shift)
    ##c1<-as.double(as.character(sample_grep$BD_shift))
    ##otu<-as.character(sample_grep$OTU)
   ##plot(c1, ylab="BD_shift",xlab="OTUs",main=sample_grep$Sample[1]) ##+text(c1, labels=otu,pos=4,cex=0.9)
   ##}
