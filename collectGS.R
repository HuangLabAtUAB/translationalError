#collect the golden standard peptides, which are the ones searched by using mouse+house reference
#mq = maxQuant
#combined means human and mouse combined fasta reference.
library(data.table)

####### maxquant conventional search mouse human combined for all samples ###
mqCombined = fread('./mqCombined/evidence.txt')
colnames(mqCombined) = gsub(colnames(mqCombined), pattern='\\s+', replacement='_')
length(unique(mqCombined$Raw_file)) #24 raw files
dim(mqCombined) #913,616

##spectral count uniq-human vs uniq-mouse
mqCombined[, species:=ifelse(grepl(Proteins, pattern='taxon\\|9606'),'human','mouse')]
mqCombined[, species:=ifelse(grepl(Proteins, pattern='taxon\\|10090.*taxon\\|9606'),'both',species)]
mqCombined[, species:=ifelse(grepl(Proteins, pattern='taxon\\|9606.*taxon\\|10090'),'both',species)]
mqCombined[, msSample:=gsub(Raw_file, pattern='180316(_[0-9]+_).*', replacement='\\1')]
table(mqCombined$species)
#both  human  mouse 
#438,576 283,288 191,752 
##spectral count
ggplot(mqCombined, aes(msSample,fill=species))+
  geom_bar(stat = 'count')+
  theme(axis.text.x = element_text(angle = -90))
mqCombinedSptrCount = mqCombined[,.N, by=c('species','msSample')]
mqCombinedSptrCount[,perc:=N/sum(N),by='msSample']

ggplot(mqCombinedSptrCount, aes(species,perc))+
  geom_boxplot(outliers = F)+geom_jitter(width = 0.05)+
  theme(axis.text.x = element_text(angle = -90))+
  geom_line(aes(group=msSample),col='grey',alpha=0.5)
range(mqCombinedSptrCount[species=='mouse',perc])#0.1000719 0.3144621
median(mqCombinedSptrCount[species=='mouse',perc])#0.2317519

ggplot(mqCombinedSptrCount, aes(msSample, N,fill = species))+
  geom_bar(stat = 'identity',position = 'fill')+
  scale_y_continuous(labels = scales::percent) +  # Convert y-axis to percentage format
  labs(y = "Percentage (%)", x = "Category")
##
#need to focus only on 'discoverable' mouse peptides, i.e., have single SAV
#ground truth in the spectral level:put spectrum ID, seq, and sample together as column
mqMouseDt = mqCombined[species=='mouse'&Sequence%in%hmPairs$mousePep,
                         .(msmsScanNum = `MS/MS_scan_number`, 
                                           msmsScanNums = `MS/MS_scan_numbers`,
                                           msSample,
                                           Sequence, Modified_sequence, 
                                          PEP, Score, Retention_time)]
table(grepl(mqMouseDt$msmsScanNums, pattern=';')) #48,739  T4625 #one row has more than one scans
##for very few raw case, mq report two sequence for a ms/ms scan
mqMouseDt[,table(duplicated(msmsScanNum)),by='msSample']

mqMouseDt = mqMouseDt[, .(msmsScanNum = unlist(tstrsplit(msmsScanNums, ";", type.convert = TRUE))), 
                      by=c('msSample','Sequence','msmsScanNums',
                           'Modified_sequence', 'PEP', 'Score','Retention_time')] #process ; rows
dim(mqMouseDt) #58539      8
mqMouseDt[,(sum(duplicated(msmsScanNum))/.N),by='msSample'] #very few, all less than 0.01
mqMouseDt = mqMouseDt[,.(Sequence=Sequence[which.min(PEP)],
                         Modified_sequence=Modified_sequence[which.min(PEP)],
                         PEP=PEP[which.min(PEP)],
                         Score=Score[which.min(PEP)],
                         msmsScanNums=msmsScanNums[which.min(PEP)],
                         Retention_time=Retention_time[which.min(PEP)]) ,
                      by='msSample,msmsScanNum'] 
dim(mqMouseDt) #58429 
mqMouseDt[,(sum(duplicated(msmsScanNum))/.N),by='msSample']  #all zero now
#need to add sav info to mqMouseDt
table(duplicated(hmPairsShort$humanPep)) #284186   478 

hmPairs_uniqMouse = hmPairsShort[,.(humanPep=paste0(humanPep, collapse = ';'),
                                    savIdx = paste0(savIdx, collapse = ';'),
                                    humanAA = paste0(humanAA, collapse = ';'),
                                    endAA = paste0(endAA, collapse = ';'),
                                    sav = paste0(sav, collapse = ';')),
                                 by='mousePep']

table(duplicated(hmPairs_uniqMouse$mousePep)) #all F
hmPairs_uniqMouse[grepl(savIdx,pattern=';')][1:5]



mqMouseDt = cbind(mqMouseDt,
                  hmPairs_uniqMouse[mqMouseDt$Sequence,.(humanPep,mousePep,savIdx, humanAA, endAA, sav),
                               on='mousePep'])
table(mqMouseDt$Sequence == mqMouseDt$mousePep) #all T
mqMouseDt$mousePep=NULL #duplicated column with "Sequence"


dim(mqMouseDt) #58429    
mqMouseDt = tidyr::separate_rows(mqMouseDt, humanPep, savIdx, humanAA,
                          endAA,
                          sav, sep = ';')#process ; rows
dim(mqMouseDt) #58762 
mqMouseDt = as.data.table(mqMouseDt)

mqMouseDt[, spectrumID := paste0(msSample,msmsScanNum)]
mqMouseDt[,table(substr(Sequence,savIdx, savIdx)==endAA)]
mqMouseDt[endAA!='X',table(substr(Sequence,savIdx, savIdx)==endAA)] #all T
s=mqMouseDt$Sequence
for(idx in 1:nrow(mqMouseDt)){
  if(mqMouseDt[idx,endAA]=='X')
    substr(s[idx],mqMouseDt[idx,savIdx],mqMouseDt[idx,savIdx]) = 'X'
}
mqMouseDt$seq4Eva=s
mqMouseDt[, spectrumSeq:=paste0(spectrumID,'_',seq4Eva)]

###
mqMouse = split(mqMouseDt, f = mqMouseDt$msSample)
#the theoretical human counterpart is detected or not in the sample ms sample?
#for dp, only the detected ones can be the GS
for(ms_sample in names(mqMouse)){
  mqMouse[[ms_sample]]$humanDetected = ifelse(mqMouse[[ms_sample]]$humanPep%in%
                                                 mqCombined[msSample==ms_sample,Sequence],
                                              'Yes','No')
}
lapply(mqMouse, function(x)table(x$humanDetected))
mqMouseDt = rbindlist(mqMouse)


#########fragpipe conventional search mouse human combined for all samples ###
fpCombinedList = list()

#fpCombined = fread('./fpcombined/psm.tsv')
fpCombined = fread('/Volumes/ChenDrive/fpSearchMSNew_pancreatic/fpCombine/psm.tsv')

fpPsmCount = fpCombined[,.N,by='Spectrum File']
fwrite(fpPsmCount, file = 'fpcombined/fpPsmCount.tsv')
# fpCombined = fread('./fpcombined/psm_trypsin.tsv')
fpCombined[,grep(Peptide, pattern='^P',value = T)]
fpCombined[,table(grepl(Peptide, pattern='^P'))] #0.007482805 ->#0.001553727??
colnames(fpCombined) = gsub(colnames(fpCombined), 
                            pattern='\\s+', 
                            replacement = '_')
dim(fpCombined) #2,312,572 39 ->2307717

fpCombined[, species:=ifelse(grepl(Mapped_Proteins, pattern='taxon\\|9606'),'human','mouse')]
fpCombined[, species:=ifelse(grepl(Mapped_Proteins, pattern='taxon\\|10090.*taxon\\|9606'),'both',species)]
fpCombined[, species:=ifelse(grepl(Mapped_Proteins, pattern='taxon\\|9606.*taxon\\|10090'),'both',species)]
fpCombined[, msSample:=gsub(Spectrum, pattern='^[0-9]+(_[0-9]{2}_).*', replacement='\\1')]
fpCombined = fpCombined[as.numeric(substring(msSample,2,3)) < 25] #only need the first twenty-four samples
fpCombined[,table(species)]
fpCombined[,msmsScanNums := gsub(Spectrum, pattern=".*?\\.([0-9]+)\\..*", replacement = '\\1')]
fpCombined[, msmsScanNum := as.numeric(msmsScanNums)]
# both  human  mouse 
#395,502 296,659 469,956 
#->394196 293985 471574 

##spectral count for each sample
ggplot(fpCombined, aes(msSample,fill=species))+
  geom_bar(stat = 'count')+
  theme(axis.text.x = element_text(angle = -90))
fpCombinedSptrCount = fpCombined[,.N, by=c('species','msSample')]
fpCombinedSptrCount[,perc:=N/sum(N),by='msSample']

ggplot(fpCombinedSptrCount, aes(species,perc))+
  geom_boxplot(outliers = F)+geom_jitter(width = 0.05)+
  theme(axis.text.x = element_text(angle = -90))
##
#ground truth in the spectral level
table(fpCombined[species=='mouse',unique(Peptide) %in% hmPairs$mousePep])
# FALSE  TRUE 
# 29130  5943 -> 28884  5947 
# table(fpCombined[species=='both',unique(Peptide) %in% hmPairs$mousePep])
# table(hmPairs$mousePep%in%fpCombined[species=='mouse',unique(Peptide)])

fpMouseDt = fpCombined[species=='mouse'&Peptide%in%hmPairs$mousePep, 
                       .(msmsScanNums, msmsScanNum,msSample,Peptide,Expectation,Hyperscore,
                         Modified_Peptide,Assigned_Modifications,Retention)]
table(grepl(fpMouseDt$msmsScanNums, pattern=';')) #all F

fpMouseDt[, spectrumID := paste0(msSample,msmsScanNum)]

fpMouseDt = cbind(fpMouseDt, hmPairs_uniqMouse[fpMouseDt$Peptide,.(humanPep,mousePep,savIdx, humanAA, endAA, sav),
                                               on='mousePep'])


dim(fpMouseDt) #67760  -> 68250
fpMouseDt = tidyr::separate_rows(fpMouseDt, humanPep, savIdx, humanAA,
                                 endAA,
                                 sav,sep = ';')
dim(fpMouseDt) #68129 ->68621
table(fpMouseDt$Peptide == fpMouseDt$mousePep)
fpMouseDt$mousePep = NULL
fpMouseDt = as.data.table(fpMouseDt)

fpMouseDt[endAA!='X',table(substr(Peptide,savIdx, savIdx)==endAA)] #all T
s=fpMouseDt$Peptide
for(idx in 1:nrow(fpMouseDt)){
  if(fpMouseDt[idx,endAA]=='X')
    substr(s[idx],fpMouseDt[idx,savIdx],fpMouseDt[idx,savIdx]) = 'X'
}
fpMouseDt$seq4Eva=s
fpMouseDt[, spectrumSeq:=paste0(spectrumID,'_',seq4Eva)]

tmp= hmPairs[!duplicated(mousePep)]
fpMouseDt[,sav:=tmp[fpMouseDt$Peptide, sav,on='mousePep']]
rm(tmp)

fpMouse = split(fpMouseDt, f = fpMouseDt$msSample)
#the theoretical human counterpart is detected or not in the sample ms sample?
for(ms_sample in names(fpMouse)){
  fpMouse[[ms_sample]]$humanDetected = ifelse(fpMouse[[ms_sample]]$humanPep%in%
                                                fpCombined[msSample==ms_sample,Peptide],
                                              'Yes','No')
}
#lapply(fpMouse, function(x)table(x$humanDetected))

fpMouseDt = rbindlist(fpMouse)

fpMouseDt[spectrumSeq=='_01_94_QATENAKEEVKR'] #the retention time is second
mqMouseDt[spectrumSeq=='_01_94_QATENAKEEVKR'] #the retention time is minute
rm(fpCombined)

#the overlap between fp and mq search results in the spectrum level
dim(fpMouseDt) #68129 ->68621
dim(mqMouseDt) #58762
table(fpMouseDt$spectrumSeq %in% mqMouseDt$spectrumSeq) 
#F17507 T51114 pos rate 74.5%

#17,002 51,127 -> 17,311 TRUE51,310 
table(fpMouseDt$spectrumID %in% mqMouseDt$spectrumID) #17311 51310 

table(mqMouseDt$spectrumSeq %in% fpMouseDt$spectrumSeq) #7648 T51114 87.0%
table(mqMouseDt$spectrumID %in% fpMouseDt$spectrumID) #7451 T51,311 

##the overlap between fp and mq search results in the peptide level
all(names(fpMouse) == names(mqMouse)) #TRUE
sapply(names(fpMouse), function(x)table(unique(fpMouse[[x]]$Peptide) %in% mqMouse[[x]]$Sequence))
sapply(names(fpMouse), function(x)table(unique(mqMouse[[x]]$Sequence) %in% fpMouse[[x]]$Peptide))



##plot out the overlap
commonSearchOverlap = list()
for(ms_names in names(fpMouse)){
  commonSearchOverlap[[ms_names]] = data.table(fp_in_mq_sp = as.numeric(table(fpMouse[[ms_names]]$spectrumSeq %in% 
                                            mqMouse[[ms_names]]$spectrumSeq)),
                                            mp_in_fq_sp = as.numeric(table(mqMouse[[ms_names]]$spectrumSeq %in% 
                                                    fpMouse[[ms_names]]$spectrumSeq)),
                                            fp_in_mq_pep = as.numeric(table(unique(fpMouse[[ms_names]]$Peptide) %in% 
                                                                   mqMouse[[ms_names]]$Sequence)),
                                            mp_in_fq_pep = as.numeric(table(unique(mqMouse[[ms_names]]$Sequence) %in% 
                                                                  fpMouse[[ms_names]]$Peptide))
                                            )
                                            
}
commonSearchOverlap = rbindlist(commonSearchOverlap, idcol = 'sample_name')
commonSearchOverlap[, type:= rep(c('Fa','Tr'),length(fpMouse))]
commonSearchOverlap = data.table::melt(commonSearchOverlap, 
                                       id.vars = c('sample_name','type'),
                                       variable.name = 'comparison', 
                                       value.name = 'count')
commonSearchOverlap[, groupByType := gsub(comparison, pattern='.*_', replacement='')]
commonSearchOverlap[, uniq:=gsub(comparison, pattern='_.*', replacement='')]
commonSearchOverlap[, type:=ifelse(type=='Tr', 'overlap', uniq)]
commonSearchOverlap = commonSearchOverlap[,.(count=count[1]),by=c('sample_name','type','groupByType')]
commonSearchOverlap[, perc:=count/sum(count), by=c('sample_name','groupByType')]

ggplot(commonSearchOverlap, aes(sample_name,count,fill=type))+
  geom_bar(stat = 'identity', position='fill')+
  facet_wrap(.~groupByType)+
  theme(axis.text.x = element_text(angle = -90))

ggplot(commonSearchOverlap[groupByType=='sp'], aes(type,perc))+
  geom_boxplot(outliers = F)+geom_jitter(width = 0.05)+
  #facet_wrap(.~groupByType)+
  scale_y_continuous(labels = scales::percent,limits = c(0,1))+
  theme(axis.text.x = element_text(angle = -90))
  
median(commonSearchOverlap[type=='overlap'&groupByType=='sp',perc])  
  
###collect good quality peptides for pepquery (autoRT!!)
mqCombined[,spSpectSeq:=paste0(msSample,`MS/MS_scan_number`,'_',Sequence)]
fpCombined[, spSpectSeq:=paste0(msSample,msmsScanNum,'_',Peptide)]
mqCombined[,spSpect:=paste0(msSample,`MS/MS_scan_number`)]
fpCombined[, spSpect:=paste0(msSample,msmsScanNum)]

table(mqCombined[species=='human',spSpectSeq]%in%fpCombined[species=='human',spSpectSeq]) 
#FALSE   TRUE 
#121,220 162,068 ->123032 160256 
table(fpCombined[species=='human',spSpectSeq]%in%mqCombined[species=='human',spSpectSeq]) 
# FALSE   TRUE 
# 134,591 162,068 -> 133729 160256 
table(mqCombined[species=='human',spSpect]%in%fpCombined[species=='human',spSpect]) 
#FALSE   TRUE 
#117547 165741  ->119335 163953 
table(fpCombined[species=='human',spSpect]%in%mqCombined[species=='human',spSpect]) 
# FALSE   TRUE 
# 132290 164369 -> 131388 162597 
#as long as the spectrum can be identified, 
#the sequence are mostly the same between fp and mq
fpMqHmCommon = intersect(fpCombined[species=='human',spSpectSeq],
                         mqCombined[species=='human',spSpectSeq]) 

fpMqHmCommon = fpCombined[species=='human'&spSpectSeq%in%fpMqHmCommon,
                          .(msSample, Peptide,Modified_Peptide,Assigned_Modifications,Retention)]
table(fpMqHmCommon$Assigned_Modifications)

fpMqHmCommon[, pepQueryFormat:=forMatPep(Peptide,Assigned_Modifications)]
table(fpMqHmCommon$msSample)
###output for autoRT
for(msSampleName in unique(fpMqHmCommon$msSample)){
  fwrite(fpMqHmCommon[msSample==msSampleName,.(x=pepQueryFormat,y=Retention)],
         file = paste0('./GS4autoRT/',msSampleName,'training.tsv'), sep = '\t')
  fwrite(fpMqHmCommon[Modified_Peptide==''][msSample==msSampleName,.(x=pepQueryFormat,y=Retention)],
         file = paste0('./GS4autoRT/',msSampleName,'training_noModi.tsv'), sep = '\t')
}

###

#remove these matrices to save memory
rm(fpCombined, mqCombined)
rm(fpMouse, mqMouse)
