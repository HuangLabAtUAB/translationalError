#fpOpen = fread('./fpOpenPointFive/psm.tsv', sep = '\t', select = 1:34, fill=T)
fpOpen = fread('/Volumes/ChenDrive/fpSearchMSNew_pancreatic/fpOpen/psm.tsv', sep = '\t', select = 1:34, fill=T)
colnames(fpOpen) = gsub(colnames(fpOpen), pattern='\\s+', replacement = '_')
dim(fpOpen) # 2602321      34
fpOpen = fpOpen[MSFragger_Localization!=''] #only keep the location is available
dim(fpOpen) # 621,338      34
fpOpen[, msSample:=gsub(Spectrum, pattern='^[0-9]+(_[0-9]{2}_).*', replacement='\\1')]
fpOpen = fpOpen[as.numeric(substring(msSample,2,3)) < 25] #only need the first twenty-four samples

fpOpen[,msmsScanNums := gsub(Spectrum, pattern=".*?\\.([0-9]+)\\..*", replacement = '\\1')]
fpOpen[, msmsScanNum := as.numeric(msmsScanNums)]

dim(fpOpen) #313,516
fpOpen[,spetrNO:=paste0(msSample,msmsScanNums)]

# unique(unlist(strsplit(unique(fpOpen$Assigned_Modifications), ', | ')))
#assigned_modifications are: M(15.9949) N-term(42.0106) C(57.0214)


#below are the results for one sample
# fpOpenEva = fpOpenResEvaFuc(fpOpenRes = fpOpen, gtSeqs = mqMouse$`_01_`$Sequence, 
#                 gtSpectr = mqMouse$`_01_`$spectrumSeq,
#                 minValThres = 0.01) #101   T107  ||2865   111 
# fpOpenResEvaFuc(fpOpenRes = fpOpen, gtSeqs = mqMouse$`_01_`$Sequence, 
#          gtSpectr = mqMouse$`_01_`$spectrumSeq,
#          minValThres = 0.05) #F412   T418  ||2552   424 
# 
# fpOpenResEvaFuc(fpOpenRes = fpOpen, gtSeqs = fpMouse$`_01_`$Peptide, 
#                 gtSpectr = fpMouse$`_01_`$spectrumSeq,
#                 minValThres = 0.05) #365   T465  ||2847   471 

#only evaluate the single sav
fpOpen[,numSAV := stringr::str_count(MSFragger_Localization, pattern="[a-z]")]
table(fpOpen$numSAV)
fpOpen = fpOpen[numSAV==1]
dim(fpOpen) #100,692     38
fpOpenEva = fpOpenResEvaFuc(fpOpenRes = fpOpen, 
                            gtSpectr = mqMouseDt$spectrumSeq,
                            gtSpectr2 = fpMouseDt$spectrumSeq,
                            returnDt = TRUE)

# [1] "-----precision(spectrum)-----"
# 
# FALSE  TRUE 
#  15518  24747 
# [1] "-----sensitivity(spectrum)-----"
# 
# FALSE  TRUE 
# 33909 24853 
table(fpOpenEva$isinAnyGS) #12967 27298 ->12964 27301 
fpOpenEvaCount = fpOpenEva[,.(totalSample=.N,
                              totalTrueSample1=sum(isinGS1),
                              totalTrueSample2=sum(isinGS2),
                              totalTrueAny=sum(isinGS1|isinGS2)),
                           by='msSample']
range(fpOpenEvaRatio$totalTrueAny)#0.5094891 0.7266541
median(fpOpenEvaRatio$totalTrueAny)#0.6921221
fpOpenEvaRatio = fpOpenEvaCount[,.(msSample, totalTrueSample1=totalTrueSample1/totalSample,
                                   totalTrueSample2=totalTrueSample2/totalSample,
                                   totalTrueAny=totalTrueAny/totalSample)]
  
fpOpenEvaRatioMelt = melt.data.table(fpOpenEvaRatio,id.vars = 'msSample',
                                 variable.name = 'type', 
                                 value.name = 'frac')
ggplot(fpOpenEvaRatioMelt, aes(type,frac,col=type))+
  geom_boxplot(outliers = FALSE)+
  geom_jitter(width = 0.1)+ylim(c(0.1,0.8))+ggtitle('precision')

fpMouseDt[,ID_by_fpOpen:=spectrumSeq%in%fpOpenEva[isinGS2==TRUE,spectrSeq]]
mqMouseDt[,ID_by_fpOpen:=spectrumSeq%in%fpOpenEva[isinGS1==TRUE,spectrSeq]]
table(mqMouseDt$ID_by_fpOpen) #identified by fpOpen search
table(fpMouseDt$ID_by_fpOpen) #identified by fpOpen search

#sensitivity evaluation
#build up count matrix
fpOpenEvaSen = list()
fpOpenEvaSen[['GSFromMqCombine']] = mqMouseDt[,.(totalSample=.N,
                                      totalID=sum(ID_by_fpOpen)),
                                    by='msSample']
fpOpenEvaSen[['GSFromFpCombine']] = fpMouseDt[,.(totalSample=.N,
                                       totalID=sum(ID_by_fpOpen)),
                                    by='msSample']
fpOpenEvaSen[['GSFromBoth']] = fpMouseDt[spectrumSeq%in%mqMouseDt$spectrumSeq][,.(totalSample=.N,
                                                 totalID=sum(ID_by_fpOpen)),
                                              by='msSample']
fpOpenEvaSenDt=rbindlist(fpOpenEvaSen,idcol = 'type')
fpOpenEvaSenDt[,frac:=totalID/totalSample]
#plot the sensitivity results
fpOpenEvaSenDt[,type:=factor(type,levels = c('GSFromMqCombine','GSFromFpCombine','GSFromBoth'))]
ggplot(fpOpenEvaSenDt, aes(type,frac,col=type))+
  geom_boxplot(outliers = FALSE)+
  geom_jitter(width = 0.1)+ylim(c(0.2,0.8))+
  ggtitle('sensitivity')

table(fpMouseDt$Modified_Peptide!='',fpMouseDt$ID_by_fpOpen)
range(fpOpenEvaSenDt[type=='GSFromBoth',frac])
# [1] 0.4732544 0.5113895
median(fpOpenEvaSenDt[type=='GSFromBoth',frac])
#0.4835638


#sensitivity with only no modification
fpOpenEvaSen = list()
fpOpenEvaSen[['GSFromMqCombine']] = mqMouseDt[!grepl(Modified_sequence, pattern='\\('),
                                              .(totalSample=.N,
                                                 totalID=sum(ID_by_fpOpen)),
                                              by='msSample']
fpOpenEvaSen[['GSFromFpCombine']] = fpMouseDt[Modified_Peptide=='',.(totalSample=.N,
                                                 totalID=sum(ID_by_fpOpen)),
                                              by='msSample']
fpOpenEvaSen[['GSFromBoth']] = fpMouseDt[spectrumSeq%in%mqMouseDt$spectrumSeq &Modified_Peptide==''][,.(totalSample=.N,
                                                                                  totalID=sum(ID_by_fpOpen)),
                                                                               by='msSample']
fpOpenEvaSenDt2=rbindlist(fpOpenEvaSen,idcol = 'type')
fpOpenEvaSenDt2[,frac:=totalID/totalSample]
ggplot(fpOpenEvaSenDt2, aes(type,frac,col=type))+
  geom_boxplot(outliers = FALSE)+
  geom_jitter(width = 0.1)+ylim(c(0.1,0.5))+ggtitle('without_modi_sensitivity')

##export the fpOpen searched results for autoRT
table(fpOpenEva$Assigned_Modi)
fpOpenEva[, pepSAV4autoRT:=forMatPep(pep = pepSAV, modi = Assigned_Modi)]
fpOpenEva[, pepOrg4autoRT:=forMatPep(pep = PepOrg, modi = Assigned_Modi)]
dim(fpOpenEva)#40265    21
fpOpenEva$pepLength=nchar(fpOpenEva$PepOrg)
max(fpOpenEva$pepLength) #47

for(msSampleName in unique(fpOpenEva$msSample)){
  fwrite(fpOpenEva[msSample==msSampleName,
                   .(x=pepSAV4autoRT,y=Retention)],
         file = paste0('./fpOpenRes4autoRT/',msSampleName,'predict.tsv'), 
         sep = '\t')
  fwrite(fpOpenEva[Assigned_Modi==''][msSample==msSampleName,
                                   .(x=pepSAV4autoRT,y=Retention)],
         file = paste0('./fpOpenRes4autoRT/',msSampleName,'predict_noModi.tsv'), sep = '\t')
  fwrite(fpOpenEva[msSample==msSampleName,
                   .(x=pepOrg4autoRT,y=Retention)],
         file = paste0('./fpOpenRes4autoRT/',msSampleName,'predict_org.tsv'), sep = '\t')
  fwrite(fpOpenEva[Assigned_Modi==''][msSample==msSampleName,
                                   .(x=pepOrg4autoRT,y=Retention)],
         file = paste0('./fpOpenRes4autoRT/',msSampleName,'predict_org_noModi.tsv'), sep = '\t')
}

#for pdeep2
##add charge information
table(fpOpenEva$spetrNO %in% fpOpen$spetrNO)
fpOpenEva[,charge:=fpOpen[fpOpenEva$spetrNO,Charge,on='spetrNO']]
fpOpenEva[,mod4pDeep:=gsub(Assigned_Modi,pattern=',',replacement=';')]
fpOpenEva[,mod4pDeep:=gsub(mod4pDeep,pattern='C\\(57.0214\\)',replacement=',Carbamidomethyl[C]')]
fpOpenEva[,mod4pDeep:=gsub(mod4pDeep,pattern='C\\(57.0214\\)',replacement=',Carbamidomethyl[C]')]
fpOpenEva[,mod4pDeep:=gsub(mod4pDeep,pattern='M\\(15.9949\\)',replacement=',Oxidation[M]')]
fpOpenEva[,mod4pDeep:=gsub(mod4pDeep,pattern='; ',replacement=';')]
fpOpenEva[,mod4pDeep:=gsub(mod4pDeep,pattern='N-term\\(42.0106\\)',replacement='')]
fpOpenEva[,mod4pDeep:=gsub(mod4pDeep,pattern=';$',replacement='')]
fpOpenEva[,mod4pDeep:=gsub(mod4pDeep,pattern='^;',replacement='')]
fpOpenEva[,mod4pDeep:=gsub(mod4pDeep,pattern='\\"',replacement='')]

for(msSampleName in unique(fpOpenEva$msSample)){
  fwrite(fpOpenEva[msSample==msSampleName,
                   .(peptide=gsub(pepSAV,pattern='X',replacement='L'), 
                     modification=mod4pDeep,
                     charge=charge)],
         col.names = F,
         file = paste0('./fpOpenRes4pDeep2/',msSampleName,'predict.tsv'), 
         quote = F, sep = '\t')
}
##add back spect
fpOpenEva[,spectrum:=fpOpen[fpOpenEva$spetrNO,Spectrum,on='spetrNO']]
#to fit for the mgf file
fpOpenEva[,spectrum:=gsub(spectrum, pattern='\\.[0]+',replacement='.')]
fpOpenEva[,pDeepTitle := paste0(gsub(pepSAV,pattern='X',replacement='L'),'|',mod4pDeep,'|',charge)]
for(msSampleName in unique(fpOpenEva$msSample)){
  fwrite(fpOpenEva[msSample==msSampleName,
                   .(index=spectrum,
                     peptide=gsub(pepSAV,pattern='X',replacement='L'), 
                     modification=mod4pDeep,
                     charge=charge,
                     pdeep_title=pDeepTitle)],
         file = paste0('./fpOpenRes4pDeep2/',msSampleName,'format_title.tsv'), 
         quote = F, col.names = T,
         sep = '\t')
}

###autoRT res 
autortRes = list()
for (rtRes in dir('./fpOpenResPred', pattern='tsv')){
  predName = gsub(rtRes, pattern='\\.tsv', replacement='')
  predName = gsub(predName, pattern='([0-9])$', replacement='\\1_')
  autortRes[[predName]] = fread(paste0('./fpOpenResPred/', rtRes))
}
names(autortRes)
predRes = rbindlist(autortRes[unique(fpOpenEva$msSample)],idcol = 'msSampleName')
all(predRes$x == fpOpenEva$pepSAV4autoRT) #TRUE

predResOrg = rbindlist(autortRes[paste0(unique(fpOpenEva$msSample),'org')],
                       idcol = 'msSampleName')
all(predResOrg$x == fpOpenEva$pepOrg4autoRT) #TRUE

fpOpenEva$predictRT = predRes$y_pred * 60
fpOpenEva$predictRT_orgPep = predResOrg$y_pred * 60

fpOpenEva[,overlapWithGS:=sapply(1:nrow(fpOpenEva), 
                                 function(x)sum(isinGS1[x],isinGS2[x]))]

ggplot(fpOpenEva, aes(isinGS1, abs(predictRT-Retention)))+geom_boxplot()
ggplot(fpOpenEva, aes(isinGS2, abs(predictRT-Retention)))+geom_boxplot()
ggplot(fpOpenEva, aes(isinAnyGS, abs(predictRT-Retention)))+geom_boxplot()
ggplot(fpOpenEva, aes(as.character(overlapWithGS), abs(predictRT-Retention)))+
  geom_boxplot()
fpOpenEva[,wilcox.test(abs(predictRT-Retention)~isinAnyGS)]

#######
##precision and sensitivity across different SAVs
fpOpenEva[, sav:=paste0(aaFrom,'->', aaTo)]
fpOpenSavCount = fpOpenEva[,.(totalN=.N,
                              aaFrom=unique(aaFrom),
                              aaTo=unique(aaTo),
                              isinGS1Ratio=sum(isinGS1)/.N,
                             isinGS2Ratio=sum(isinGS2)/.N,
                             isinAnyGSRatio=sum(isinAnyGS)/.N),
                           by=c('sav')]
fpOpenSavCount[,aaFromMass:=aaMass[fpOpenSavCount$aaFrom,Average,on='one_letter_code']]
fpOpenSavCount[,aaToMass:=aaMass[fpOpenSavCount$aaTo,Average,on='one_letter_code']]
fpOpenSavCount[, aaDeltaMass := abs(aaToMass - aaFromMass)]
ggplot(fpOpenSavCount, aes(aaDeltaMass, isinAnyGSRatio))+geom_point()+
  geom_smooth(method = "lm",se = FALSE)
cor.test(fpOpenSavCount$aaDeltaMass, fpOpenSavCount$isinAnyGSRatio) #p=0.2743
cor.test(fpOpenSavCount[totalN>20,aaDeltaMass], 
         fpOpenSavCount[totalN>20,isinAnyGSRatio])#p=0.3092

table(fpOpenSavCount$totalN>20)
table(fpOpenSavCount$totalN>100)

fpOpenSavCountBySample = fpOpenEva[,.(totalN=.N,
                              aaFrom=unique(aaFrom),aaTo=unique(aaTo),
                              isinGS1Ratio=sum(isinGS1)/.N,
                              isinGS2Ratio=sum(isinGS2)/.N,
                              isinAnyGSRatio=sum(isinAnyGS)/.N),
                           by=c('sav','msSample')]
fpOpenSavCountBySample=fpOpenSavCountBySample[totalN>5]
fpSavCount = fpOpenSavCountBySample[,.N,by='sav']# number in each sample
sav2display = fpSavCount[N>10, sav]
length(sav2display)



###
table(fpMouseDt$spectrumSeq %in% mqMouseDt$spectrumSeq)
#sensitivity by using overlapped GS
#mqIntersectFp = mqMouseDt[spectrumSeq%in%fpMouseDt$spectrumSeq]

fpCount4Sen = mqIntersectFp[,.(totalN=.N,
                     idRatio=sum(ID_by_fpOpen)/.N), 
                  by=c('msSample','sav')]
fpSavCount4sen = mqIntersectFp[,.(totalN=.N,
                        idRatio=sum(ID_by_fpOpen)/.N),by='sav']

# ggplot(fpCount4Sen, aes(sav, idRatio))+
#   geom_boxplot()+geom_jitter(width = 0.1)+
#   theme(axis.text.x = element_text(angle = -90))
#fpCount4Sen[, totalNbySAV:=sum(totalN), by='sav']
fpSavCount4sen[,aaFrom:=gsub(sav,pattern='->.*',replacement='')]
fpSavCount4sen[,aaTo:=gsub(sav,pattern='.*->',replacement='')]
fpSavCount4sen[,aaFromMass:=aaMass[fpSavCount4sen$aaFrom,Average,on='one_letter_code']]
fpSavCount4sen[,aaToMass:=aaMass[fpSavCount4sen$aaTo,Average,on='one_letter_code']]
fpSavCount4sen[, aaDeltaMass := abs(aaToMass - aaFromMass)]


sav2display = intersect(sav2display,fpSavCount4sen[totalN>10, sav])
length(sav2display) #60

# ggplot(fpCount4Sen[totalNbySAV>100], aes(sav, idRatio))+
#   geom_boxplot()+geom_jitter(width = 0.1)+
#   theme(axis.text.x = element_text(angle = -90))

fpCount4Sen[sav%in%sav2display,length(unique(sav))]
ggplot(fpCount4Sen[sav%in%sav2display], aes(sav, idRatio))+
  geom_boxplot()+
  #geom_jitter(width = 0.1)+
  theme(axis.text.x = element_text(angle = -90))

#precision boxplot
fpOpenSavCountBySample[sav%in%sav2display,length(unique(sav))]
ggplot(fpOpenSavCountBySample[sav%in%sav2display], 
       aes(sav, isinAnyGSRatio))+geom_boxplot()+
  #geom_jitter(width = 0.1)+
  theme(axis.text.x = element_text(angle = -90))

fpMouseDt[spectrumSeq%in%mqMouseDt$spectrumSeq][,.(totalSample=.N,
                                                   totalID=sum(ID_by_fpOpen)),
                                                by='msSample']

ggplot(fpOpenSavCount[sav%in%sav2display], aes(sav,aaDeltaMass))+
   geom_bar(stat = 'identity')+theme(axis.text.x = element_text(angle = -90))

ggplot(fpOpenSavCount[sav%in%sav2display], aes(aaDeltaMass, isinAnyGSRatio))+
  geom_point()+
  geom_smooth(method = "lm",se = FALSE)+scale_y_continuous(labels = scales::percent)
cor.test(fpOpenSavCount[sav%in%sav2display,aaDeltaMass], 
         fpOpenSavCount[sav%in%sav2display,isinAnyGSRatio])#p0.1305 cor-0.1974257 




cor.test(fpSavCount4sen[sav%in%sav2display,aaDeltaMass], 
         fpSavCount4sen[sav%in%sav2display,idRatio])#p=0.3519
ggplot(fpSavCount4sen[sav%in%sav2display], 
       aes(aaDeltaMass, idRatio))+geom_point()+
  geom_smooth(method = "lm",se = FALSE)+
  scale_y_continuous(labels = scales::percent)
