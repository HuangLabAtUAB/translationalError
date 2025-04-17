pFind = fread('./pFind/pFind-Filtered.spectra')
#the "sequence" column is the before the SAV identified by pFind
pFind[,MsSample:=gsub(File_Name, pattern='180316(_[0-9]+_).*', replacement='\\1')]
table(pFind$MsSample)
dim(pFind) #1,590,153    20
pFind[order(-`Mass_Shift(Exp.-Calc.)`)]

pFind=pFind[grepl(pFind$Modification, pattern='->')] #-> means substitution
dim(pFind) #129,579   20

#"3,Glu->Asp[E];"
pFind=pFind[grepl(pFind$Modification, pattern='->[A-Z][a-z][a-z]\\[[A-Z]\\]')] #subset the SAVs
dim(pFind) #97418   20

pFind[,numOfSub := stringr::str_count(Modification, pattern='[A-Z][a-z][a-z]->[A-Z][a-z][a-z]')]
table(pFind$numOfSub) #1: 97412     2: 6 
pFind = pFind[numOfSub==1]
pFind[, pos:= gsub(".*?([0-9]+),[A-Z][a-z][a-z]->[A-Z][a-z][a-z].*", '\\1', Modification)]
pFind[, pos:= gsub(".*?([0-9]+),[^;]+[A-Z][a-z][a-z]->[A-Z][a-z][a-z].*", '\\1', pos)]

# pFind[, pos:= gsub(".*?([0-9]+),[A-Z]*[a-z]+\\[.+\\]\\([A-Z][a-z][a-z]->.*", '\\1', pos)]
# pFind[, pos:= gsub(".*?([0-9]+),[A-Z]*[a-z]+\\+.*", '\\1', pos)]
# pFind[, pos:= gsub(".*?([0-9]+),.*", '\\1', pos)]#21,Delta_H(4)C(2)O(-1)S(1)[S](Ser->Met[S]);
pFind[, pos:=as.numeric(pos)]
pFind[is.na(pos), Modification]

pFind[,startAA:=gsub('.*([A-Z][a-z][a-z])->.*', '\\1', Modification)]
pFind[,endAA:=gsub('.*->([A-Z][a-z][a-z]).*', '\\1', Modification)]

pFind[, startAA:=aaMass[pFind$startAA,one_letter_code,on='three_letter_code']]
pFind[, endAA:=aaMass[pFind$endAA,one_letter_code,on='three_letter_code']]
table(is.na(pFind$startAA))
table(is.na(pFind$endAA)) #F97116   T296so of the endAA is Orn, Dha
pFind = pFind[!is.na(endAA)]
dim(pFind) #97116
table(pFind$endAA) #no I or L, but X
table(pFind$startAA) #no I or L, but X

pFind[startAA!='X',table(substr(Sequence,pos,pos)==startAA)] #All TRUE
pFind[startAA=='X',table(substr(Sequence,pos,pos))] #only I and L
pFind[grepl(Modification, pattern='Xle->'), table(startAA)] #all X
pFind[grepl(Modification, pattern='->Xle'), table(endAA)] #all X

pFind[,pepSAV:=Sequence]

pFind[startAA=='X',startAA:=substr(Sequence,pos,pos)] 

substr(pFind$pepSAV,pFind$pos, pFind$pos) <- pFind$endAA
#now, pepSAV is the sequence after pFind-based SAV

##add this for spectrum-level evaluation
pFind[, spectrumID:=paste0(MsSample,Scan_No)]
pFind[, spectrumSeq:=paste0(spectrumID,'_',pepSAV)]
####


###spectral level evaluation
table(pFind$spectrumSeq %in% mqMouseDt$spectrumSeq) #F57,047 T40,365 |#65487 T31925 
table(mqMouseDt$spectrumSeq %in% pFind$spectrumSeq) #F174,651,  T40,365|#28,857 T32,435 

pFind[, isinGS1:=spectrumSeq%in%mqMouseDt$spectrumSeq]
pFind[, isinGS2:=spectrumSeq%in%fpMouseDt$spectrumSeq]
pFind[, isinAnyGS:= isinGS1|isinGS2]
#pfind precision
table(pFind$isinGS1) #60866 36250 
table(pFind$isinGS2) #60167 36949
table(pFind$isinAnyGS)#57670 39446 
mqMouseDt[,ID_by_pFind:=spectrumSeq%in%pFind$spectrumSeq]
fpMouseDt[, ID_by_pFind:=spectrumSeq%in%pFind$spectrumSeq]
table(mqMouseDt$ID_by_pFind) #22383 T36379 
table(fpMouseDt$ID_by_pFind) #31548 T37073 

pFindRatio = pFind[,.(totalTrueSample1=sum(isinGS1)/.N,
                      totalTrueSample2=sum(isinGS1)/.N,
                      totalTrueAny=sum(isinAnyGS)/.N), by='MsSample']
range(pFindRatio$totalTrueAny)
median(pFindRatio$totalTrueAny)
pFindRatioMelt = melt.data.table(pFindRatio,id.vars = 'MsSample',
                                     variable.name = 'type', 
                                     value.name = 'frac')
ggplot(pFindRatioMelt, aes(type,frac,col=type))+
  geom_boxplot(outliers = FALSE)+
  geom_jitter(width = 0.1)+ylim(c(0,0.7))+ggtitle('precision')

#sens
pfindEvaSen = list()
pfindEvaSen[['GSFromMqCombine']] = mqMouseDt[,.(totalSample=.N,
                                                 totalID=sum(ID_by_pFind)),
                                              by='msSample']
pfindEvaSen[['GSFromFpCombine']] = fpMouseDt[,.(totalSample=.N,
                                                 totalID=sum(ID_by_pFind)),
                                              by='msSample']
pfindEvaSen[['GSFromBoth']] = fpMouseDt[spectrumSeq%in%mqMouseDt$spectrumSeq][,.(totalSample=.N,
                                                                                  totalID=sum(ID_by_pFind)),
                                                                               by='msSample']
pfindEvaSenDt=rbindlist(pfindEvaSen,idcol = 'type')
pfindEvaSenDt[,frac:=totalID/totalSample]
#plot the sensitivity results
pfindEvaSenDt[,type:=factor(type,levels = c('GSFromMqCombine','GSFromFpCombine','GSFromBoth'))]
ggplot(pfindEvaSenDt, aes(type,frac,col=type))+
  geom_boxplot(outliers = FALSE)+
  geom_jitter(width = 0.)+ylim(c(0,0.7))+
  ggtitle('sensitivity')
##
range(pfindEvaSenDt[type=='GSFromBoth',frac])
# [1] 0.5994513 0.6843783
median(pfindEvaSenDt[type=='GSFromBoth',frac])
#0.6623611

##precision and sensitivity across different SAVs
pFind[, sav:=paste0(startAA,'->', endAA)]
pFindSavCount = pFind[,.(totalN=.N,
                              aaFrom=unique(startAA),aaTo=unique(endAA),
                              isinGS1Ratio=sum(isinGS1)/.N,
                              isinGS2Ratio=sum(isinGS2)/.N,
                              isinAnyGSRatio=sum(isinAnyGS)/.N),
                           by=c('sav')]
pFindSavCount[,aaFromMass:=aaMass[pFindSavCount$aaFrom,Average,on='one_letter_code']]
pFindSavCount[,aaToMass:=aaMass[pFindSavCount$aaTo,Average,on='one_letter_code']]
pFindSavCount[, aaDeltaMass := abs(aaToMass - aaFromMass)]




ggplot(pFindSavCount, aes(aaDeltaMass, isinAnyGSRatio))+geom_point()
cor.test(pFindSavCount$aaDeltaMass, pFindSavCount$isinAnyGSRatio) #p=0.03015
#cor.test(fpOpenSavCount[totalN>20,aaDeltaMass], fpOpenSavCount[totalN>20,isinAnyGSRatio])#p=0.5874

table(fpOpenSavCount$totalN>20)# 107   120 
table(fpOpenSavCount$totalN>100) #157    70 

##comparison of precision between pFind and fpOpen
commonSAVs = intersect(pFindSavCount$sav, 
                       fpOpenSavCount$sav)
cor.test(pFindSavCount[commonSAVs,isinAnyGSRatio,on='sav'], 
         fpOpenSavCount[commonSAVs,isinAnyGSRatio,on='sav'])#pearson0.637 #sp 0.677
cor.test(pFindSavCount[commonSAVs,isinAnyGSRatio,on='sav'], 
         fpOpenSavCount[commonSAVs,isinAnyGSRatio,on='sav'])
precisionDoubleDt = data.table(pFindRatio=pFindSavCount[commonSAVs,isinAnyGSRatio,on='sav'],
                               fpRatio=fpOpenSavCount[commonSAVs,isinAnyGSRatio,on='sav'],
                               sav=commonSAVs)
precisionDoubleDt[,label:=ifelse(pFindRatio>0.65&fpRatio>0.65,sav,'')]
#precisionDoubleDt[,label:=ifelse(pFindRatio<0.1&fpRatio<0.1,sav,label)]

table(sav2display %in% sav2display_pFind)
table(sav2display_pFind %in% sav2display)

ggplot(precisionDoubleDt, aes(pFindRatio,fpRatio))+geom_point()+
  ggrepel::geom_text_repel(mapping = aes(label=label),max.overlaps=100)
ggplot(precisionDoubleDt[sav%in%sav2display], aes(pFindRatio,fpRatio))+geom_point()+
  ggrepel::geom_label_repel(mapping = aes(label=label),max.overlaps=200, alpha=.7)+
  geom_smooth(method = "lm", se = FALSE)+
  scale_y_continuous(labels = scales::percent)+scale_x_continuous(labels = scales::percent)
cor.test(pFindSavCount[sav2display,isinAnyGSRatio,on='sav'], 
         fpOpenSavCount[sav2display,isinAnyGSRatio,on='sav'])
##comparision in another way
pFindSavCount_byS = pFind[,.(totalN=.N,
                         isinGS1Ratio=sum(isinGS1)/.N,
                         isinGS2Ratio=sum(isinGS2)/.N,
                         isinAnyGSRatio=sum(isinAnyGS)/.N),
                      by=c('MsSample')]
fpOpenSavCount_byS = fpOpenEva[,.(totalN=.N,
                             isinGS1Ratio=sum(isinGS1)/.N,
                             isinGS2Ratio=sum(isinGS2)/.N,
                             isinAnyGSRatio=sum(isinAnyGS)/.N),
                          by=c('msSample')]
precisionDoubleDt_byS = data.table(pfind=pFindSavCount_byS$isinAnyGSRatio,
                                   fp=fpOpenSavCount_byS[pFindSavCount_byS$MsSample,isinAnyGSRatio,on='msSample'],
                                   msSample=pFindSavCount_byS$MsSample)
cor.test(precisionDoubleDt_byS$pfind, precisionDoubleDt_byS$fp) #0.8602556
ggplot(precisionDoubleDt_byS, aes(fp,pfind))+geom_point()
##

pFindSavCountBySample = pFind[,.(totalN=.N,
                                      aaFrom=unique(startAA),aaTo=unique(endAA),
                                      isinGS1Ratio=sum(isinGS1)/.N,
                                      isinGS2Ratio=sum(isinGS2)/.N,
                                      isinAnyGSRatio=sum(isinAnyGS)/.N),
                                   by=c('sav','MsSample')]
pFindSavCountBySample=pFindSavCountBySample[totalN>5]

tmpSampleCount = pFindSavCountBySample[,.N,by='sav']
sav2display = tmpSampleCount[N>3, sav]
#precision boxplot
ggplot(pFindSavCountBySample[sav%in%sav2display], 
       aes(sav, isinAnyGSRatio))+geom_boxplot(outliers = F)+geom_jitter(width = 0.1)+
  theme(axis.text.x = element_text(angle = -90))

ggplot(fpOpenSavCount[sav2display,,on='sav'], aes(sav,aaDeltaMass))+
  geom_bar(stat = 'identity')+theme(axis.text.x = element_text(angle = -90))


# fpMouseDt[spectrumSeq%in%mqMouseDt$spectrumSeq][,.(totalSample=.N,
#                                                    totalID=sum(ID_by_fpOpen)),
#                                                 by='msSample']


###sensitivity
#sensitivity by using overlapped GS
tmp = fpMouseDt[spectrumSeq %in% mqMouseDt$spectrumSeq]
pFindCount4Sen = tmp[,.(totalN=.N,
                     idRatio=sum(ID_by_pFind)/.N), 
                  by=c('msSample','sav')]
ggplot(pFindCount4Sen, aes(sav, idRatio))+
  geom_boxplot()+geom_jitter(width = 0.1)+
  theme(axis.text.x = element_text(angle = -90))

pFindCount4Sen[,aveIDratio := mean(idRatio), by='sav']
pFindCount4Sen[, totalNbySAV:=sum(totalN), by='sav']
pFindCount4Sen = pFindCount4Sen[order(aveIDratio)]
View(fpCount4Sen[!duplicated(sav),.(sav,aveIDratio,totalNbySAV)])

ggplot(pFindCount4Sen[totalNbySAV>50], aes(sav, idRatio))+
  geom_boxplot()+geom_jitter(width = 0.1)+
  theme(axis.text.x = element_text(angle = -90))

#sensitivity comparison between fpOpen and pFind
#correlation of sensitivity between pFind and fpOpen
tmp = fpMouseDt[spectrumSeq %in% mqMouseDt$spectrumSeq]
senDoubleDt = tmp[,.(by_fpOpen=sum(ID_by_fpOpen)/.N, 
                     by_pFind=sum(ID_by_pFind)/.N,
                     totalN=.N),
                  by=c("sav")]
table(senDoubleDt$totalN>100)
ggplot(senDoubleDt, aes(by_fpOpen,by_pFind))+geom_point()
ggplot(senDoubleDt[sav%in%sav2display], aes(by_pFind,by_fpOpen))+
  geom_point()+
  geom_smooth(method = "lm", se = FALSE)+
  scale_y_continuous(labels = scales::percent)+
  scale_x_continuous(labels = scales::percent)
cor.test(senDoubleDt[sav%in%sav2display,by_fpOpen], 
         senDoubleDt[sav%in%sav2display, by_pFind])

senDoubleDt_byS = tmp[,.(totalN=.N,
                         ID_by_fpOpen=sum(ID_by_fpOpen)/.N,
                         ID_by_pFind=sum(ID_by_pFind)/.N),
                                                 by=c('msSample')]
cor.test(senDoubleDt_byS$ID_by_fpOpen,senDoubleDt_byS$ID_by_pFind)
ggplot(senDoubleDt_byS, aes(ID_by_fpOpen,ID_by_pFind))+geom_point()
