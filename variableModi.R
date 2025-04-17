library(data.table)

#obtain the most abundant SAAVs for fixed search
savTotalCount = mqMouseDt[,.N,by=c('msSample','sav')]
savTotalCountBySample = savTotalCount[,.(totalCount=sum(N),cancer='pancreatic'),by='msSample']
saveRDS(savTotalCountBySample,file = 'savTotalCountBySamplePanc.rds')
savTotalCount[,.(medianCount = median(N)),by='sav']
savTotalCountMedian = savTotalCount[,.(medianCount = median(N)),by='sav']
savTotalCountMedian = savTotalCountMedian[order(-medianCount)]

savTotalCountMedian[1:10]
#"V->X",""D->E"",I->V"  

savTotalCountSub=savTotalCount[sav%in%savTotalCountMedian[1:20,sav]]

savTotalCountSub[, fromAA := gsub(sav, pattern='->.*', replacement='')]
savTotalCountSub[, toAA := gsub(sav, pattern='.*->', replacement='')]
savTotalCountSub[, fromAAtype:=aaMass[savTotalCountSub$fromAA,AAtype,on='one_letter_code']]
savTotalCountSub[, toAAtype:=aaMass[savTotalCountSub$toAA,AAtype,on='one_letter_code']]
savTotalCountSub[, sametypeConvert:=ifelse(fromAAtype==toAAtype,'sameType','diffType')]

savTotalCountSub[,sav:=factor(sav,levels = savTotalCountMedian[1:20,sav])]
ggplot(savTotalCountSub,aes(sav,N,col=sametypeConvert))+geom_boxplot(outliers = F)+
  geom_jitter(width = 0.1)




#D->E
mqVXED_ED = fread('/Volumes/ChenDrive/mqSearchFixed/mqVXED_ED.tsv')
table(mqVXED_ED$spectrSeq %in% mqMouseDt$spectrumSeq) #1845   T1308  #precision
table(mqVXED_ED$spectrSeq %in% fpMouseDt$spectrumSeq) # 1913   T1240
table(fpMouseDt[sav=='D->E',spectrumSeq] %in% mqVXED_ED$spectrSeq)
# FALSE  TRUE 
# 1297   1232 
table(mqMouseDt[sav=='D->E',spectrumSeq] %in% mqVXED_ED$spectrSeq)
# FALSE  TRUE 
# 988   1300 
plotFixedMqPre(mqRes = mqVXED_ED,gtDt1 = mqMouseDt, 
               gtDt2 = fpMouseDt,ylimV = c(0,1))
plotFixedMqSen(mqRes = mqVXED_ED,gtDt1 = mqMouseDt, 
               gtDt2 = fpMouseDt,saavStr = 'D->E',ylimV = c(0,1))
#add fpOpen for comparison
plotFixedMqPre(mqRes = fpOpenEva[sav=='E->D'],gtDt1 = mqMouseDt, gtDt2 = fpMouseDt)
plotFixedMqSen(mqRes = fpOpenEva[sav=='E->D'],gtDt1 = mqMouseDt, 
               gtDt2 = fpMouseDt,saavStr = 'E->D')

#V->X
mqVXED_VX = fread('/Volumes/ChenDrive/mqSearchFixed/mqVXED_VX.tsv')
table(mqVXED_VX$spectrSeq %in% mqMouseDt$spectrumSeq) #1799  2318   #precision
table(mqVXED_VX$spectrSeq %in% fpMouseDt$spectrumSeq) # 1892  2225 
table(fpMouseDt[sav=='V->X',spectrumSeq] %in% mqVXED_VX$spectrSeq)
# FALSE  TRUE 
# 2719  2225 
table(mqMouseDt[sav=='V->X',spectrumSeq] %in% mqVXED_VX$spectrSeq)
# FALSE  TRUE 
# 1949  2317 
###V->X
mqVX_VX = fread('/Volumes/ChenDrive/mqSearchFixed/mqVX_VX.tsv')
table(mqVX_VX$spectrSeq %in% mqMouseDt$spectrumSeq) #2835  2390  #precision
table(mqVX_VX$spectrSeq %in% fpMouseDt$spectrumSeq) # 2938  2287 
table(fpMouseDt[sav=='V->X',spectrumSeq] %in% mqVX_VX$spectrSeq)
# FALSE  TRUE 
# 2657  2287 
table(mqMouseDt[sav=='V->X',spectrumSeq] %in% mqVX_VX$spectrSeq)
# FALSE  TRUE 
# 1877  2389 
plotFixedMqPre(mqRes = mqVX_VX,gtDt1 = mqMouseDt, 
               gtDt2 = fpMouseDt,ylimV = c(0,1))
plotFixedMqSen(mqRes = mqVX_VX,gtDt1 = mqMouseDt, 
               gtDt2 = fpMouseDt,saavStr = 'V->X',ylimV = c(0,1))
plotFixedMqPre(fpOpenEva[sav=='V->X'],gtDt1 = mqMouseDt, gtDt2 = fpMouseDt,
                    ylimV = c(0,1))
plotFixedMqSen(fpOpenEva[sav=='V->X'],gtDt1 = mqMouseDt, gtDt2 = fpMouseDt,
                ylimV = c(0,1),saavStr = 'V->X')

#mqIV from mqXV
mqXV = fread('/Volumes/ChenDrive/mqSearchFixed/mqXV_XV.tsv')
dim(mqXV)#6167   13
table(mqXV$startAA) #I    L 
#3707 2460 
plotFixedMqPre(mqXV[startAA=='I'],gtDt1 = mqMouseDt, gtDt2 = fpMouseDt,
               ylimV = c(0,1))
plotFixedMqSen(mqXV[startAA=='I'],gtDt1 = mqMouseDt, gtDt2 = fpMouseDt,
               ylimV = c(0,1),saavStr = 'I->V')

##fpED
fpED[,sav:=ifelse(grepl(Assigned_Modifications,pattern='D'),'D->E','E->D')]
table(grepl(fpED$Assigned_Modifications,pattern='D'))
# D->E1 1793 E->D 17870 

plotFixedMqPre(fpED[sav=='D->E'],gtDt1 = mqMouseDt, gtDt2 = fpMouseDt,
               ylimV = c(0,1))
plotFixedMqSen(fpED[sav=='D->E'],gtDt1 = mqMouseDt, gtDt2 = fpMouseDt,
               ylimV = c(0,1),saavStr = 'D->E')
#
plotFixedMqPre(fpOpenEva[sav=='D->E'],gtDt1 = mqMouseDt, gtDt2 = fpMouseDt,
               ylimV = c(0,1))#0.74537037037037"
plotFixedMqSen(fpOpenEva[sav=='D->E'],gtDt1 = mqMouseDt, gtDt2 = fpMouseDt,
               ylimV = c(0,1),saavStr = 'D->E')

##fp fixed search identified many 'false positives'
table(fpED[sav=='D->E',spectrumID] %in% fpMouseDt$spectrumID)
# FALSE  TRUE 
# 10841  7029 

#fpVX
fpVX[,startAA:=gsub(Assigned_Modifications, pattern='.*(I|L|V)\\(.*',replacement='\\1')]
table(fpVX$startAA)
fpVX[,sav:=ifelse(startAA=='V','V->X',paste0(startAA,'->V'))]
table(fpVX$sav)
# V->X  X->V 
# 19514 17657 
#I->V  L->V  V->X 
#9915  7742 19514 
plotFixedMqPre(fpVX[sav=='V->X'],gtDt1 = mqMouseDt, gtDt2 = fpMouseDt,
               ylimV = c(0,1))
plotFixedMqSen(fpVX[sav=='V->X'],gtDt1 = mqMouseDt, gtDt2 = fpMouseDt,
               ylimV = c(0,1), saavStr = 'V->X')



#I->V
plotFixedMqPre(fpVX[sav=='I->V'],gtDt1 = mqMouseDt, gtDt2 = fpMouseDt,
               ylimV = c(0,1))#0.381524914786561
plotFixedMqSen(fpVX[sav=='I->V'],gtDt1 = mqMouseDt, gtDt2 = fpMouseDt,
               ylimV = c(0,1),saavStr = 'I->V')
plotFixedMqPre(fpOpenEva[sav=='I->V'],gtDt1 = mqMouseDt, gtDt2 = fpMouseDt,
               ylimV = c(0,1))
plotFixedMqSen(fpOpenEva[sav=='I->V'],gtDt1 = mqMouseDt, gtDt2 = fpMouseDt,
               ylimV = c(0,1),saavStr = 'I->V')

#why the fixed search performed so poorly in precision?
commonTmp = intersect(fpED[sav=='D->E',spectrumID],fpMouseDt$spectrumID)
debugTmp[spectrumID=='_01_1672']
debugTmp = fpMouseDt[spectrumID%in%commonTmp]
debugSaavTypeDt =data.table(saav=names(table(debugTmp$sav)), 
                            count=as.numeric(table(debugTmp$sav))) 
debugSaavTypeDt[, startAA := gsub(saav,pattern='->.*',replacement='')]
debugSaavTypeDt[, endAA := gsub(saav,pattern='.*->',replacement='')]
debugSaavTypeDt[, deltaMass:=aaMass[endAA,Monoisotopic,on='one_letter_code']-aaMass[startAA,Monoisotopic,on='one_letter_code']]
debugSaavTypeDt=debugSaavTypeDt[order(-count)]
debugSaavTypeDt = debugSaavTypeDt[1:10]
debugSaavTypeDt$saav=factor(debugSaavTypeDt$saav,levels = debugSaavTypeDt$saav)
debugSaavTypeDt[, similarity:=abs(deltaMass-14.01565)]
#top ones are V->X I->V 
ggplot(debugSaavTypeDt[1:10], aes(saav, count))+geom_bar(stat = 'identity')
ggplot(debugSaavTypeDt[1:10], aes(saav, similarity))+geom_bar(stat = 'identity')
#conclusion - the mass difference ED(14.01565) VX(14.01569), 
#are very close or identical, and the location of fp fixed is still poor

##same debug for mqDE
mqVXED_ED$spectrumID = paste0(mqVXED_ED$msSample, mqVXED_ED$msmsScanNum)
commonTmp2 = intersect(mqVXED_ED$spectrumID,fpMouseDt$spectrumID)
debugTmp2 = fpMouseDt[spectrumID%in%commonTmp2]
debugSaavTypeDt2 =data.table(saav=names(table(debugTmp2$sav)), 
                            count=as.numeric(table(debugTmp2$sav))) 
debugSaavTypeDt2[, startAA := gsub(saav,pattern='->.*',replacement='')]
debugSaavTypeDt2[, endAA := gsub(saav,pattern='.*->',replacement='')]
debugSaavTypeDt2[, deltaMass:=aaMass[endAA,Monoisotopic,on='one_letter_code']-aaMass[startAA,Monoisotopic,on='one_letter_code']]
debugSaavTypeDt2=debugSaavTypeDt2[order(-count)]
debugSaavTypeDt2 = debugSaavTypeDt2[1:10]
debugSaavTypeDt2$saav=factor(debugSaavTypeDt2$saav,levels = debugSaavTypeDt2$saav)
debugSaavTypeDt2[, similarity:=abs(deltaMass-14.01565)]

ggplot(debugSaavTypeDt2[1:10], aes(saav, count))+geom_bar(stat = 'identity')
ggplot(debugSaavTypeDt2[1:10], aes(saav, similarity))+geom_bar(stat = 'identity')




