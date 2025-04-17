library(data.table)

dp = fread('./Dp/allPeptides.txt')
colnames(dp) = gsub(colnames(dp), pattern='\\s+', replacement='_')
dp = dp[DP_base_sequence!='']
##
dp[, msSample:=gsub(Raw_file, pattern='180316(_[0-9]+_).*', replacement='\\1')]


dpEva = dpEvaFuc(resDt = dp, gtSpectr = mqMouseDt[humanDetected=='Yes', spectrumSeq], 
               gtSpectr2 = fpMouseDt[humanDetected=='Yes', spectrumSeq],
               minValThres = 0.04, 
               minProb = 0.5, 
               returnDt = TRUE) #False746 TRUE 453 || F924 T453


dpEvaCount = dpEva[,.(totalSample=.N,
                              totalTrueSample1=sum(isinGS1),
                              totalTrueSample2=sum(isinGS2),
                              totalTrueAny=sum(isinGS1|isinGS2)),
                           by='msSample']
dpEvaRatio = dpEvaCount[,.(msSample, totalTrueSample1=totalTrueSample1/totalSample,
                                   totalTrueSample2=totalTrueSample2/totalSample,
                                   totalTrueAny=totalTrueAny/totalSample)]
range(dpEvaRatio$totalTrueAny) #0.2222222 0.4513834
median(dpEvaRatio$totalTrueAny) #0.4293548
dpEvaRatioMelt = melt.data.table(dpEvaRatio,id.vars = 'msSample',
                                     variable.name = 'type', 
                                     value.name = 'frac')
#precison
ggplot(dpEvaRatioMelt, aes(type,frac,col=type))+
  geom_boxplot(outliers = FALSE)+
  geom_jitter(width = 0.1)+ylim(c(0.1,0.8))

fpMouseDt[,ID_by_dp:=spectrumSeq%in%dpEva[isinGS2==TRUE,spectrSeq]]
mqMouseDt[,ID_by_dp:=spectrumSeq%in%dpEva[isinGS1==TRUE,spectrSeq]]
mqMouseDt[humanDetected=='Yes', table(ID_by_dp)] #16872  10268 
fpMouseDt[humanDetected=='Yes', table(ID_by_dp)] #21735 10119 

#sensitivity evaluation
#build up count matrix
dpEvaSen = list()
dpEvaSen[['GSFromMqCombine']] = mqMouseDt[humanDetected=='Yes',.(totalSample=.N,
                                                 totalID=sum(ID_by_dp)),
                                              by='msSample']
dpEvaSen[['GSFromFpCombine']] = fpMouseDt[humanDetected=='Yes',.(totalSample=.N,
                                                 totalID=sum(ID_by_dp)),
                                              by='msSample']
dpEvaSen[['GSFromBoth']] = fpMouseDt[spectrumSeq%in%mqMouseDt$spectrumSeq][humanDetected=='Yes',.(totalSample=.N,
                                                                                  totalID=sum(ID_by_dp)),
                                                                               by='msSample']
dpEvaSenDt=rbindlist(dpEvaSen,idcol = 'type')
dpEvaSenDt[,frac:=totalID/totalSample]

range(dpEvaSenDt[type=='GSFromBoth',frac])
# [1] 0.4732544 0.5113895
median(dpEvaSenDt[type=='GSFromBoth',frac])
#0.4835638
#plot the sensitivity results
dpEvaSenDt[,type:=factor(type,levels = c('GSFromMqCombine','GSFromFpCombine','GSFromBoth'))]

ggplot(dpEvaSenDt, aes(type,frac,col=type))+
  geom_boxplot(outliers = FALSE)+
  geom_jitter(width = 0.1)+ylim(c(0,0.7))
