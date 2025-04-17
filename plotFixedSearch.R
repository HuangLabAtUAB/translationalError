#functions for plotting mq fixed search results, detailing each sample
library(ggplot2)
plotFixedMqPre <- function(mqRes, gtDt1, gtDt2,ylimV=NULL, vec='spectrSeq'){
  #mqXM$spectrSeq %in% mqMouseDt$spectrumSeq
  preDt1 = mqRes[[vec]] %in% gtDt1[["spectrumSeq"]]
  preDt2 = mqRes[[vec]] %in% gtDt2[["spectrumSeq"]]
  preDtU = preDt1|preDt2
  mqRes[['preDt1']] = preDt1 
  mqRes[['preDt2']] = preDt2
  mqRes[['preDtU']] = preDtU #U means union
  mqCount = mqRes[,.(preDt1N=sum(preDt1)/.N,
                     preDt2N=sum(preDt2)/.N,
                     preDtUN=sum(preDtU)/.N),
                  by='msSample']
  mqCount = data.table::melt.data.table(mqCount, 
                                        id.vars = 'msSample',
                                        variable.name = 'group',
                                        value.name = 'ratio')
  print(paste0('median of pre=',mqCount[group=='preDtUN',median(ratio)]))
  gplotRes = ggplot(mqCount,aes(group, ratio))+geom_boxplot(outliers = FALSE)+
    geom_jitter(width = 0.1)
  if(!is.null(ylimV)) gplotRes+ylim(ylimV)
  else gplotRes
}
plotFixedMqSen <- function(mqRes, gtDt1, gtDt2, saavStr,ylimV=NULL, vec='spectrSeq'){
  gt1 = gtDt1[sav%in%saavStr, spectrumSeq]
  gt2 = gtDt2[sav%in%saavStr, spectrumSeq]
  gtU = unique(c(gt1,gt2))
  senDt = data.table(spectrumSeq = gtU,
    isInGt1 = ifelse(gtU%in%gt1, T,F),
    isInGt2 = ifelse(gtU%in%gt2, T,F),
    isIden = ifelse(gtU%in%mqRes[[vec]], T,F))
  senDt$isInBothGt = senDt$isInGt1 & senDt$isInGt2
  senDt$msSample <- gsub(senDt$spectrumSeq,  pattern='^(_[0-9]+_).*',replacement='\\1')
  
  mqCount = senDt[,.(sen4Gt1=sum(isInGt1&isIden)/sum(isInGt1),
                     sen4Gt2=sum(isInGt2&isIden)/sum(isInGt2),
                     sen4Gt=sum(isInBothGt&isIden)/sum(isInBothGt)),
                  by='msSample']
  mqCount = data.table::melt.data.table(mqCount, 
                                        id.vars = 'msSample',
                                        variable.name = 'group',
                                        value.name = 'ratio')
  print(paste0('median of sen=',mqCount[group=='sen4Gt',median(ratio)]))
  gplotRes = ggplot(mqCount,aes(group, ratio))+geom_boxplot(outliers = FALSE)+
    geom_jitter(width = 0.1)+scale_y_continuous(labels = scales::percent)
  if(!is.null(ylimV)) gplotRes+ylim(ylimV)
  else gplotRes
}


