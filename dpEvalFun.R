dpEvaFuc = function(resDt,
                    gtSpectr, 
                    gtSpectr2=NULL, #can use two gs dataset
                    minValThres=0.04,
                    minProb=0, returnDt =F){
  
  #example: 'E(0.487)P(0.487)P(0.024)HLMSIFK'
  locVal = stringr::str_extract_all(resDt$DP_probabilities, pattern='(0|1)(\\.|)[0-9]*')
  locVal = sapply(locVal, as.numeric)
  locVal = sapply(locVal, function(x)which.max(x)[1]) #which has the max prob 1,2,3
  tempStr = gsub(resDt$DP_probabilities, pattern='\\((0|1)(\\.|)[0-9]*', replacement='') #"VEM)P)SLK"
  locPos = stringr::str_locate_all(tempStr, pattern = '\\)')
  locPos = lapply(locPos, function(x)x[,1])
  locPos = lapply(locPos, function(x)x-(1:length(x))) #which pos in terms of the aa sequence
  locPos = sapply(1:length(locPos), function(ii)locPos[[ii]][locVal[ii]])
  startAA = substr(resDt$DP_base_sequence, locPos, locPos)
  
  filter = startAA%in%unique(aaDelta$aaFrom)
  startAA = startAA[filter]
  resDt = resDt[filter]
  
  minValue = rep(NA, nrow(resDt))
  whichMin = rep(NA, nrow(resDt))
  endAA = rep(NA, nrow(resDt))
  deltaTheory = rep(NA, nrow(resDt))
  for (idx in 1:nrow(resDt)){
    aaDeltaSub = aaDelta[aaFrom==startAA[idx]]
    minValue[idx] = min(abs(resDt$DP_mass_difference[idx]-aaDeltaSub$detaMass))
    whichMin[idx] = which.min(abs(resDt$DP_mass_difference[idx]-aaDeltaSub$detaMass))
    endAA[idx] = aaDeltaSub[whichMin[idx],aaTo]
    deltaTheory[idx] = aaDeltaSub[whichMin[idx],detaMass]
  } 
  
  dpEva = data.table(index = 1:nrow(resDt), 
                     msSample= resDt$msSample,
                     whichMin = whichMin,
                     aaFrom = startAA,
                     aaTo = endAA,
                     deltaTheory = deltaTheory,
                     minValue = minValue,
                     
                     deltaMass = resDt$DP_mass_difference,
                     posProb = resDt$DP_positional_probability,
                     PepOrg = resDt$DP_base_sequence,
                     scanNum = resDt$MSMS_scan_numbers,
                     Retention_time = resDt$Retention_time)
  
  dpEva[, pepSAV:=PepOrg]
  substr(dpEva$pepSAV, locPos, locPos) = dpEva$aaTo
  dpEva[, spectrSeq := paste0(msSample,scanNum,'_',pepSAV)]
  
  ##filtered by common knowledge and data sense
  dpEvaThres = dpEva[minValue < minValThres & posProb>minProb]
  
  # print("-----precision(spectrum)-----")
  # print(table(dpEvaThres$spectrSeq %in% gtSpectr))
  # print("-----sensitivity(spectrum)-----")
  # print(table(gtSpectr %in% dpEvaThres$spectrSeq)) 
  
  if(returnDt == TRUE) {
    dpEvaThres[, isinGS1 := spectrSeq%in%gtSpectr]
    if(!is.null(gtSpectr2)){
      dpEvaThres[, isinGS2 := spectrSeq%in%gtSpectr2]
      dpEvaThres[, isinAnyGS := spectrSeq%in%c(gtSpectr,gtSpectr2)]
    }
    return(dpEvaThres)
  }
}


