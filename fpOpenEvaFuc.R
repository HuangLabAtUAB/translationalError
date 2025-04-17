##Need to do: the s/p rate for the modified and non modified

fpOpenResEvaFuc <- function(fpOpenRes,
                            minValThres=0.04,
                            #gtSeqs, 
                            gtSpectr,
                            #gtSeqs2=NULL,
                            gtSpectr2=NULL, #can use two gs dataset
                            returnDt =F){
  
  #the abs deviations between open search to all the sav pairs, what is the minum?
  #update: the start AA need to be fixed first
  #zero means there is a sav record for that
  startAA = toupper(gsub(fpOpenRes$MSFragger_Localization, 
                         pattern='.*([a-z]).*', replacement='\\1'))
  filter = startAA%in%unique(aaDelta$aaFrom)
  startAA = startAA[filter]
  fpOpenRes = fpOpenRes[filter]
  
  minValue = rep(NA, nrow(fpOpenRes))
  whichMin = rep(NA, nrow(fpOpenRes))
  endAA = rep(NA, nrow(fpOpenRes))
  deltaTheory = rep(NA, nrow(fpOpenRes))
  for (idx in 1:nrow(fpOpenRes)){
    aaDeltaSub = aaDelta[aaFrom==startAA[idx]]
    minValue[idx] = min(abs(fpOpenRes$Delta_Mass[idx]-aaDeltaSub$detaMass))
    whichMin[idx] = which.min(abs(fpOpenRes$Delta_Mass[idx]-aaDeltaSub$detaMass))
    endAA[idx] = aaDeltaSub[whichMin[idx],aaTo]
    deltaTheory[idx] = aaDeltaSub[whichMin[idx],detaMass]
    
  } 
  
  fpEvaDt = data.table(spetrNO = fpOpenRes$spetrNO, 
                   msSample=fpOpenRes$msSample,
                   msmsScanNums = fpOpenRes$msmsScanNum,
                   Retention = fpOpenRes$Retention,
                   Assigned_Modi = fpOpenRes$Assigned_Modifications,
                   whichMin = whichMin,
                   aaFrom = startAA,
                   aaTo = endAA,
                   deltaTheory = deltaTheory,
                   minValue = minValue,
                   deltaMass = fpOpenRes$Delta_Mass,
                   PepOrg = fpOpenRes$Peptide,
                   PepLoc = fpOpenRes$MSFragger_Localization)#the modified AA is lower case
                   
  #pepSAV: the peptide sequence after deltaMass-inferred SAV
  fpEvaDt[, pepSAV:=sapply(1:nrow(fpEvaDt), function(ii)gsub(PepLoc[ii], pattern='[a-z]', 
                                                       replacement = aaTo[ii]))]
  fpEvaDt[,spectrSeq:=paste0(msSample,msmsScanNums,'_',pepSAV)]
  
  #the users themselves can do filtering by whether delta mass fit to the theoretical and 
  #and whether the minimum mass change corresponds to the starting AA
  fpEvaDt = fpEvaDt[minValue< minValThres]
  
  print("-----precision(spectrum)-----")
  print(table(fpEvaDt$spectrSeq %in% gtSpectr))
  print("-----sensitivity(spectrum)-----")
  print(table(gtSpectr %in% fpEvaDt$spectrSeq)) 

  if(returnDt == TRUE) {
    fpEvaDt[, isinGS1 := spectrSeq%in%gtSpectr]
    if(!is.null(gtSpectr2)){
      fpEvaDt[, isinGS2 := spectrSeq%in%gtSpectr2]
      fpEvaDt[, isinAnyGS := spectrSeq%in%c(gtSpectr,gtSpectr2)]
    }
    return(fpEvaDt)
  }
  
  
}