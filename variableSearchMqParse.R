args <- commandArgs(T)
evidenceDt <- fread(args[1])
#evidenceDt <- fread('/Volumes/ChenDrive/mqXV_evidence.txt')
modiStr <- args[2]
modiSingle <- args[3]
output <- args[4]
colnames(evidenceDt) = gsub(colnames(evidenceDt), pattern='\\s+', replacement='_')
evidenceDt[, msSample:=gsub(Raw_file, pattern='180316(_[0-9]+_).*', replacement='\\1')]
evidenceDt = evidenceDt[grepl(Modifications, 
                              pattern=paste0('(^|,)',modiStr))]

# mqXM[Modifications=='Xle->Met']
# mqXM[Modifications=='Oxidation (M),2 Xle->Met']
##mqXM
evidenceDt = evidenceDt[,.(msSample, msmsScanNum = `MS/MS_scan_number`, 
               msmsScanNums = `MS/MS_scan_numbers`,
               Sequence, 
               Modifications,
               Modified_sequence,
               PEP, Score, Retention_time)]
evidenceDt[['targetPro']] = evidenceDt[[paste0(modiStr,"_Probabilities")]]


evaVarModi = function(resDt, varModi){
  #if(is.null(varModiGS)) varModiGS=varModi 
  #varModiGS can be I/L while varModi can be X
  resAA = gsub(varModi, pattern='.*->',replacement='')
  #example: 'ADRDQYEL(0.001)L(0.006)CL(0.993)DNTR'
  locVal = stringr::str_extract_all(resDt$targetPro, pattern='(0|1)(\\.|)[0-9]*')
  locVal = sapply(locVal, as.numeric)
  locVal = sapply(locVal, function(x)which.max(x)[1]) #which has the max prob 1,2,3
  tempStr = gsub(resDt$targetPro, 
                 pattern='\\((0|1)(\\.|)[0-9]*', replacement='') #"ADRDQYEL)L)CL)DNTR"
  locPos = stringr::str_locate_all(tempStr, pattern = '\\)')
  locPos = lapply(locPos, function(x)x[,1])
  locPos = lapply(locPos, function(x)x-(1:length(x))) #which pos in terms of the aa sequence
  locPos = sapply(1:length(locPos), function(ii)locPos[[ii]][locVal[ii]])
  resDt$startAA = substr(resDt$Sequence, locPos, locPos)
  resDt$pepSAV = paste0(substr(resDt$Sequence,1,(locPos-1)),
                        resAA,
                        substr(resDt$Sequence,locPos+1,nchar(resDt$Sequence)))
  
  resDt$spectrSeq = paste0(resDt$msSample,resDt$msmsScanNum,'_',resDt$pepSAV)
  return(resDt)
}


evidenceDt = evaVarModi(evidenceDt,varModi = modiSingle)
fwrite(evidenceDt, sep = '\t', file = out)

