library(data.table)
hmPairs = fread('~/translationError/allTrypsinPeptidePairsWithOneSAV_2maxMissCleavge.csv')

#sav in all possible human and mouse pairs
humanPepAA = strsplit(hmPairs$humanPep, '') #str -> vec
mousePepAA = strsplit(hmPairs$mousePep, '') #str -> vec
diffIdx = sapply(1:length(humanPepAA), function(x)which(humanPepAA[[x]] != mousePepAA[[x]]))
hmPairs[, savIdx := diffIdx]
hmPairs[, humanAA := sapply(1:length(diffIdx), function(x)humanPepAA[[x]][diffIdx[x]])]
hmPairs[, mouseAA := sapply(1:length(diffIdx), function(x)mousePepAA[[x]][diffIdx[x]])]
hmPairs[, sav:=paste0(humanAA, '->', mouseAA)]
sort(table(hmPairs$sav), decreasing = T)
rm(humanPepAA, humanPepAA)

dim(hmPairs) #294,513      9
table(hmPairs$humanPep %in% hmPairs$mousePep)
hmCommonPep = intersect(hmPairs$humanPep,hmPairs$mousePep)
length(hmCommonPep) #1419
# hmPairs[humanPep=='VLMLLYSSKK' | mousePep=='VLMLLYSSKK']
# hmPairs[humanPep=='MLMRREQVLK' | mousePep=='MLMRREQVLK']
hmPairs = hmPairs[!(humanPep%in%hmCommonPep) & !(mousePep%in%hmCommonPep)]
dim(hmPairs) #291,925 
table(hmPairs$sav%in%c('I->L','L->I'))#285655   6270 #these are not discoverable
hmPairs = hmPairs[!hmPairs$sav%in%c('I->L','L->I')]
table(hmPairs$mouseAA%in%c('I','L'))#253,486  32,169
hmPairs[, endAA:=ifelse(mouseAA%in%c('I','L'), 'X',mouseAA )]
hmPairs[, sav:=paste0(humanAA, '->', endAA)]
##for endAA of I or L, credit will be given as long as tools can identify X

##for this evaluation, we don't need to consider the genes, therefore removing duplicated rows
# QEVPSWLENMAYEHHYK QEVPSWLENMAFEHHYK     DDX3X    D1Pas1     12       Y       F   Y->F
# QEVPSWLENMAYEHHYK QEVPSWLENMAFEHHYK     DDX3X     Ddx3x     12       Y       F   Y->F
hmPairsShort = hmPairs[,.(dbKey=dbKey[1],humanGene=humanGene[1],mouseGene=mouseGene[1]),
                          by=c('humanPep','mousePep','savIdx','humanAA','mouseAA','endAA',
                               'sav','humanAAtype','mouseAAtype','convertType','sametypeConvert')]
dim(hmPairsShort)#284,664
table(duplicated(hmPairs$mousePep))#1437
table(duplicated(hmPairsShort$mousePep))#446
hmPairsShort[mousePep=='MEISELNR']

aaMass = fread('../aa.csv')
aaMass = rbind(aaMass, 
               data.table(one_letter_code='X',three_letter_code='Xle',
                          chemical_formula='C6H11ON',Monoisotopic=113.0841, 
                          Average=113.1594,AAtype='unpolar'))


#aaDelta = fread('~/translationError/aaDelta.tsv') #deltaMass is first column minus the second column
aaDelta = list()

for (i in aaMass$one_letter_code){
  aaDelta[[i]] = aaMass[one_letter_code==i, Average] - aaMass$Average
  names(aaDelta[[i]]) = aaMass$one_letter_code
}

aaDelta = do.call('rbind',aaDelta)
aaDelta = as.data.table(aaDelta)
aaDelta$aaTo = colnames(aaDelta)
aaDelta = melt(aaDelta, id.vars = 'aaTo', variable.name = 'aaFrom', value.name = 'detaMass')
aaDelta = aaDelta[detaMass!=0]
aaDelta = aaDelta[!aaTo%in%c('I','L')]
aaDelta = aaDelta[aaFrom!='X']
aaDelta$aaFrom = as.character(aaDelta$aaFrom)
dim(aaDelta) #360   3
fwrite(aaDelta, file = '~/translationError/aaDelta.tsv', sep = '\t')

aaDelta[,min(dist(detaMass)),by='aaFrom'] #0.0434

