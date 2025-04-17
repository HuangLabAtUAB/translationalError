library(ggalluvial)
##to analyze the variation across different SAVs

## AA annotation
aaMass[, AAtype := '']
aaMass[one_letter_code%in%c('R','H','K'), AAtype := 'positive']
#lysine, arginine, and histidine
aaMass[one_letter_code%in%c('D','E'), AAtype := 'negative']
#aspartic acid and glutamic acid
aaMass[one_letter_code%in%c('S','T','N','Q','C','P'), AAtype := 'polar']
#asparagine, cysteine, glutamine, proline, serine, and threonine
aaMass[one_letter_code%in%c('G','A','V','I','L','M','F'), AAtype := 'unpolar']
#glycine, alanine, isoleucine, leucine, methionine, and valine
aaMass[one_letter_code%in%c('F','W','Y'), AAtype := 'aromatic']
#phenylalanine, tyrosine, and tryptophan

#AA sav freq
hmPairs[,humanAAtype:=aaMass[hmPairs$humanAA,AAtype,on='one_letter_code']]
hmPairs[,mouseAAtype:=aaMass[hmPairs$mouseAA,AAtype,on='one_letter_code']]
hmPairs[, sametypeConvert:=ifelse(humanAAtype==mouseAAtype,'sameType','diffType')]
table(hmPairs$sametypeConvert)
# diffType sameType 
# 121,540     170,385 
hmPairs[, convertType:=paste0(humanAAtype,'->',mouseAAtype)]
table(hmPairs$convertType)
# 164115/nrow(hmPairs)
# [1] 0.5745217

  
  
#tcga mutation profile
tcgaMut = fread('./frequent-mutations.2024-06-03.tsv')
dim(tcgaMut)#117,196
tcgaMut = tcgaMut[grepl(protein_change, pattern='\\s[A-Z][0-9]+[A-Z]$')]
dim(tcgaMut)#114722
tcgaMut[, fromAA := gsub(protein_change, pattern='.*\\s+([A-Z])[0-9]+.*', replacement='\\1')]
tcgaMut[, toAA := gsub(protein_change, pattern='.*\\s+[A-Z][0-9]+([A-Z])$', replacement='\\1')]
tcgaMut[, fromAAtype:=aaMass[tcgaMut$fromAA,AAtype,on='one_letter_code']]
tcgaMut[, SAV := paste0(fromAA, '->', toAA)]
tcgaMut[, toAAtype:=aaMass[tcgaMut$toAA,AAtype,on='one_letter_code']]
tcgaMut[, sametypeConvert:=ifelse(fromAAtype==toAAtype,'sameType','diffType')]
tcgaMut[, convertType:=paste0(fromAAtype,'->',toAAtype)]
table(tcgaMut$convertType)
table(tcgaMut$sametypeConvert)
# diffType sameType 
# 79809      34913 
#34913/nrow(tcgaMut) 0.304327

#count matrix for band plot
hmPairsCount = hmPairs[,.N,by=c('humanAAtype','mouseAAtype','convertType','sametypeConvert')]
tcgaMutCount = tcgaMut[,.N,by=c('fromAAtype','toAAtype','convertType','sametypeConvert')]
hmPairsCount[,perc:=N/sum(N)]
tcgaMutCount[,perc:=N/sum(N)]

ggplot(hmPairsCount,aes(y = perc, axis1 = humanAAtype, axis2 = mouseAAtype)) +
  geom_alluvium(aes(fill = sametypeConvert), width = 1/12) +
  geom_stratum(width = 1/12) +
  scale_x_discrete(limits = c("human AA", "mouse AA"), expand = c(.1, .1))+
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_fill_brewer(type = "qual", palette = "Set1") +
  ggtitle("human mouse SAV")


ggplot(tcgaMutCount,aes(y = perc, axis1 = fromAAtype, axis2 = toAAtype)) +
  geom_alluvium(aes(fill = sametypeConvert), width = 1/12) +
  geom_stratum(width = 1/12) +
  scale_x_discrete(limits = c("human AA", "mouse AA"), expand = c(.1, .1))+
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_fill_brewer(type = "qual", palette = "Set1") +
  ggtitle("TCGA mut SAV")

fisher.test(matrix(c(table(tcgaMut$sametypeConvert)[1],
                     nrow(tcgaMut),
                     table(hmPairs$sametypeConvert)[1],
                     nrow(hmPairs)),nrow = 2))
###
length(unique(hmPairs$sav)) #366
length(unique(tcgaMut$SAV)) #150

#plot the top SAV types
hmPairsCountBySAv = hmPairs[,.N,by=c('sav','humanAAtype','mouseAAtype','convertType','sametypeConvert')]
hmPairsCountBySAv = hmPairsCountBySAv[order(-N)]

tcgaMutCountBySAv = tcgaMut[,.N,by=c('SAV','fromAAtype','toAAtype','convertType','sametypeConvert')]
tcgaMutCountBySAv = tcgaMutCountBySAv[order(-N)]

ggplot(hmPairsCountBySAv,aes(sav))+
  stat_bin(aes(y=cumsum(..N..)),geom="step")
  
ggplot(hmPairsCountBySAv[1:20],aes(sav,N))+
  geom_bar(stat = 'identity', aes(
    #colour =convertType, 
                                                                           fill = sametypeConvert))+
  scale_x_discrete(limit=hmPairsCountBySAv$sav[1:20])+
  theme(axis.text.x = element_text(angle = -90))

ggplot(tcgaMutCountBySAv[1:20],aes(SAV,N))+
  geom_bar(stat = 'identity', 
           aes(
             #colour =convertType, 
                                                                           fill = sametypeConvert))+
  scale_x_discrete(limit=tcgaMutCountBySAv$SAV[1:20])+
  theme(axis.text.x = element_text(angle = -90))
                   