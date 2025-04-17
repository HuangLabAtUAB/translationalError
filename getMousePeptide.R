library(data.table)
allSymbol = fread('./HOM_MouseHumanSequence.rpt')
#table(allSymbol$`Common Organism Name`)
colnames(allSymbol) = gsub(colnames(allSymbol), pattern='\\s+', replacement='_')
allSymbol[Common_Organism_Name =='mouse, laboratory',Common_Organism_Name:='mouse']
allSymbol=allSymbol[SWISS_PROT_IDs !='']
allSymbol[, combineID := paste0(Common_Organism_Name,',',
                                `Symbol`,',',SWISS_PROT_IDs)]
#built mouse-human corresponding relationship
allSymbolMerge = allSymbol[,.(Symbol=paste(combineID, collapse = ':')), by='DB_Class_Key']
allSymbolMerge = allSymbolMerge[grepl(Symbol, pattern='human')&grepl(Symbol, pattern='mouse')]
dim(allSymbolMerge)#16603     2
allSymbol = allSymbol[DB_Class_Key%in%allSymbolMerge$DB_Class_Key]
allSymbol[grepl(Symbol, pattern=',')]
allSymbol[grepl(SWISS_PROT_IDs, pattern=',')]

humanProList = list()
mouseProList = list()

for (dbKey in allSymbolMerge$DB_Class_Key){
  mSymbols = allSymbol[DB_Class_Key==dbKey&Common_Organism_Name=='mouse',unique(Symbol)]
  hSymbols = allSymbol[DB_Class_Key==dbKey&Common_Organism_Name=='human',unique(Symbol)]
  humanIDs = allSymbol[DB_Class_Key==dbKey&Common_Organism_Name=='human',SWISS_PROT_IDs]
  mouseIDs = allSymbol[DB_Class_Key==dbKey&Common_Organism_Name=='mouse',SWISS_PROT_IDs]
  if(length(humanIDs)>1) humanIDs=paste(humanIDs, collapse = ',')
  if(length(mouseIDs)>1) mouseIDs=paste(mouseIDs, collapse = ',')
  if(grepl(humanIDs, pattern=',')) humanIDs=unlist(strsplit(humanIDs, split=','))
  if(grepl(mouseIDs, pattern=','))  mouseIDs=unlist(strsplit(mouseIDs, split=','))
    
  for (humanID in humanIDs){
    tmp = proteins(edb, 
             filter = UniprotFilter(humanID),
             return.type = "AAStringSet")
    humanProList[[hSymbols[1]]] = c(humanProList[[hSymbols[1]]],tmp)
  }
  for (mouseID in mouseIDs){
    tmp = proteins(medb, 
             filter = UniprotFilter(mouseID),
             return.type = "AAStringSet")
    mouseProList[[mSymbols[1]]] = c(mouseProList[[mSymbols[1]]],tmp)
  }
  if(length(hSymbols)>1) for(sym in hSymbols[2:length(hSymbols)]) {humanProList[[sym]]=humanProList[[hSymbols[1]]]}
  if(length(mSymbols)>1) for(sym in mSymbols[2:length(mSymbols)]) {mouseProList[[sym]]=mouseProList[[mSymbols[1]]]}
}
length(humanProList)#16794

humanPro = lapply(humanProList, function(x)lapply(x, as.character))
humanPro = lapply(humanPro, function(x){data.table(proID=names(unlist(x)), seqs= unlist(x))})
humanPro = humanPro[sapply(humanPro, nrow)>0]
length(humanPro) #16605
humanPro = rbindlist(humanPro,idcol = 'symbol')
humanPro = humanPro[!duplicated(proID)]
dim(humanPro) #39919
table(humanPro$symbol %in% allSymbol$Symbol)

mousePro = lapply(mouseProList, function(x)lapply(x, as.character))
mousePro = lapply(mousePro, function(x){data.table(proID=names(unlist(x)), seqs= unlist(x))})
mousePro = mousePro[sapply(mousePro, nrow)>0]
length(mousePro) #15527
mousePro = rbindlist(mousePro,idcol = 'symbol')
mousePro = mousePro[!duplicated(proID)]
dim(mousePro) #23130


####
mousePro[,peptides:=sapply(seqs, function(x)cleaver::cleave(x, enzym = "trypsin", 
                                                 missedCleavages = 0:2,
                                                 custom = NULL, unique = TRUE))]

mousePro[,peptides:=sapply(peptides, function(x){aaNo=sapply(x, nchar); x[aaNo>6 & aaNo < 61]})]


humanPro[,peptides:=sapply(seqs, function(x)cleaver::cleave(x, enzym = "trypsin", 
                                                            missedCleavages = 0:2,
                                                            custom = NULL, unique = TRUE))]
humanPro[,peptides:=sapply(peptides, function(x){aaNo=sapply(x, nchar); x[aaNo>6 & aaNo < 61]})]



##
#table(humanPro$symbol %in% allSymbol$Symbol)
#table(mousePro$symbol %in% allSymbol$Symbol)
##
allSAVpep = list()
for (dbKey in allSymbolMerge$DB_Class_Key){
  hS = allSymbol[DB_Class_Key==dbKey&Common_Organism_Name=='human', Symbol]
  mS = allSymbol[DB_Class_Key==dbKey&Common_Organism_Name=='mouse', Symbol]
  if( any(hS%in% humanPro$symbol) & 
      any(mS%in% mousePro$symbol)){
    hS = intersect(hS, humanPro$symbol)
    mS = intersect(mS, mousePro$symbol)
    #by using %in%, i want to be more inclusive
    hPeps = humanPro[symbol%in%hS, unique(unlist(peptides))]
    mPeps = mousePro[symbol%in%mS, unique(unlist(peptides))]
    pepDis = stringdist::stringdistmatrix(a=hPeps,b=mPeps,weight = c(d=1,i=1,s=0.1,t=1))
    idices = which(pepDis==0.1,arr.ind = T)
    allSAVpep[[as.character(dbKey)]]=data.table(humanPep=hPeps[idices[,'row']], mousePep=mPeps[idices[,'col']], 
                                  humanGene=paste(hS, collapse = ':'), mouseGene=paste(mS, collapse = ':'))
  }
}

allSAVpepDt = rbindlist(allSAVpep, idcol = 'dbKey')
allSAVpepDt = allSAVpepDt[!is.na(humanPep)]
dim(allSAVpepDt)

allSAVpepDt[,table(grepl(humanGene,pattern=':'))]
fwrite(allSAVpepDt, file = 'allTrypsinPeptidePairsWithOneSAV_2maxMissCleavge.csv')

