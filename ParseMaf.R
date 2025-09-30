library(maftools)
library(data.table)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38)
# 1. Load MAF file
maf <- read.maf("./PDAC_somatic_mutation.maf.tsv")

# 2. Extract relevant mutation data
mutations <- maf@data[, .(Hugo_Symbol, Chromosome,tx, Start_Position, Reference_Allele, 
                          Tumor_Seq_Allele2, Variant_Classification, Protein_Change)]
mutations[,tx:=gsub(pattern = '\\..*',replacement = '',tx)]
table(mutations$Variant_Classification)
mutations = mutations[Variant_Classification%in%
                       c('Nonsense_Mutation','Missense_Mutation','In_Frame_Ins','In_Frame_Del')]
dim(mutations) #4252    8
# 3. Retrieve reference protein sequences (example using Ensembl or UniProt API)
# For simplicity, assume a reference FASTA file with protein sequences
ref_proteins <- readAAStringSet("./Homo_sapiens.GRCh38.pep.all.fa")  # e.g., from UniProt
# Example format: >TP53|ENSP00000269305, >BRCA1|ENSP00000493795, etc.

mutations[,txMut :=paste0(tx, '_', Protein_Change)]
length(unique(mutations$txMut)) #4130



# 5. Generate altered sequences
altered_sequences <- list()
for (i in 1:nrow(mutations)) {
  gene = mutations$Hugo_Symbol[i]
  tx <- mutations$tx[i]
  protein_change <- mutations$Protein_Change[i]
  variant_class <- mutations$Variant_Classification[i]
  
  # Find reference protein sequence (match by gene name)
  protein_seq <- ref_proteins[grep(tx, names(ref_proteins))]
  if(length(protein_seq)>1) stop('duplicated tx id!')
  if (length(protein_seq) == 0) {
    warning("No reference sequence for ", gene)
    next
  }
  # Apply mutation
  altered_seq <- apply_mutation(protein_seq, protein_change, variant_class)
  
  # Store with identifier
  altered_sequences[[paste0(gene, "_",tx,"_", protein_change)]] <- altered_seq
}

fil = sapply(altered_sequences, function(x)any(is.na(x)))
altered_sequence_total = unlist(altered_sequences[!fil])
for (idx in 1:length(altered_sequence_total)){
  if(class(altered_sequence_total[[idx]])=='AAStringSet')
    altered_sequence_total[[idx]] = altered_sequence_total[[idx]][[1]]
}
# 6. Save altered sequences as FASTA
writeXStringSet(AAStringSet(altered_sequence_total), 
                "pdac_cptac_altered_proteins.fasta")
