input_dir <- "/pfs/downloadGBMData/"
cell_dir <- "/pfs/getGBMCellData/"
out_dir <- "/pfs/out/" 
# functions <- "/pfs/getGBMCellData/functions.R"

# input_dir <- "~/Documents/pfs/downloadGBMData/"
# cell_dir <- "~/Documents/pfs/getGBMCellData/"
# out_dir <- "~/Documents/pfs/getGBMMutation/" 
functions <- "https://github.com/BHKLAB-Pachyderm/getGBMCellData/raw/main/functions.R"

source(functions)
load(paste0(input_dir, "Ensembl.v99.annotation.RData"))
cell_annotation <- readRDS(paste0(cell_dir, "cell.rds"))

# ============= Mutation data =============
#Assay data
mutation<- read.delim(paste(input_dir, "HGCC_WES_mutations_variants.txt", sep=""), header=T, sep="\t", stringsAsFactors = FALSE,na.strings=c("", " ", "NA","NaN"))

mutation$CELL_LINE<- paste("U", mutation $CELL_LINE, sep="")
assay_mut<-data.frame(matrix("wt", nrow = length(unique(mutation$Gene.refGene[!is.na(mutation$Gene.refGene)])), ncol = length(unique(mutation$CELL_LINE))))
rownames(assay_mut)<-unique(mutation$Gene.refGene[!is.na(mutation$Gene.refGene)])
colnames(assay_mut)<-unique(mutation$CELL_LINE)

#Labeling wild type cells as "wt"
for (i in 1:nrow(assay_mut)){
  if(i%%100 ==0 ){print(i)}
  gene<-rownames(assay_mut)[i]
  
  indS<-which(mutation$Gene.refGene == gene)
  for (ind in indS){
    cell <- mutation$CELL_LINE [ind]
    func <- mutation$ExonicFunc.refGene [ind]
    
    if (!is.na(func)){
      if(assay_mut[gene, cell]=="wt"){assay_mut[gene, cell]=func}
      else if(grepl(func, assay_mut[gene, cell], perl=TRUE)) {next} #If the mutation effect is already reported for a gene then we skip
      else{
        assay_mut[gene, cell] = paste(assay_mut[gene,cell],func, sep="///")
      }
    } 
  }
}

assay_mut<-assay_mut[,-2]

feat_mutation<-fdata_builder(annotation=features_gene, assay=assay_mut)
phen_mutation<-ph_data_builder(annotation=cell_annotation, assay=assay_mut)
phen_mutation$Patient_id[phen_mutation$cellid=="U3067"]<-"U3067MG"# Based on "GSM_map" data frame

#Creating ExpressionSet 
assay_mut<-assay_mut[,rownames(phen_mutation)]#rearranging the colnames so it is similar to pheno data
mutation_eSet<- ExpressionSet(assayData = as.matrix(assay_mut), phenoData = AnnotatedDataFrame(phen_mutation), featureData = AnnotatedDataFrame(feat_mutation)) 
print("Mutation: done")

mutation_SE <- eSetToSE(mutation_eSet,annot_name="mut")
saveRDS(mutation_SE, paste0(out_dir, "mutation_SE.rds"))
saveRDS(phen_mutation, paste0(out_dir, "phen_mutation.rds"))