library(biomaRt)
library(data.table)
library(GenomicRanges)

# setup Ensembl biomaRt
ensemblHost = "asia.ensembl.org"
listMarts()
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
listDatasets()
filters = listFilters(ensembl) # band_start, band_end
head(filters)
attribs <- listAttributes(ensembl)
attribs[grep("band", attribs$description),]

# local functions
cytobandSearch <- function(chromosome, band_start, band_end, mart, ret = "gr"){
  dt <- as.data.table(getBM(mart = mart, attributes=c('ensembl_gene_id', 'hgnc_symbol','chromosome_name',
                                   'start_position','end_position', 'band'),
                      filters = c('chromosome_name', 'band_start', 'band_end'),
                      values = list(chromosome, band_start, band_end)))
  if (ret == "gr") {
    gr <- makeGRangesFromDataFrame(dt, 
                                   start.field = 'start_position',
                                   end.field = 'end_position',
                                   keep.extra.columns = T)
    return(gr) } else {
      return(dt)
    }
}

# load DNA Damage Reponse genes
gois <- fread("https://raw.githubusercontent.com/skurscheid/GeneSets/master/Literature/Pearl_et_al_2015_gois/gois_Pearl_et_al_2015.csv")
length(unique(gois$GeneID))

gois_ensembl <- as.data.table(getBM(mart = ensembl, 
                                         attributes=c('ensembl_gene_id', 'hgnc_symbol','chromosome_name',
                                                      'start_position','end_position', 'band'),
                                         filter = "hgnc_symbol",
                                         values = unique(gois$GeneID)))

gr.gois <- makeGRangesFromDataFrame(gois_ensembl, 
                                         start.field = 'start_position', 
                                         end.field = 'end_position', 
                                         keep.extra.columns = T)

histone_genes <- fread("https://raw.githubusercontent.com/skurscheid/GeneSets/master/DBSearches/histoneGenes.txt")
histone_genes_ensembl <- as.data.table(getBM(mart = ensembl, 
                                         attributes=c('ensembl_gene_id', 'hgnc_symbol','chromosome_name',
                                                      'start_position','end_position', 'band'),
                                         filter = "hgnc_symbol",
                                         values = unique(histone_genes$ApprovedSymbol)))
gr.histone_genes <- makeGRangesFromDataFrame(histone_genes_ensembl, 
                                         start.field = 'start_position', 
                                         end.field = 'end_position', 
                                         keep.extra.columns = T)

epigentic_modifier_genes <- fread("https://raw.githubusercontent.com/skurscheid/GeneSets/master/Literature/Nanda_et_al_2016/exp_cnv_final.csv")
epigentic_modifier_genes_ensembl <- as.data.table(getBM(mart = ensembl, 
                                             attributes=c('ensembl_gene_id', 'hgnc_symbol','chromosome_name',
                                                          'start_position','end_position', 'band'),
                                             filter = "hgnc_symbol",
                                             values = unique(epigentic_modifier_genes$V1)))
gr.epigentic_modifier_genes <- makeGRangesFromDataFrame(epigentic_modifier_genes_ensembl, 
                                             start.field = 'start_position', 
                                             end.field = 'end_position', 
                                             keep.extra.columns = T)

#chromatin_associated_gene_products <- fread("https://raw.githubusercontent.com/LixinyuLiu/Tissue-and-cell-type-specific-expression-patterns-of-histone-and-other-chromatin-associated-genes/master/Biogridderived.csv")

gr.gois <- c(gr.gois, gr.histone_genes, gr.epigentic_modifier_genes)

# cytoband information extracted from Marella et al. 2009 doi:10.1158/0008-5472.CAN-09-0420
# MCF10A gains
mcf10a.gains <- GRangesList()
mcf10a.gains[['mcf10a.chr5gains']] <- cytobandSearch('5', 'q23.1', 'q35.3', ensembl)
mcf10a.gains[['mcf10a.chr8gains']] <- cytobandSearch('8', 'p23.3', 'q24.3', ensembl) 
mcf10a.gains[['mcf10a.chr13gains']] <- cytobandSearch('13', 'q32.1', 'q32.2', ensembl)
mcf10a.gains[['mcf10a.chr19gains']] <- cytobandSearch('19', 'q13.11', 'q13.43', ensembl)
mcf10a.gains

# MCF10A losses
mcf10a.losses <- GRangesList()
mcf10a.losses[['mcf10a.chr3loss']] <- cytobandSearch('3', 'p26.3', 'p26.3', ensembl)
mcf10a.losses[['mcf10a.chr9loss']] <- cytobandSearch('9', 'p21.3', 'p21.3', ensembl)
mcf10a.losses[['mcf10a.chr16loss']] <- cytobandSearch('16', 'p11.2', 'p11.2', ensembl)
mcf10a.losses[['mcf10a.chr21loss']] <- cytobandSearch('21', 'p11.2', 'q11.2', ensembl)
mcf10a.losses[['mcf10a.chr22loss']] <- cytobandSearch('22', 'q11.1', 'q11.1', ensembl)
mcf10a.losses

# MCF10AT1 gains
mcf10at1.gains <- GRangesList()
mcf10at1.gains[['mcf10at1.chr3gains']] <- c(cytobandSearch('3', 'p14.3', 'p14.3', ensembl), cytobandSearch('3', 'q13.31', 'q13.31', ensembl))
mcf10at1.gains[['mcf10at1.chr9gains']] <- c(cytobandSearch('9', 'p24.3', 'p11.2', ensembl), cytobandSearch('9', 'q12', 'q34.3', ensembl))
mcf10at1.gains[['mcf10at1.chr10gains']] <- cytobandSearch('10', 'q22.1', 'q22.2', ensembl)
mcf10at1.gains[['mcf10at1.chr16gains']] <- cytobandSearch('16', 'q23.3', 'q23.3', ensembl)
mcf10at1.gains[['mcf10at1.chr17gains']] <- cytobandSearch('17', 'p11.2', 'p11.2', ensembl)
mcf10at1.gains

# MCF10AT1 losses
mcf10at1.losses <- GRangesList()
mcf10at1.losses[['mcf10at1.chr5losses']] <- c(cytobandSearch('5', 'q12.1', 'q12.1', ensembl), cytobandSearch('5', 'q14.3', 'q15', ensembl))
mcf10at1.losses[['mcf10at1.chr8losses']] <- cytobandSearch('8', 'p23.3', 'q24.3', ensembl)
mcf10at1.losses[['mcf10at1.chr15losses']] <- cytobandSearch('15', 'q21.1', 'q21.1', ensembl)
mcf10at1.losses

# MCF10CA1a gains
mcf10ca1a.gains <- GRangesList()
mcf10ca1a.gains[['mcf10ca1a.chr2gains']] <- cytobandSearch('2', 'p25.3','q21.2', ensembl)
mcf10ca1a.gains[['mcf10ca1a.chr3gains']] <- cytobandSearch('3', 'p14.1', 'q29', ensembl)
mcf10ca1a.gains[['mcf10ca1a.chr9gains']] <- c(cytobandSearch('9', 'p24.3', 'p11.2', ensembl), cytobandSearch('9', 'q34.3', 'q34.13', ensembl))
mcf10ca1a.gains[['mcf10ca1a.chr10gains']] <- cytobandSearch('10', 'q11.1', 'q26.3', ensembl)
mcf10ca1a.gains[['mcf10ca1a.chr17gains']] <- cytobandSearch('17', 'p11.2', 'p11.2', ensembl)
mcf10ca1a.gains[['mcf10ca1a.chr20gains']] <- cytobandSearch('20', 'p13', 'q13.33', ensembl)
mcf10ca1a.gains

# MCF10CA1a losses
mcf10ca1a.losses <- GRangesList()
mcf10ca1a.losses[['mcf10ca1a.chr2losses']] <- cytobandSearch('2', 'q21.2', 'q22.1', ensembl, ret = "gr")
mcf10ca1a.losses[['mcf10ca1a.chr5losses']] <- c(cytobandSearch('5', 'q12.1', 'q12.1', ensembl), cytobandSearch('5', 'q14.3', 'q15', ensembl))
mcf10ca1a.losses[['mcf10ca1a.chr8losses']] <- cytobandSearch('8', 'p23.3', 'q24.3', ensembl) 
mcf10ca1a.losses[['mcf10ca1a.chr16losses']] <- cytobandSearch('16', 'q23.1', 'q23.1', ensembl)
mcf10ca1a.losses

# overlap between gois genes and CNV
gois.gains <- c(subsetByOverlaps(gr.gois, mcf10a.gains), 
               subsetByOverlaps(gr.gois, mcf10at1.gains),
               subsetByOverlaps(gr.gois, mcf10ca1a.gains))
gois.gains.GeneSymbols <- unique(gois.gains$hgnc_symbol)
gois.gains.EnsemblIDs <- unique(gois.gains$ensembl_gene_id)

gois.losses <- c(subsetByOverlaps(gr.gois, mcf10a.losses),
                subsetByOverlaps(gr.gois, mcf10at1.losses),
                subsetByOverlaps(gr.gois, mcf10ca1a.losses))
gois.losses.GeneSymbols <- unique(gois.losses$hgnc_symbol)
gois.losses.EnsemblIDs <- unique(gois.losses$ensembl_gene_id)



