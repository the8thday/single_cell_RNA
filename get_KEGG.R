

library("KEGGREST")
library("EnrichmentBrowser") #这个包里面的一些函数会调用KEGGREST里面的函数


### download the pathways
hsapathway <- downloadPathways("hsa") #只有在第一次运行这句代码时，耗时较长

### retrieve gene sets
hsa <- getGenesets(org = "hsa", db = "kegg", 
                   gene.id.type = "SYMBOL",cache = TRUE, 
                   return.type="list") ##只有在第一次运行这句代码时，耗时较长

writeGMT(hsa, gmt.file = "20230205_kegg_hsa.gmt")



# REST --------------------------------------------------------------------


gsInfo = KEGGREST::keggGet('hsa04110')[[1]]
names(gsInfo)

geneSetRaw = sapply(strsplit(gsInfo$GENE, ";"), function(x) x[1])
geneSet = list(geneSetRaw[seq(2, length(geneSetRaw), 2)])
names(geneSet) = gsInfo$NAME


getKEGG <- function(ID){
  
  library("KEGGREST")
  
  gsList = list()
  for(xID in ID){
    
    gsInfo = keggGet(xID)[[1]]
    if(!is.null(gsInfo$GENE)){
      geneSetRaw = sapply(strsplit(gsInfo$GENE, ";"), function(x) x[1])
      xgeneSet = list(geneSetRaw[seq(2, length(geneSetRaw), 2)])
      NAME = sapply(strsplit(gsInfo$NAME, " - "), function(x) x[1])
      names(xgeneSet) = NAME
      gsList[NAME] = xgeneSet
    } else{
      cat(" ", xID, "No corresponding gene set in specific database.\n")
    }
  }
  return(gsList)
}



library("rjson")
download.file("http://togows.dbcls.jp/entry/pathway/hsa04930/genes.json", "hsa04930.json")
json = fromJSON(file ="hsa04930.json")
geneSet = list(as.character(sapply(json[[1]], function(x) sapply(strsplit(x[1], ";"), 
                                                                 function(x) x[1]))))




### by bash

curl -s https://rest.kegg.jp/list/pathway/hsa -o hsa.ko2pathway
sed -i 's# - Homo sapiens (human)##' hsa.ko2pathway
head -5 hsa.ko2pathway



curl -s http://togows.dbcls.jp/entry/pathway/hsa00010/genes.json \
| awk 'BEGIN{OFS="\t"}$0~/KO:/{
        match($0,/"([^"]+)": "([^;[]+);?(.+?) \[(KO:[^\]]+)/,a);
        a[3]=a[3]==""?"-":a[3];
        print a[1],a[2],a[3],a[4]
        }' \
| head -10


cat hsa.ko2pathway| while IFS=$'\t' read -r id name
do
gene_set=`curl -s http://togows.dbcls.jp/entry/pathway/${id}/genes.json \
| awk 'BEGIN{OFS="\t"}$0~/KO:/{match($0,/"([^"]+)": "([^;[]+);?(.+?) \[(KO:[^\]]+)/,a); print a[1]}' \
| paste -d"\t" -s`
echo -e "${id}\t${name}\t${gene_set}" >> kegg.hsa.gmt
done


# ggKEGG ------------------------------------------------------------------

library(ggkegg)



# biomart -----------------------------------------------------------------


library(biomaRt)


ensembl <- useEnsembl(biomart = "ensembl")
datasets <- listDatasets(ensembl)
head(datasets)
searchDatasets(mart = ensembl, pattern = "btaurus")


# useast
# https://jul2023.archive.ensembl.org
# listEnsemblArchives()
# Bovine, Bovidae
human <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", 
                    # mirror = 'www',
                    # host = 'https://jan2019.archive.ensembl.org'
                    )
mouse <- useEnsembl("ensembl", dataset = "mmusculus_gene_ensembl",
                    host = 'https://jan2019.archive.ensembl.org'
                    )

# (ARS-UCD1.2) (29472)
cow <- useEnsembl('genes', dataset = 'btaurus_gene_ensembl',
                  # mirror = 'www',
                  # host = 'https://jan2019.archive.ensembl.org'
)

searchAttributes(mart = cow, pattern = "_gene_name")

# getBM(c("ensembl_gene_id", 'chromosome_name'), 'hgnc_symbol','BTG1', cow)
getBM(c('chromosome_name', 'external_gene_name'), 
      'ensembl_gene_id', 
      'ENSBTAG00000001069', cow)

conv <- getLDS(attributes = c("hgnc_symbol","entrezgene_id", "ensembl_gene_id"),
       filters = "hgnc_symbol", 
       values = c("TP53", 'BTG1'),
       mart = human,
       attributesL = c("external_gene_name"), 
       martL = human,
       verbose = FALSE
)

conv



#' Gene between Human and Mouse
#'
#' convert human genes to mouse genes by LDS, or vice verse.
#'
#' @param genes genes to convert.
#' @param geneid the gene id type, one of symbol, ensembl and entrez, default is symbol.
#' @param invert bool value for conversion, if TRUE, convert human to mouse; if FALSE, convert mouse to human.
#' @param host Default is the current site for ensembl. For archives, add the prefix.
#' the latest GRCh38/GRCm38 is https://aug2020.archive.ensembl.org (version 101)
#' @importFrom biomaRt useMart getLDS
#' @return A vector of genes
#' @export
convertHumanMouse <- function(genes, geneid = c("symbol", "ensembl", "entrez"), 
                              invert = FALSE, 
                              host = "https://www.ensembl.org") {
  geneid <- match.arg(geneid)
  human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = host)
  mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl",host = host)
  if (geneid == "symbol") {
    mouseid <- "mgi_symbol"
    humanid <- "hgnc_symbol"
  } else if (geneid == "ensembl") {
    mouseid <- humanid <- "ensembl_gene_id"
  } else if (geneid == "entrez") {
    mouseid <- humanid <- "entrezgene_id"
  }
  
  if (invert) {
    conversion <- getLDS(
      attributes = humanid, filters = humanid, values = genes, mart = human,
      attributesL = mouseid, martL = mouse, uniqueRows = TRUE
    )
  } else {
    conversion <- getLDS(
      attributes = mouseid, filters = mouseid, values = genes, mart = mouse,
      attributesL = humanid, martL = human, uniqueRows = TRUE
    )
  }
  
  #unique(conversion[, 2])
  conversion
}


convertHumanMouse(genes = 'TP53',
                  geneid = 'symbol',
                  invert = TRUE
                  )











