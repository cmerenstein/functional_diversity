library(tidyverse)
library(curatedMetagenomicData)
library(SummarizedExperiment)
library(vegan)


## --- list studies & pick -----
studies = curatedMetagenomicData()
studies = studies[grepl("pathabundance", studies)]

case_control = sapply(studies, function(stud){

        d = curatedMetagenomicData(stud, dryrun = F)[[1]]
        meta = pData(d)
        condition = meta$study_condition

        if ( length(unique(condition)) > 1) { return(stud)}
        else {return(NA)}
})
case_control = case_control[ !(is.na(case_control)) ]


get_assigned_taxa = function(pathways) {

    pathways_df =data.frame( pathway = rownames(pathways), pathways)
    long = gather(pathways_df, "sample", "abundance", -pathway)

    # everything assigned to a taxa
    assigned = long[ !(grepl("UNINTEGRATED", long$pathway)) & long$pathway!= "UNMAPPED",]
    assigned = separate(assigned, pathway, c("PWY", "description"), sep = ": ") %>%
                    separate(description, c("descirption", "taxa"), sep = "\\|", fill = "right") %>%
                    filter(taxa != "unclassified" & !(is.na(taxa)))
    return(assigned)
}

get_pathway_totals = function(pathways){

    pathways_df =data.frame( pathway = rownames(pathways), pathways)
    long = gather(pathways_df, "sample", "abundance", -pathway)

    # everything assigned to a taxa
    assigned = long[ !(grepl("UNINTEGRATED", long$pathway)) & long$pathway!= "UNMAPPED",]
    totals = separate(assigned, pathway, c("PWY", "description"), sep = ": ") %>%
                    separate(description, c("descirption", "taxa"), sep = "\\|", fill = "right") %>%
                    filter(is.na(taxa))
    return(totals)

}

## takes df of assigned taxa and produces p value of bray dissimilarity
pathway_bray = function(pwy, assigned, condition){

    print(pwy)
    
    pathway_long = filter(assigned, PWY == pwy) %>%
                    select(taxa, sample, abundance) %>%
                    spread( sample, abundance)    
    pathway_mat = as.matrix(pathway_long[, 2:ncol(pathway_long)])
    if (nrow(pathway_mat) == 1) { return(NA)}

    ## remove columns with no values, return NA if we no longer have each condition
    pathway_mat = pathway_mat[, colSums(pathway_mat) > 0 ]
    if (length( unique( condition[ colnames(pathway_mat)])) < 2 |
         ncol(pathway_mat) < 10 ) { return(NA)}
    
    bray = as.matrix(vegdist( t(pathway_mat) ))

    set.seed(19104)
    permanova = adonis2( bray ~ condition[rownames(bray)], permutations = 9999 )
    return( permanova[1, "Pr(>F)"] )    
        
}


## takes df of pathway totals  and produces p value from wilcox test
pathway_wilcox = function(pwy, totals, condition){

    pathway_abundance = filter(totals, PWY == pwy) 

    if (length( unique( condition[pathway_abundance$sample])) < 2 |
        nrow(pathway_abundance) < 10) {return(NA)} 

    wx = wilcox.test(pathway_abundance$abundance ~ condition[ pathway_abundance$sample ])
    return(wx$p.value)
}

## --- Main loop of studies ----
for (stud in case_control[1:5]){
    
    d = curatedMetagenomicData(stud, dryrun = F)[[1]]

    meta = pData(d)
    pathways = exprs(d)

    ## because R gets funky sometimes, we have to make sure names are 'proper'
    rownames(meta) = make.names(rownames(meta))
    colnames(pathways) = make.names(colnames(pathways))
    condition = meta$study_condition 
    names(condition) = rownames(meta)

    assigned = get_assigned_taxa(pathways)
    totals = get_pathway_totals(pathways)

    adonis_p = sapply(unique(assigned$PWY), pathway_bray, assigned = assigned, condition = condition) 
    wilcox_p = sapply(unique(totals$PWY), pathway_wilcox, totals = totals, condition = condition)
   
    adonis_df = data.frame(pathway = names(adonis_p), adonis_p, adonis_FDR = p.adjust(adonis_p, "BH")) 
    wilcox_df = data.frame(pathway = names(wilcox_p), wilcox_p, wilcox_FDR = p.adjust(wilcox_p, "BH")) 

    print(head(adonis_df))
    print(wilcox_df %>% head())

    pathways_sig = full_join(adonis_df, wilcox_df, by = "pathway")

    cat("adonis: ",  sum(adonis_df$adonis_FDR < .1, na.rm = T), "\n")
    cat("wilcox:", sum(wilcox_df$wilcox_FDR < .1, na.rm = T), "\n")
}






