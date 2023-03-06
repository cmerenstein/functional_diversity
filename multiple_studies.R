library(tidyverse)
library(curatedMetagenomicData)
library(SummarizedExperiment)
library(vegan)
library(parallel)


## --- list studies & pick -----
studies = curatedMetagenomicData()
studies = studies[grepl("pathabundance", studies)]

case_control = sapply(studies, function(stud){

        d = curatedMetagenomicData(stud, dryrun = F)[[1]]
        meta = pData(d)
        condition = meta$study_condition

        if ( length(unique(condition)) > 1 & nrow(meta) > 50) { return(stud)}
        else {return(NA)}
})
case_control = case_control[ !(is.na(case_control)) ]
write.table(unname(case_control), "case_control.txt", row.names = F, col.names = F)

## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## -------------- Necessary functions -------------------

clean_pathways= function(pathways) {
## get the dataframes of pathways/taxa/samples

    pathways_df =data.frame( pathway = rownames(pathways), pathways)
    tidy = gather(pathways_df, "sample", "abundance", -pathway)

    # everything assigned to a taxa
    tidy_clean = tidy[ !(grepl("UNINTEGRATED", tidy$pathway)) & tidy$pathway!= "UNMAPPED",] %>%
                       separate(pathway, c("PWY", "description", "taxa"), 
                        sep = ": |\\|", fill = "right")

    assigned = filter(tidy_clean, taxa != "unclassified" & !(is.na(taxa)))
    totals = filter(tidy_clean, is.na(taxa))

    return(list(assigned, totals))
}

pathway_bray = function(pwy, assigned, condition){
## takes df of assigned taxa and produces p value based on permanova w/ bray distances

    pathway_long = filter(assigned, PWY == pwy) %>%
                    select(taxa, sample, abundance) %>%
                    spread( sample, abundance)    
    pathway_mat = as.matrix(pathway_long[, 2:ncol(pathway_long)])
    if (nrow(pathway_mat) == 1) { return(NA)}

    ## remove columns with no values, return NA if we no longer have each condition
    pathway_mat = pathway_mat[, colSums(pathway_mat) > 0 , drop = F]
    if (length( unique( condition[ colnames(pathway_mat)])) < 2 |
         ncol(pathway_mat) < 10 ) { return(NA)}
    
    bray = as.matrix(vegdist( t(pathway_mat) ))

    set.seed(19104)
    permanova = adonis2( bray ~ condition[rownames(bray)], permutations = 999 )
    return( permanova[1, "Pr(>F)"] )    
        
}


pathway_kruskal = function(pwy, totals, condition){
## takes df of pathway totals  and produces p value from kruskal test

    pathway_abundance = filter(totals, PWY == pwy) 

    if (length( unique( condition[pathway_abundance$sample])) < 2 |
        nrow(pathway_abundance) < 10) {return(NA)} 

    wx = kruskal.test(pathway_abundance$abundance ~ condition[ pathway_abundance$sample ])
    return(wx$p.value)
}

## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## --- Main loop of studies ----
significance_list = lapply( case_control[1:10], function(stud){
    
    d = curatedMetagenomicData(stud, dryrun = F)[[1]]

    meta = pData(d)
    pathways = exprs(d)

    ## first remove samples where condition is NA
    meta = meta[ !(is.na(meta$study_condition)),]
    pathways = pathways[, colnames(pathways) %in% rownames(meta)]

    ## because R gets funky sometimes, we have to make sure names are 'proper'
    rownames(meta) = make.names(rownames(meta))
    colnames(pathways) = make.names(colnames(pathways))
    condition = meta$study_condition 
    names(condition) = rownames(meta)
    

    clean = clean_pathways(pathways)
    assigned = clean[[1]]
    totals = clean[[2]]

    ## get p values
    adonis_p = sapply(unique(assigned$PWY), pathway_bray, assigned = assigned, condition = condition) 
    kruskal_p = sapply(unique(totals$PWY), pathway_kruskal, totals = totals, condition = condition)
   
    ## make into dataframe and merge
    adonis_df = data.frame(pathway = names(adonis_p), adonis_p, 
                            adonis_FDR = p.adjust(adonis_p, "BH")) 
    kruskal_df = data.frame(pathway = names(kruskal_p), kruskal_p, 
                            kruskal_FDR = p.adjust(kruskal_p, "BH")) 
    pathways_sig = full_join(adonis_df, kruskal_df, by = "pathway")

    pathways_sig$study = stud
    return(pathways_sig)
})

significance_df = do.call("rbind", significance_list)
saveRDS(significance_df, "sig_1_10.rds")

head(arrange(significance_df, adonis_p), 20)

adonis_sig = significance_df %>% filter(adonis_FDR < 0.05)
rownames(adonis_sig) = NULL
table(adonis_sig$pathway) %>% sort() %>% tail(20)

kruskal_sig = significance_df %>% filter(kruskal_FDR < 0.05)
table(kruskal_sig$pathway) %>% sort() %>% tail(20)

