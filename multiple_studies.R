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
    permanova = adonis2( bray ~ condition[rownames(bray)], permutations = 9999 )
    return( permanova[1, "Pr(>F)"] )    
        
}


## takes df of pathway totals  and produces p value from wilcox test
pathway_wilcox = function(pwy, totals, condition){

    print(pwy)

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

    print( sum(adnois_df$FDR < .1))
    print( sum(wilcox_df$FDR < .1))
}


d = curatedMetagenomicData("YeZ_2018.pathabundance_relab.stool", dryrun = F)[[1]] 
meta = pData(d)
pathways = exprs(d)

## because R gets funky sometimes, we have to make sure names are 'proper'
rownames(meta) = make.names(rownames(meta))
colnames(pathways) = make.names(colnames(pathways))


pathways_df =data.frame( pathway = rownames(pathways), pathways)
long = gather(pathways_df, "sample", "abundance", -pathway)

# everything assigned to a taxa
assigned = long[ !(grepl("UNINTEGRATED", long$pathway)) & long$pathway!= "UNMAPPED",]
assigned = separate(assigned, pathway, c("PWY", "description"), sep = ": ") %>%
                separate(description, c("descirption", "taxa"), sep = "\\|", fill = "right") %>%
                filter(taxa != "unclassified" & !(is.na(taxa)))

## get study condition
condition = meta$study_condition
names(condition) = rownames(meta)

adonis_p = sapply(unique(assigned$PWY), pathway_bray, assigned = assigned, condition = condition) 
adonis_p = adonis_p[ !(is.na(adonis_p)) ]
adonis_df = data.frame(pathway = names(adonis_p), p = adonis_p, fdr = p.adjust(adonis_p, "BH"))


wilcox_p = sapply(unique(totals$PWY), pathway_wilcox, totals = totals, condition = condition)



totals = assigned[ is.na(assigned$taxa), ]

## taxa per pathway per sample
taxa_richness = filter(assigned, !(is.na(taxa) & taxa == "unclassified")) %>%
                    group_by(sample, PWY) %>% summarize(n_taxa = n())

## taxa per pathway across samples
unq = unique(assigned[ !(is.na(assigned$taxa)) & assigned$taxa != "unclassified" ,c("PWY", "taxa")])
unq_richness = group_by(unq, PWY) %>% summarize(n_taxa = n()) %>% ungroup() %>% as.data.frame()

pdf("figures/taxa_per_pathway.pdf")

hist(taxa_richness$n_taxa, breaks=45,  main = "taxa per pathway per sample")

hist(unq_richness$n_taxa, breaks = 45, main ="taxa per pathway across all samples")

dev.off()



## significant differences
all_totals = spread(totals[,c("sample", "PWY", "abundance")] , "PWY", "abundance")
totals_mat = all_totals[,2:ncol(all_totals)]

rownames(totals_mat) = gsub("\\.", "-", all_totals$sample) # fix rownames

condition = meta[rownames(totals_mat), "study_condition"] # get study condition

wilcox_p = apply(totals_mat, 2, function(d){

        if ( (sum(d > 0) / length(d) ) < 0.1){ return(NA)}
        else {
            wx = wilcox.test( d ~ condition)
            return(wx$p.value)
        }

})
wilcox_p = wilcox_p[!(is.na(wilcox_p))]
wilcox_df = data.frame(pathway = names(wilcox_p), p = wilcox_p, fdr = p.adjust(wilcox_p, "BH"))






