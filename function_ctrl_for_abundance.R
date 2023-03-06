library(tidyverse)
library(curatedMetagenomicData)
library(SummarizedExperiment)
library(vegan)
library(ade4)

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


## Load pathway abundances
d_path = curatedMetagenomicData("JieZ_2017.pathabundance_relab.stool", dryrun = F)[[1]]
meta_path = pData(d_path)
rownames(meta_path) = make.names(rownames(meta_path))
pathways = exprs(d_path)
colnames(pathways) = make.names(colnames(pathways))

pathways_cleaned = clean_pathways(pathways)
assigned = pathways_cleaned[[1]]
totals = pathways_cleaned[[2]]

## load taxa abundances
d_tax = curatedMetagenomicData("JieZ_2017.metaphlan_bugs_list.stool", dryrun = F)[[1]]
meta_tax = pData(d_tax)
rownames(meta_tax) = make.names(rownames(meta_tax))

abundance = exprs(d_tax)
colnames(abundance) = make.names(colnames(abundance))

abundance_s = abundance[ grepl( "s__", rownames(abundance) ) & !(grepl("t__", rownames(abundance))),]
rownames(abundance_s) = sapply(rownames(abundance_s), function(r){
                               split = strsplit(r, "\\|")[[1]]
                               return(paste(split[6], split[7], sep = "."))} ) %>%
                        unname()
abundance_filt = abundance_s[ rownames(abundance_s) %in% assigned$taxa ,]


condition = meta_path$study_condition
names(condition) = rownames(meta_path)
wilcox_p = sapply(unique(totals$PWY), function(pwy){

    pathway_abundance = filter(totals, PWY == pwy)

    if (length( unique( condition[pathway_abundance$sample])) < 2 |
        nrow(pathway_abundance) < 10) {return(NA)}

    wx = wilcox.test(pathway_abundance$abundance ~ condition[ pathway_abundance$sample ])
    return(wx$p.value)

})


## ------------- individual taxa functions ---------------------
common = assigned %>% filter(abundance > 0) %>%
                group_by(PWY, taxa) %>% summarize(n = n()) %>% ungroup() %>%
                filter(n > (.8 * nrow(meta_tax)) )

common$lm_p = apply(common, 1, function(row){

    pwy = row["PWY"]
    taxa = row["taxa"]

    taxa_pathway = assigned[ assigned$PWY == pwy & assigned$taxa == taxa, c("sample", "abundance")]
    taxa_pathway$taxa_total = abundance_s[taxa, taxa_pathway$sample]

    taxa_pathway$condition = meta_tax[taxa_pathway$sample, "study_condition"]

    fit = summary(lm(abundance ~ taxa_total + as.integer(as.factor(condition)), data = taxa_pathway))
    return( coef(fit)[3,4])
})
common$lm_FDR = p.adjust(common$lm_p, "BH")
arrange(common, lm_p) %>% head(20)






## -------------- individual pathways --------------------------

species_per = assigned[, c("PWY", "taxa")] %>% unique() %>%
                select("PWY") %>% table() %>% sort()
sig$species_per = species_per[rownames(sig)]

prevalence = assigned %>% filter(abundance > 1e-6) %>%
                group_by(PWY, sample) %>% summarize(n = n()) %>% ungroup() %>%
                select(PWY, sample) %>%
                group_by(PWY) %>% summarize(n = n())

median_abundance = assigned %>% group_by(PWY, sample) %>%
                    summarize(ab = sum(abundance)) %>%
                    ungroup() %>% group_by(PWY) %>% summarize(median_abundnace = median(ab)) %>%
                    ungroup() %>% as.data.frame()
left_join(sig, median_abundance, by = c("pathway" = "PWY")) %>% arrange(adonis_p) %>% head(40)
left_join(sig, prevalence, by = c("pathway" = "PWY")) %>% 
            filter(n > 20 ) %>% filter(adonis_FDR < .05 & kruskal_p > .1)

pwy = "PWY-6305"
assigned_pwy = assigned[assigned$PWY == pwy,]
assigned_pwy$condition = meta_path[assigned_pwy$sample, "study_condition"]

pdf("figures/pathway_abundance.pdf", height = 8, width = 8)

pwys = sig[ sig$species_per < 20 & sig$adonis_FDR < .01 & sig$kruskal_p > .3, "pathway"]
pwys = pwys[!(is.na(pwys))]

for (pwy in pwys){

    assigned_pwy = assigned[assigned$PWY == pwy,]
    assigned_pwy$condition = meta_path[assigned_pwy$sample, "study_condition"]

    p = filter(assigned_pwy, abundance > 0) %>% 
        ggplot( aes(x = sample, y = abundance, fill = taxa)) + 
            theme_classic() +
            geom_bar(position = "stack", stat = "identity") + 
            facet_grid(~condition, scales = "free_x", space = "free") + 
            ggtitle(pwy) + 
            theme(legend.position = "bottom") + 
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
            guides(fill=guide_legend(ncol=2)) + 
            scale_fill_manual(values = pick_colors(length(unique(assigned_pwy$taxa))))
    print(p)

    p = filter(totals, PWY == pwy) %>%
            left_join(meta_path, by = c("sample" = "subjectID")) %>%
            ggplot(aes(x = study_condition, y = abundance)) +
                theme_classic() + 
                geom_boxplot(lwd = 2) + 
                theme(text = element_text(size = 20))
    print(p)

    ecoli = assigned_pwy[ assigned_pwy$taxa == "g__Escherichia.s__Escherichia_coli",] %>% unique()
    total_ecoli = abundance_s["g__Escherichia.s__Escherichia_coli", ecoli$sample]
    if (length(total_ecoli) > 20) {
        print(cor.test( ecoli$abundance, total_ecoli, method = "spearman"))
        print(plot(ecoli$abundance, total_ecoli) )
    }
}
dev.off()



pathway_mantel = function(pwy, assigned, abundance_bray){
## takes df of assigned taxa and produces p value based on permanova w/ bray distances

    print(pwy) 
    pathway_long = filter(assigned, PWY == pwy) %>%
                    select(taxa, sample, abundance) %>%
                    spread( sample, abundance)    
    if (nrow(pathway_long) < 1){ return(NA)}
    pathway_mat = as.matrix(pathway_long[, 2:ncol(pathway_long)])
    if (nrow(pathway_mat) == 1) { return(NA)}

    ## remove columns with no values, return NA if we no longer have each condition
    pathway_mat = pathway_mat[, colSums(pathway_mat) > 0 , drop = F]
    if (ncol(pathway_mat) < 10){return(NA)}
    
    bray = as.matrix(vegdist( t(pathway_mat) ))

    ## match rownames & columns
    abundance_bray_filt = abundance_bray[rownames(bray), colnames(bray)]

    set.seed(19104)
    mantel = mantel.rtest(as.dist(bray), as.dist(abundance_bray_filt))
    return(mantel$obs )    
        
}
sig$mantel = sapply(sig$pathway, pathway_mantel, 
        abundance_bray = abundance_bray, assigned = assigned)

df = data.frame(ab = abundance_s["g__Bifidobacterium.s__Bifidobacterium_catenulatum", ], 
                cond = meta_tax[colnames(abundance_s), "study_condition"])
wilcox.test(ab ~ cond, data = df)

df = data.frame(ab = abundance_s["g__Klebsiella.s__Klebsiella_pneumoniae", ], 
                cond = meta_tax[colnames(abundance_s), "study_condition"])
wilcox.test(ab ~ cond, data = df)

