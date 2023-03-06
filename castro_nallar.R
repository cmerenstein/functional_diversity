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

pick_colors = function(n){

    n_hues = ceiling(n / 3)
    hues = seq(0, .99, 1/n_hues)
    colors = c()
    for (hue in hues){
        ## add three colors for each hue
        c1 = hsv(hue, .25, .9)
        c2 = hsv(hue, .85, .9)
        c3 = hsv(hue, .85, .5)

        colors = c(colors, c(c1, c2, c3)) ## not the right way to add to a vector, but fine here
    }
    return(colors)
}


## get pathway significance
sig = readRDS("sig_1_10.rds")
sig = sig[sig$stud == "Castro-NallarE_2015.pathabundance_relab.oralcavity",
            c("pathway", "adonis_p", "adonis_FDR", "kruskal_p", "kruskal_FDR")]
rownames(sig) = sig$pathway


## Load pathway abundances
d_path = curatedMetagenomicData("Castro-NallarE_2015.pathabundance_relab.oralcavity", dryrun = F)[[1]]
meta_path = pData(d_path)
pathways = exprs(d_path)

pathways_cleaned = clean_pathways(pathways)
assigned = pathways_cleaned[[1]]
totals = pathways_cleaned[[2]]

## load taxa abundances
d_tax = curatedMetagenomicData("Castro-NallarE_2015.metaphlan_bugs_list.oralcavity", dryrun = F)[[1]]
meta_tax = pData(d_tax)
abundance = exprs(d_tax)

abundance_s = abundance[ grepl( "s__", rownames(abundance) ) & !(grepl("t__", rownames(abundance))),]
rownames(abundance_s) = sapply(rownames(abundance_s), function(r){
                               split = strsplit(r, "\\|")[[1]]
                               return(paste(split[6], split[7], sep = "."))} ) %>%
                        unname()
abundance_filt = abundance_s[ rownames(abundance_s) %in% assigned$taxa ,]

abundance_bray = as.matrix(vegdist( t(abundance_s), "bray"))

## -------------- individual pathways --------------------------

species_per = assigned[, c("PWY", "taxa")] %>% unique() %>%
                select("PWY") %>% table() %>% sort()
sig$species_per = species_per[rownames(sig)]

prevalence = assigned %>% filter(abundance > 5e-6) %>%
                group_by(PWY, sample) %>% summarize(n = n()) %>% ungroup() %>%
                select(PWY, sample) %>%
                group_by(PWY) %>% summarize(n = n())

median_abundance = assigned %>% group_by(PWY, sample) %>%
                    summarize(ab = sum(abundance)) %>%
                    ungroup() %>% group_by(PWY) %>% summarize(median_abundnace = median(ab)) %>%
                    ungroup() %>% as.data.frame()
left_join(sig, median_abundance, by = c("pathway" = "PWY")) %>% arrange(adonis_p) %>% head(40)
left_join(sig, prevalence, by = c("pathway" = "PWY")) %>% 
            filter(n > 20 & adonis_FDR < 0.05) %>% arrange(adonis_p) %>% head(40)

pwy = "SO4ASSIM-PWY"
assigned_pwy = assigned[assigned$PWY == pwy,]
assigned_pwy$condition = meta_path[assigned_pwy$sample, "study_condition"]

pdf("figures/pathway_abundance.pdf", height = 8, width = 8)
filter(assigned_pwy, abundance > 0) %>% 
ggplot( aes(x = sample, y = abundance, fill = taxa)) + 
    theme_classic() +
    geom_bar(position = "stack", stat = "identity") + 
    facet_wrap(~condition) + 
    ggtitle(pwy) + 
    theme(legend.position = "bottom") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    guides(fill=guide_legend(ncol=2)) + 
    scale_fill_manual(values = pick_colors(length(unique(assigned_pwy$taxa))))

filter(totals, PWY == pwy) %>%
    left_join(meta_path, by = c("sample" = "subjectID")) %>%
    ggplot(aes(x = study_condition, y = abundance)) +
        theme_classic() + 
        geom_boxplot(lwd = 2) + 
        theme(text = element_text(size = 20))

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

df = data.frame(ab = abundance_s["g__Delftia.s__Delftia_acidovorans", ], 
                cond = meta_tax[colnames(abundance_s), "study_condition"])
wilcox.test(ab ~ cond, data = df)

df = data.frame(ab = abundance_s["g__Neisseria.s__Neisseria_subflava", ], 
                cond = meta_tax[colnames(abundance_s), "study_condition"])
wilcox.test(ab ~ cond, data = df)

