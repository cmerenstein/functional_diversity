library(tidyverse)
library(curatedMetagenomicData)
library(SummarizedExperiment)

d = curatedMetagenomicData("YeZ_2018.pathabundance_relab.stool", dryrun = F)[[1]] 
meta = pData(d)
pathways = exprs(d[[1]])

pathways_df =data.frame( pathway = rownames(pathways), pathways)
long = gather(pathways_df, "sample", "abundance", -pathway)

# everything assigned to a taxa
assigned = long[ !(grepl("UNINTEGRATED", long$pathway)) & long$pathway!= "UNMAPPED",]
assigned = separate(assigned, pathway, c("PWY", "description"), sep = ": ") %>%
                separate(description, c("descirption", "taxa"), sep = "\\|", fill = "right")


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






