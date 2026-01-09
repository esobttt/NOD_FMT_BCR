# [FMT Reanalysis]
## Yerim HEO, OMICS Team, Hecto Healthcare
## Last Update: Aug, 28, 2024

# Importing libraries
library(dplyr)
library(ggpicrust2)
library(rstatix)
library(KEGGREST)
library(readxl)
library(writexl)
library(insight)
library(reshape2)
library(pdftools)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(patchwork)
library(stringr)
library(grid)

source("/data1/hyr/02.script/01.Function/00.Taxon_converter.R")
source("/data1/hyr/02.script/Taxon_color_palette.R")
source("/home/hyr/03.Rfunction_script/Determine_statMethods_16s_ver3.R")

pwd <- "/data1/hyr/01.Analysis/01_Genobiome_Service/32.FMT_Re/"
paired <- FALSE
timepoint_levels <- c("Pre-FMT","Post-FMT")
group_levels <- c("NOD-FMT","PBS-T")
if(length(timepoint_levels) >= 2 & length(group_levels)>=2){
  groupxtimepoint_levels <- expand.grid(group = group_levels, timepoint = timepoint_levels) %>% as.data.frame() %>% arrange(group, timepoint) %>% mutate(groupxtimepoint = paste0(group,"_",timepoint)) %>% select(groupxtimepoint) %>% pull()
}
  
analysis_path <- paste0(pwd, "03_Analysis/")
figure_path <- paste0(pwd, "04_Figures/")
setwd(pwd)

timepoint_COLORS <- c("#EBEBEB","#333333")
group_COLORS <- c("#B2182B","#2166AC")



## 01. Data Preprocessing
# metadata
metadata <- read.table(paste0(analysis_path,"metadata.tsv"), sep="\t", header = T) #metadata.tsv

# figure_setting
{
  TEXT_SIZE=c(10)
  TITLE_SIZE=c(11.5)
  ALPHA=c(0.8)
  POINT_SIZE=c(4)
  FAMILY = "ArialMT"
  legend_title <- ""
  set.seed(12345)
}

today <- format(Sys.Date(), "%y%m%d")

if(TRUE %in% grepl("id",colnames(metadata))){
  metadata <- metadata %>% dplyr::rename(fastq=id)
} 
metadata$fastq <- factor(metadata$fastq, levels = unique(metadata$fastq))

#alpha diversity
alpha <- read.delim(paste0(analysis_path, "qiime2_alpha.diversity_rarefied.summary.tsv"))
alpha <- merge(alpha, metadata, by.x = "sample", by.y = "fastq")
alpha <- alpha %>% rename(fastq=sample) %>% select(colnames(metadata), Shannon) %>% arrange(fastq)
alpha <- alpha %>% filter(!is.na(group) & !is.na(timepoint))
alpha %>% select(timepoint, group) %>% table()
alpha %>% select(timepoint, sampleID) %>% table()

#pcoa
pcoa_dirs <- list.files(path = analysis_path, pattern = "qiime2_PCOA_all_asv.xlsx", recursive = TRUE, full.names = TRUE)
pcoa_var <- str_replace_all(pcoa_dirs,paste0(analysis_path,"/pcoa_adonis_results_"),"") %>% str_replace_all(.,"/qiime2_PCOA_all_asv.xlsx","") %>% str_replace_all(.,"/","_")

for (dir in 1:length(pcoa_dirs)){
  tryCatch({
    name_sheets <- excel_sheets(pcoa_dirs[dir])
    pc <- as.data.frame(read_excel(pcoa_dirs[dir], sheet = 1))
    rownames(pc) <- pc$comparison_no
    pc <- pc %>% select(PC1, PC2) 
    pc <- round(pc * 100, digits = 1)
    pc$p.val <- as.data.frame(read_excel(pcoa_dirs[dir], sheet = 2))[,2]
    pc %>% mutate(com=rownames(.)) %>% select(com, PC1, PC2, p.val)
    
    pcoa <- read_excel(pcoa_dirs[dir], sheet = 3)
    pcoa$com <- strsplit(name_sheets[3],".AXIES")[[1]]
    pcoa <- pcoa %>% select(com, `Fastq ID`, PC1, PC2 ) %>% dplyr::rename(sample=`Fastq ID`) %>% arrange(sample, com) %>% mutate(group=metadata$group[match(sample, metadata$fastq)] , timepoint=metadata$timepoint[match(sample, metadata$fastq)])
    if (length(name_sheets) > 3) {
      for (i in 4:length(name_sheets)) {
        temp <- read_excel(pcoa_dirs[dir], sheet = i)
        temp$com <- strsplit(name_sheets[i],".AXIES")[[1]] 
        temp <- temp %>% dplyr::rename(sample=`Fastq ID`) %>% select(com, sample, PC1, PC2 ) %>% arrange(sample, com) %>% mutate(group=metadata$group[match(sample, metadata$fastq)] , timepoint=metadata$timepoint[match(sample, metadata$fastq)])
        pcoa <- rbind(pcoa, temp)
      }
    }
    pcoa <- pcoa %>% mutate(group=metadata$group[match(pcoa$sample, metadata$fastq)] , timepoint=metadata$timepoint[match(pcoa$sample, metadata$fastq)])
    pcoa <- pcoa %>% filter(!is.na(group) & !is.na(timepoint)) %>% rename(PCoA1="PC1", PCoA2="PC2")
    assign(pcoa_var[dir],list("pc"=pc, "pcoa"=pcoa))
  }, error = function(e){
    cat("Error reading PCoA data:", conditionMessage(e), "\n")
    pc <- NULL
  })
}


# Raw abundance
# Get the full Excel sheet names in Raw_abundance_table.xlsx
name_sheets <- excel_sheets(paste0(analysis_path,"Raw_abundance_table.xlsx"))
for (i in 1:length(name_sheets)) {
  assign(strsplit(name_sheets[i],"_rf.tsv")[[1]], t(read_excel(paste0(analysis_path,"Raw_abundance_table.xlsx"), sheet=name_sheets[i], skip = 1)))
}

{
  colnames(phylum) <- sapply(phylum[1,], function(x){
    Taxon_converter(x,2)
  })
  colnames(genus) <- sapply(genus[1,], function(x){
    Taxon_converter(x,6)
  })
}

##phylum
phylum <- as.data.frame.matrix(phylum[-1,])
phylum <- phylum %>% filter(rownames(.) %in% metadata$fastq) 

phylum_df <- phylum %>% dplyr::mutate(fastq=rownames(.))%>% reshape2::melt(., id.var = "fastq", variable.name = c("taxon")) %>% merge(.,metadata, by="fastq")
phylum <- phylum %>% mutate(fastq = rownames(.)) %>% merge(., metadata, by="fastq")


##genus
phylum_genus_table <- data.frame("taxon" = as.vector(genus[1,]),"phylum" = sapply(genus[1,], function(x){
  Taxon_converter(x,2)
}),"family" = sapply(genus[1,], function(x){
  Taxon_converter(x,5)
}) ,
"genus" = sapply(genus[1,], function(x){
  Taxon_converter(x,6)
}))

genus <- as.data.frame.matrix(genus[-1,])
genus <- genus %>% filter(rownames(.) %in% metadata$fastq)

genus_df <- genus %>% dplyr::mutate(fastq=rownames(.))%>% reshape2::melt(., id.var = "fastq", variable.name = c("taxon")) %>% merge(.,metadata, by="fastq")
genus <- genus %>% mutate(fastq = rownames(.)) %>% merge(., metadata, by="fastq")


## 02. Statistatical Test
alpha_stat <- det_statMethods(input_table = alpha, final_group_levels = group_levels, final_timepoint_levels = timepoint_levels, metadata = metadata)
phylum_stat <- det_statMethods(input_table = phylum, final_group_levels = group_levels, final_timepoint_levels = timepoint_levels, metadata = metadata)
phylum_stat$group_avg <- phylum_stat$group_avg %>% mutate(count = ifelse(phylum_stat$group_avg %>% select(group_levels) > 0.01, 1, 0) %>% rowSums() , variable = factor(variable, levels = as.vector(variable)))
phylum_stat$timepoint_avg <- phylum_stat$timepoint_avg %>% mutate(count = ifelse(phylum_stat$timepoint_avg %>% select(timepoint_levels) > 0.01, 1, 0) %>% rowSums(), variable = factor(variable, levels = as.vector(variable)))
phylum_stat$groupxtimepoint_avg <- phylum_stat$groupxtimepoint_avg %>% mutate(count = ifelse(phylum_stat$groupxtimepoint_avg %>% select(groupxtimepoint_levels) > 0.01, 1, 0) %>% rowSums(), variable = factor(variable, levels = as.vector(variable)))
genus_stat <- det_statMethods(input_table = genus, final_group_levels = group_levels, final_timepoint_levels = timepoint_levels, metadata = metadata)
genus_stat$group_avg <- genus_stat$group_avg %>% mutate(count = ifelse(genus_stat$group_avg %>% select(group_levels) > 0.01, 1, 0) %>% rowSums() , variable = factor(variable, levels = as.vector(variable)))
genus_stat$timepoint_avg <- genus_stat$timepoint_avg %>% mutate(count = ifelse(genus_stat$timepoint_avg %>% select(timepoint_levels) > 0.01, 1, 0) %>% rowSums(), variable = factor(variable, levels = as.vector(variable)))
genus_stat$groupxtimepoint_avg <- genus_stat$groupxtimepoint_avg %>% mutate(count = ifelse(genus_stat$groupxtimepoint_avg %>% select(groupxtimepoint_levels) > 0.01, 1, 0) %>% rowSums(), variable = factor(variable, levels = as.vector(variable)))

write_xlsx(Filter(Negate(is.null), alpha_stat), paste0(analysis_path, "alpha_stat_all.xlsx"))
write_xlsx(Filter(Negate(is.null), phylum_stat), paste0(analysis_path, "phylum_stat_all.xlsx"))
write_xlsx(Filter(Negate(is.null), genus_stat), paste0(analysis_path, "genus_stat_all.xlsx"))


phylum_levels <- as.vector(phylum_stat$group_avg$variable)
phylum_df <- phylum_df %>% mutate(taxon=factor(taxon, levels = phylum_levels), group = factor(group, levels = group_levels), timepoint = factor(timepoint, levels =timepoint_levels)) %>% arrange(taxon)

phylum_avg <- phylum_stat$groupxtimepoint_avg

genus_levels <- as.vector(genus_stat$group_avg$variable)
genus_df <- genus_df %>% mutate(taxon=factor(taxon, levels = genus_levels), group = factor(group, levels = group_levels), timepoint = factor(timepoint, levels =timepoint_levels)) %>% arrange(taxon)

genus_avg <- genus_stat$groupxtimepoint_avg

#raw
phylum_Others <- phylum_avg %>% filter(count == 0) %>% select(variable) %>% pull()
phylum_nonOthers <- phylum_avg %>% filter(count > 0) %>% select(variable) %>% pull()
genus_Others <- genus_avg %>% filter(count == 0) %>% select(variable) %>% pull()
genus_nonOthers <- genus_avg %>% filter(count > 0) %>% select(variable) %>% pull()

## 03. Data Visualization
# alpha diversity
g_data<- alpha %>% reshape2::melt(alpha, id.vars=colnames(metadata),
                                  measure.vars=c("Shannon"), 
                                  variable.name="variable") %>% mutate(group = factor(group, levels = group_levels), timepoint = factor(timepoint, levels = timepoint_levels))

p <- ggplot(g_data, aes(x=timepoint, y = value)) + geom_boxplot(aes(fill = group, col = group), alpha = ALPHA, size = 0.3, outlier.shape = NA) +
  xlab("") + ylab("Shannon Index")+ geom_jitter(aes(color = group), position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.3)) +
  theme_bw() + 
  theme(axis.title.x = element_text(size=TITLE_SIZE, vjust=-1),
        axis.title.y = element_text(size=TITLE_SIZE,vjust=2),
        axis.text.x=element_text(size=TEXT_SIZE, angle=0, color='black',hjust=0.5,vjust=0.5),
        axis.text.y=element_text(size=TEXT_SIZE, angle=0, color='black'),
        panel.border = element_rect(colour="black", fill=NA),
        panel.background=element_blank(),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        panel.spacing = unit(1, "lines"),
        legend.text = element_text(color = "black", size=TEXT_SIZE),
        legend.position="right",
        legend.title = element_blank(),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm",),
        strip.background = element_blank(),
        strip.text.x = element_text(size =TITLE_SIZE, face = "bold"),
        strip.text.y = element_text(size =TITLE_SIZE, face = "bold")) +
  scale_fill_manual(values=c(group_COLORS)) +
  scale_color_manual(values=c(group_COLORS)) +
  stat_pvalue_manual(alpha_stat$stat_group %>% filter(variable == "Shannon")%>% filter(p<0.05) %>% mutate(variable2=variable) , label = "{p.signif}",
                     step.group.by="variable", tip.length = 0,
                     step.increase = 0.05) 
p

print_color("Saving 01.alpha_diversity.pdf\n", "white")
pdf(paste0(figure_path,"01.alpha_diversity.pdf"), width = 8.5, height = 5.5, family = FAMILY)
print(p)
dev.off()
print_color("Saved 01.alpha_diversity.pdf\n", "violet")

# beta diversity
# braycurtis
for (var in pcoa_var){
  PC_Val <- get(var)[[1]] %>% mutate(com = rownames(.))
  data <- get(var)[[2]]
  data$group <- factor(data$group, levels = group_levels)
  data$timepoint <- factor(data$timepoint, levels = timepoint_levels)
  table(data$group)
  
  for (i in unique(data$com)){
    data_sel <- data %>% filter(com == i)
    pval <- PC_Val %>% filter(com == i) %>% select(p.val) %>% pull()
    
    pcoa_p <- ggplot(data = data_sel, mapping = aes(x = PCoA1, y = PCoA2, color = group, shape = timepoint)) +
      geom_point(data = data_sel, mapping = aes(), size = POINT_SIZE, alpha = ALPHA) + 
      theme_classic() +
      theme(
        panel.border = element_rect(colour = "black", fill = NA),
        axis.line = element_line(),
        axis.text.x = element_text(size = TEXT_SIZE, angle = 0, color = 'black'),
        axis.text.y = element_text(size = TEXT_SIZE, angle = 0, color = 'black'),
        axis.title.x = element_text(size = TITLE_SIZE, vjust = -5),
        axis.title.y = element_text(size = TITLE_SIZE, vjust = 5),
        legend.key = element_rect(colour = NA, fill = NA, size = unit(0.1, 'cm')),
        legend.text = element_text(size = TEXT_SIZE),
        legend.position = "right",
        plot.margin = margin(0.5, 0.5, 0.8, 0.8, "cm")
      ) +
      xlab(paste0("PC1(",PC_Val[i,]$PC1,"%)")) +
      ylab(paste0("PC2(",PC_Val[i,]$PC2,"%)")) + 
      scale_color_manual(legend_title, values = group_COLORS) +
      scale_shape_manual(legend_title, values = c(16,17)) + 
      stat_ellipse(aes(color = group), type = 't', size = 0.5, level = 0.95, show.legend = F) + 
      annotate("text", label = paste0("p=", round(pval,4)), color = "red",  x = max(data_sel$PCoA1) +0.01, y = max(data_sel$PCoA1) + 0.15)
    
    print_color(paste0("Saving 02.", var, "_", i,".pdf\n"), "white")
    pdf(paste0(figure_path,"02.", var, "_", i,".pdf"),  width = 6.8, height = 5, family = FAMILY)
    print(pcoa_p)
    dev.off()
    print_color(paste0("Saved 02.", var, "_",i,".pdf\n"), "violet")
  }
}


# Taxon abundance
blank <- ggplot()+
  theme_minimal()

# Bar plot
Others_df <- phylum_df %>% filter(taxon %in% phylum_Others) %>% group_by(fastq) %>% dplyr::summarize(value=sum(as.numeric(value))) %>% mutate(taxon="Others") %>% merge(., metadata, by = "fastq") %>% select(taxon, colnames(metadata), value)
Bar_df <- rbind(phylum_df %>% filter(taxon %in% phylum_nonOthers), Others_df)
average <- Bar_df %>% group_by(timepoint, taxon) %>% dplyr::summarize(mean = mean(as.numeric(value))) %>% mutate(group=paste0(timepoint, "\nAVG"))


TAXO <- factor(average %>% select(taxon) %>% pull() %>% unique() ,
               levels = average %>% select(taxon) %>% pull() %>% unique())

TAXO_COLOR <- Phylum_COLOR[match(TAXO, rownames(Phylum_COLOR)),]
Phylum_order <- factor(levels(TAXO), levels = c(phylum_levels,"Others"))

Bar_df$fastq <- factor(Bar_df$fastq, levels = metadata$fastq)
bar_order <- Bar_df$fastq
levels(bar_order)


# genus (others by phylum)
phylum_order <- c((rbind(as.vector(phylum_nonOthers[order(phylum_nonOthers)]), paste0(phylum_nonOthers[order(phylum_nonOthers)], "_Others")) %>% as.vector()), "Others")

Bar_df <- genus_df %>% mutate(phylum = phylum_genus_table$phylum[match(taxon, phylum_genus_table$genus)]) %>% mutate(new_phylum = ifelse((taxon %in% genus_nonOthers), as.vector(phylum), ifelse((taxon %in% genus_Others)&(phylum %in% phylum_Others), "Others", ifelse((taxon %in% genus_Others)&(phylum %in% phylum_nonOthers),paste0(phylum,"_Others"),""))), new_taxon = ifelse((taxon %in% genus_nonOthers), as.vector(taxon),ifelse((taxon %in% genus_Others)&(phylum %in% phylum_Others), "Others",paste0(phylum,"_Others")))) %>% mutate(new_phylum = factor(new_phylum, levels = c(phylum_order))) %>% arrange(new_phylum) %>% mutate(taxon=new_taxon) 

genus_order <- Bar_df %>% mutate(new_phylum = factor(new_phylum, levels = phylum_order)) %>% arrange(new_phylum) %>% select(taxon) %>% unique() %>% pull()
selected_phylum <- Bar_df %>% mutate(new_phylum = factor(new_phylum, levels = phylum_order)) %>% arrange(new_phylum) %>% select(new_phylum) %>% unique() %>% pull() %>% str_replace_all(.,"_Others","") %>% unique()
average <- Bar_df %>% group_by(timepoint, fastq ,taxon) %>% dplyr::summarize(value = sum(as.numeric(value))) %>% group_by(timepoint, taxon) %>% summarize(mean = mean(as.numeric(value))) %>% mutate(group=paste0(timepoint, "\nAVG"))

TAXO <- factor(unique(Bar_df$taxon), 
               levels = c(genus_order)) 

COLOR_list <- data.frame("genus" = TAXO, 
                         "phylum" = phylum_genus_table$phylum[match(TAXO,phylum_genus_table$genus)])
COLOR_list[is.na(COLOR_list$phylum),"phylum"] <- str_replace_all(COLOR_list[is.na(COLOR_list$phylum),"genus"],"_Others","")

COLOR_list <- COLOR_list %>% mutate(genus = factor(genus, levels = c(genus_order))) %>% arrange(genus)

normalized_count <- table(COLOR_list$phylum)/max(table(COLOR_list$phylum))
legends <- rev(paste0(unique(COLOR_list$phylum), "_legend"))
nrow <- length(unique(COLOR_list$phylum))
rel_heights <- rev(normalized_count[unique(COLOR_list$phylum)] + c(0,rep(0.2, length(unique(COLOR_list$phylum))-1)))

Bar_df$sample <- factor(Bar_df$sample, levels = metadata$fastq)

color <- c()
for (i in unique(COLOR_list$phylum)){
  color <- c(color,colorRampPalette(c(Phylum_COLOR[i,],hex_lighter(Phylum_COLOR[i,])))(nrow(COLOR_list[COLOR_list$phylum == i,])))
}
COLOR_list$color <- color

rownames(COLOR_list) <- COLOR_list$genus

COLOR_list <- t(COLOR_list) 
TAXO_COLOR <- COLOR_list["color",]

for (i in selected_phylum){
  select <- Bar_df %>% filter(new_phylum == i | new_phylum == paste0(i,"_Others")) 
  barp <-  ggplot(select, aes(x=sample, y = (as.numeric(value) * 100),
                              fill=factor(taxon,levels=rev(unique(genus_order)))),
                  group=group)+
    geom_bar(position=position_stack(reverse = FALSE)
             ,stat="identity"
             ,width=1) +
    
    theme(
      legend.title = element_text(family = FAMILY, face = "bold", size=TITLE_SIZE),
      legend.text = element_text(size=TEXT_SIZE),
      legend.key.height = unit(1, "lines"),  
      legend.key.width = unit(1, "lines"),  
      legend.position = "right")+
    scale_fill_manual(name=i,
                      values = TAXO_COLOR)+
    scale_y_continuous(expand = expansion(mult = c(0, 0.01))) +
    guides(fill = guide_legend(ncol = 1, keyheight = 1.1 ))
  assign(paste0(i,"_legend"), get_legend(barp))
}

# Normalize species count for relative legend height
complete_legend <- plot_grid(plotlist = lapply(legends, get), ncol=1, nrow = nrow, align = "v", rel_heights = rel_heights)

bar_p <-  ggplot(Bar_df, aes(x=fastq, y =as.numeric(value) * 100,
                             fill=factor(taxon,levels=rev(genus_order)), group=group)) + 
  xlab("")+
  ylab("Relative abundance(%)")+
  geom_bar(position=position_stack(reverse = FALSE)
           ,stat="identity"
           ,width=1) +
  facet_grid(.~ factor(timepoint,levels=timepoint_levels) + factor(group, levels = group_levels),
             scales = 'free_x',
             space = 'free',
             switch='x',
             drop=T) +
  theme_minimal()+
  theme(
    axis.title = element_text(size=TITLE_SIZE,color = "black"),
    axis.text.x=element_blank(),
    axis.text.y=element_text(size=TEXT_SIZE,color = "black"),
    axis.line = element_line(size = 0.1),
    panel.grid.minor=element_blank(), 
    panel.grid.major=element_blank(), 
    panel.background=element_blank(),
    panel.spacing = unit(0.1, "lines"),
    legend.text = element_text(size=TEXT_SIZE),
    legend.position = "none",
    strip.placement = "outside",
    strip.text.x = element_text(size=TEXT_SIZE,color = "black",angle = 0, face = "bold"),
    plot.margin = unit(c(0.5, 0, 0.5, 0.5), "cm",)
  )+
  scale_fill_manual(name="",
                    values = TAXO_COLOR)+
  scale_y_continuous(expand = expansion(mult = c(0, 0.01)))
bar_p

avg_p <- ggplot(average, aes(x=group, y = (mean*100),
                             fill=factor(taxon, levels=rev(genus_order))))+
  xlab("")+
  ylab("")+
  geom_bar(position=position_stack(reverse = FALSE)
           ,stat="identity"
           ,width=1) +
  facet_grid(. ~ factor(group,levels=paste0(timepoint_levels,"\nAVG")),
             scales = 'free_x',
             space = 'free',
             switch='x',
             drop=T
  )+
  theme_minimal()+
  theme(
    axis.title = element_text(size=TITLE_SIZE,color = "black"),
    axis.text.x=element_blank(),
    
    axis.text.y=element_text(size=TEXT_SIZE,color = "black"),
    axis.line = element_line(size = 0.1),
    panel.grid.minor=element_blank(),
    panel.grid.major=element_blank(),
    panel.background=element_blank(),
    panel.spacing = unit(0.1, "lines"),
    legend.title = element_blank(),
    legend.text = element_text(size=TEXT_SIZE),
    legend.position = "right",
    strip.placement = "outside",
    strip.text.x = element_text(size=TEXT_SIZE,color = "black",angle = 0, face = "bold"),
    plot.margin = unit(c(0.5, 0, 0.5, 0.5), "cm"))+  
  scale_fill_manual(name="", values = TAXO_COLOR)+
  scale_y_continuous(expand = expansion(mult = c(0, 0.01))) + 
  guides(fill = guide_legend(ncol = 1 , keyheight = 1.1))

p_legend <- get_legend(avg_p)

pdf(paste0(figure_path,"03_3.Bar_plot_genus_ver2_timepoint.pdf"), width = 16.5, height = 10, family = FAMILY)
print(bar_p+blank+complete_legend+plot_layout(guides = 'collect', ncol = 3, widths = c(8,0.2,2.5)))
dev.off()

genus_box_selected <- c("Bacteroides","Blautia","Escherichia-Shigella","Unclassified_family_Lachnospiraceae","Muribaculaceae","Clostridia_vadinBB60_group")

genus_box <- genus_df %>% filter(taxon %in% genus_box_selected) %>% mutate(taxon = sapply(as.vector(taxon), function(x){
  ifelse((grepl("Unclassified", x) |grepl("Uncultured", x)), genus_replace_to_newline(as.character(x)), x)
}) %>% as.vector()) %>% ggplot(., aes(x=timepoint, y = as.numeric(value) * 100)) + geom_boxplot(aes(fill = group, col = group), alpha = ALPHA, size = 0.3, outlier.shape = NA) +
  facet_wrap(.~ taxon, scales = "free_y",ncol = 2) +
  xlab("") + ylab("Relative abundance(%) ")+geom_jitter(aes(color = group), position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.1)) +
  theme_pubr() +
  theme(axis.title.x = element_text(size=TITLE_SIZE,vjust=-1),
        axis.title.y = element_text(size=TITLE_SIZE,vjust=2),
        axis.text.x=element_text(size=TEXT_SIZE, angle=0, color='black',hjust=0.5,vjust=1),
        axis.text.y=element_text(size=TEXT_SIZE, angle=0, color='black'),
        panel.border = element_rect(colour="black", fill=NA),
        panel.background=element_blank(),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        panel.spacing = unit(1, "lines"),
        legend.text = element_text(color = "black", size=TEXT_SIZE),
        legend.position="right",
        legend.title = element_blank(),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm",),
        strip.background = element_blank(),
        strip.text.x = element_text(size =TITLE_SIZE-1.5, face = "bold"),
        strip.text.y = element_text(size =TITLE_SIZE-1.5, face = "bold")) +
  scale_fill_manual(values=c(group_COLORS)) +
  scale_color_manual(values=c(group_COLORS)) + 
  stat_pvalue_manual(genus_stat$stat_group %>% filter(variable %in% genus_box_selected) %>% 
                       filter(p < 0.05) %>% mutate(taxon= sapply(as.vector(variable), function(x){
                         ifelse((grepl("Unclassified", x) |grepl("Uncultured", x)), genus_replace_to_newline(as.character(x)), x)
                       }) %>% as.vector())  %>% mutate(y.position=y.position * 100), label = "{p.signif}",
                     step.group.by="taxon", tip.length = 0, step.increase = 0.07) + 
  scale_y_continuous(expand = expansion(mult = 0.099))
genus_box

height <- ifelse(length(genus_box_selected) %% 2 > 0, ((length(genus_box_selected) +1) %/% 2) * 2.5, length(genus_box_selected)%/%2*2.5)
width <- ifelse(length(genus_box_selected)/2 < 1,   length(genus_box_selected) %% 2 * 4, 8)

pdf(paste0(figure_path,"04_2.Box_plot_signif_genus_filtered_ncol2.pdf"), width = width , height = height, family = FAMILY)
print(genus_box)
dev.off()




#ANCOM
ancom <- read.table(paste0(analysis_path, "differentials-ancombc2/ancom_NODvsPBS.tsv"),header = T)
ancom$LFC <- -ancom$LFC
ancom$id <- sapply(ancom$id, function(x){
  Taxon_converter(x,6)
})

ancom$Class <- ifelse(ancom$LFC > 0, "NOD-FMT", "PBS-T")
ancom <- ancom %>% filter(pvalue < 0.05) %>% arrange(desc(LFC))
ancom <- ancom %>% filter(abs(LFC) > 2)
label_y = ifelse(ancom$LFC < 0, 0.2,-0.2)
label_hjust = ifelse(ancom$LFC > 0, 1, 0)

p <- ggplot(ancom, aes(x=rev(id), 
                       y=as.numeric(LFC), fill = Class)) +
  geom_bar(stat = "identity", linewidth=0.1, width = 0.8,alpha=0.8,
           position=position_dodge2(preserve = "single")
           ,color="#4B4A4A"
  ) + 
  geom_text(aes(y = label_y, label = id, hjust = label_hjust),size=4) +
  coord_flip() +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(colour="black",size=TEXT_SIZE),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size=TITLE_SIZE),
        axis.title.y = element_blank(),
        legend.position = "right",
        legend.justification=0.9, 
        legend.box.spacing = unit(-1, "cm"),
        legend.text = element_text(size=TEXT_SIZE),
        legend.spacing.x = unit(0.3, 'cm'),
        legend.spacing.y = unit(0.3, 'cm'),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(colour = "grey80", linetype = "dashed"),
        panel.grid.minor.x = element_blank(),
        panel.spacing = unit(1, "lines"),
        legend.key.size = unit(0.3, "cm"),
        strip.text.x = element_text(size=TEXT_SIZE)) +
  scale_x_discrete(limits = ancom$id) + 
  scale_y_continuous("Log Fold Change(LFC)",
                     breaks = -7:7,
                     limits = c(-7,7)
  ) + guides(fill = guide_legend(byrow = TRUE))
legend_title <- ""

ancom_plot <- p + scale_fill_manual(legend_title,
                                    values=group_COLORS)
ancom_plot


pdf(paste0(figure_path, "05.ancom_p0.05_LFC2_new.pdf"), width = 12, height = length(ancom$id)/4.5, family = "ArialMT")
print(ancom_plot)
dev.off()

write.table(ancom, paste0(analysis_path, "differentials-ancombc2/ancom_NODvsPBS_renamed.tsv"), quote = F, sep = "\t", row.names = F)

##### Convert all figures to TIFF
{
  setwd(figure_path)
  lapply(list.files(path = figure_path, pattern = "\\.pdf$", full.names = TRUE)[file.info(list.files(path = figure_path, pattern = "\\.pdf$", full.names = TRUE))$size > 0], function(x){
    pdf_convert(x, format = "tiff", dpi = 600)
  })
  
  dir.create(paste0(figure_path,"tiff"))
  system(paste0("mv ",figure_path,"*.tiff ",figure_path,"tiff"))
}



save.image(paste0(figure_path, "16s_Analysis_Results_",format(today , format="%y%m%d"),".RData"))

