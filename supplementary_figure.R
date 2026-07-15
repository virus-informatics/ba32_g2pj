library(tidyverse)
library(stats)
library(base)
library(circlize)
library(ComplexHeatmap)
library(data.table)
library(datasets)
library(grDevices)
library(patchwork)
library(RColorBrewer)
library(cmdstanr)

args = commandArgs(trailingOnly=T)

set_cmdstan_path("cmdstan-2.38.0")
cmdstan_path()

########## args ##########
#Change when using new input
out_prefix <- "2026_07_05"
download_date <- gsub("_", "-", out_prefix)
date_w_space <- gsub("_", "", out_prefix)

metadata.name <- paste(out_prefix,"/metadata_tsv_",out_prefix,"/metadata.tsv",sep = "")
mut.info.name <- paste(out_prefix,"/metadata_tsv_",out_prefix,"/metadata.mut_long.tsv",sep = "")
nextclade.name <- paste(out_prefix,"/nextclade.tsv",sep = "") 
stan_f.name <- 'multinomial_independent.ver2.stan'
country.name <- 'country_info.txt'

#output
mutation_figure.name <- paste(date_w_space,'.mut_comparision.NB181.pdf', sep = "")

dir <- paste(out_prefix,sep = "")
setwd(dir)

## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##########parameters##########
##general
core.num <- 4
variant.ref <-"XEC"

##period to be analyzed
date.start <- as.Date("2025-08-01")
date.end <- as.Date("2026-07-01")

##min numbers
limit.count.analyzed <- 30

##Transmissibility
bin.size <- 1
generation_time <- 2.1

##model
multi_nomial_model <- cmdstan_model(stan_f.name)

## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
########## Filtering metadata #########

########## Added script to use nextclade 

nextclade <- fread(nextclade.name,header=T,sep="\t",quote="",check.names=T)
metadata <- fread(metadata.name,header=T,sep="\t",quote="",check.names=T)

country.info <- read.table(country.name,header=T,sep="\t",quote="")
country.info <- country.info %>% select(-region)

temp <- nextclade %>% left_join(metadata,by=c('seqName'='Virus.name'),relationship = "many-to-many") # Merge the nextclade and Gisaid data
temp <- temp %>% rename("Virus.name" = "seqName") # %>% select(-Pango_n_clade) %>% rename("Pango_n_clade"="Nextclade_pango")
metadata <- temp


##########

## metadata

# i) only ‘original passage’ sequences
# iii) host labelled as ‘Human’
# iv) sequence length above 28,000 base pairs
# v) proportion of ambiguous bases below 2%.


metadata.filtered <- metadata %>%
  distinct(Accession.ID,.keep_all=T) %>%
  filter(Host == "Human",
         #!N.Content > 0.02 | is.na(N.Content),
         str_length(Collection.date) == 10,
         Sequence.length > 28000,
         Passage.details.history == "Original",
         Pango.lineage != "",
         Pango.lineage != "None",
         Pango.lineage != "Unassigned",
         !str_detect(Additional.location.information,"[Qq]uarantine")
  )

metadata.filtered <- metadata.filtered %>%
  mutate(Collection.date = as.Date(Collection.date),
         region = str_split(Location," / ",simplify = T)[,1],
         country = str_split(Location," / ",simplify = T)[,2],
         state = str_split(Location," / ",simplify = T)[,3])

metadata.filtered <- merge(metadata.filtered,country.info,by="country")
metadata.filtered <- metadata.filtered[!duplicated(metadata.filtered$Virus.name),]

nrow(metadata.filtered)

##########
data_add.name <- paste(out_prefix,"/add/add.metadata.filtered.tsv",sep = "")

data_add <- fread(data_add.name,header=T,sep="\t",quote="",check.names=T)
data_add <- data_add %>%
  mutate(
    region = ifelse(!is.na(region.y), region.y, region.x)
  )
data_add <- data_add %>%
  select(-region.x, -region.y)
data_add <- data_add %>% select(Accession.ID, Pango.lineage, region, Collection.date, clade, country, partiallyAliased)

data_prev <- metadata.filtered
nrow(data_prev)
data_prev <- data_prev %>% select(Accession.ID, Pango.lineage, region, Collection.date, clade, country, partiallyAliased)
data_prev <- data_prev %>% filter(Collection.date >= date.start, Collection.date <= date.end)

data_combind <- rbind(data_prev, data_add) %>% distinct(Accession.ID,.keep_all=T)
metadata.filtered <- data_combind
metadata.filtered <- metadata.filtered %>% filter(Collection.date >= date.start, Collection.date <= date.end)
nrow(data_combind)

min(metadata.filtered$Collection.date)
max(metadata.filtered$Collection.date)



## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##use this only when clade is needed 
target_clades <- c("24A", "24C", "24E", "24F", "25A", "25B", "25C", "25E", "25I", "25J")
clade.interest.v <- c("24A (JN.1)","24C (KP.3)","24E (KP.3.1.1)","24F (XEC)","25A (LP.8.1)","NB.1.8.1")

metadata.filtered.clade <- metadata.filtered %>%
  mutate(Pango_n_clade = if_else(clade %in% target_clades,clade,Pango.lineage))

metadata.filtered.clade <- metadata.filtered.clade %>%
  mutate(Pango_n_clade = recode(Pango_n_clade,
                                "24A" = "24A (JN.1)",
                                "24C" = "24C (KP.3)",
                                "24E" = "24E (KP.3.1.1)",
                                "24F" = "24F (XEC)",
                                "25A" = "25A (LP.8.1)",
                                "25B" = "25B (NB.1.8.1)",
                                "25C" = "25C (XFG)",
                                "25E" = "25E (XFJ)",
                                "25I" = "25I (BA.3.2)",
                                "25J" = "25J (PY.1)"
  ))

metadata.filtered.clade <- metadata.filtered.clade %>%
  mutate(Pango_n_clade = if_else(startsWith(partiallyAliased, "BA.3.2.1"),"RE.1",
                                 if_else(startsWith(partiallyAliased, "BA.3.2.2"),"RE.2", 
                                         if_else(partiallyAliased == "BA.3.2", "BA.3.2", Pango_n_clade))))

# backup <- metadata.filtered
metadata.filtered <- metadata.filtered.clade


metadata.filtered %>% filter(startsWith(partiallyAliased, "BA.3.2.2")) %>% nrow()
metadata.filtered %>% filter(startsWith(partiallyAliased, "BA.3.2.2"), country == "United Kingdom") %>% nrow()
metadata.filtered %>% filter(startsWith(partiallyAliased, "BA.3.2.2"), country == "Germany") %>% nrow()
metadata.filtered.clade %>% group_by(Pango.lineage) %>% summarize(count = n()) %>% arrange(desc(count)) %>% head(count.pango.df,n=20)

## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
metadata.filtered.clade %>% group_by(Pango_n_clade) %>% summarize(count = n()) %>% arrange(desc(count)) %>% head(count.pango.df,n=20)
a <- metadata.filtered %>% group_by(Pango_n_clade,country) %>% 
  summarize(count = n())  %>% 
  arrange(desc(count))
a %>% 
  filter(startsWith(Pango_n_clade, "RE.2")) %>% 
  head(n = 20)

metadata.filtered %>% group_by(Pango_n_clade, country) %>% summarize(count = n()) %>% arrange(desc(count)) %>% 
  filter(country == "United Kingdom") %>%
  head(count.pango.df,n=20)

metadata.filtered %>% filter(startsWith(Pango_n_clade, "RE.1"))%>% group_by(country)%>% summarize(count = n()) %>% arrange(desc(count)) 
metadata.filtered %>% filter(startsWith(Pango_n_clade, "RE.2"))%>% group_by(country)%>% summarize(count = n()) %>% arrange(desc(count)) 

count.pango.df <- metadata.filtered.clade %>% group_by(Pango_n_clade, country) %>% summarize(count = n()) %>% arrange(desc(count)) 
count.pango.df %>% filter(country=="Australia", count >= limit.count.analyzed) %>% arrange(desc(count)) %>% head(n=15)
count.pango.df %>% filter(country=="Germany", count >= limit.count.analyzed) %>% arrange(desc(count)) %>% head(n=15)
count.pango.df %>% filter(country=="USA", count >= limit.count.analyzed) %>% arrange(desc(count)) %>% head(n=15)
count.pango.df %>% filter(country=="Netherlands", count >= limit.count.analyzed) %>% arrange(desc(count)) %>% head(n=15)
count.pango.df %>% filter(country=="United Kingdom", count >= limit.count.analyzed) %>% arrange(desc(count)) %>% head(n=15)
count.pango.df %>% filter(country=="Spain", count >= limit.count.analyzed) %>% arrange(desc(count)) %>% head(n=15)
count.pango.df %>% filter(country=="Singapore", count >= limit.count.analyzed) %>% arrange(desc(count)) %>% head(n=15)
count.pango.df %>% filter(country=="New Zealand", count >= limit.count.analyzed) %>% arrange(desc(count)) %>% head(n=15)
count.pango.df %>% filter(country=="France", count >= limit.count.analyzed) %>% arrange(desc(count)) %>% head(n=15)
count.pango.df %>% filter(country=="South Korea", count >= limit.count.analyzed) %>% arrange(desc(count)) %>% head(n=15)
count.pango.df %>% filter(country=="France", count >= limit.count.analyzed) %>% arrange(desc(count)) %>% head(n=15)
count.pango.df %>% filter(country=="Japan", count >= limit.count.analyzed) %>% arrange(desc(count)) %>% head(n=15)
count.pango.df %>% filter(country=="Canada", count >= limit.count.analyzed) %>% arrange(desc(count)) %>% head(n=15)
count.pango.df %>% filter(country=="Italy", count >= limit.count.analyzed) %>% arrange(desc(count)) %>% head(n=15)

count.pango.df %>% filter(Pango_n_clade=="RE.2") %>% arrange(desc(count)) %>% head(n=15)
count.pango.df %>% filter(Pango_n_clade=="RE.1") %>% arrange(desc(count)) %>% head(n=15)
count.pango.df %>% filter(Pango_n_clade=="BA.3.2") %>% arrange(desc(count)) %>% head(n=15)


variant.ref <-"25C (XFG)" 
lineage.interest.v <- c("25C (XFG)","25B (NB.1.8.1)","25E (XFJ)",
                        "RE.2","RE.1")


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


########## Predict effective reproductive number (Re) ##########

country.interest <- c("Australia","Germany","Canada","Spain","United Kingdom","Netherlands","Japan")

##min numbers
limit.count.analyzed <- 40

plot.l <- list()

for (a in country.interest) {
#a <- "USA" #"Hong Kong" #Australia
#a <- args[2]

metadata.filtered.interest <- metadata.filtered %>% filter(country == a)
metadata.filtered.interest <- metadata.filtered.interest %>% mutate(date.num = as.numeric(Collection.date) - min(as.numeric(Collection.date))  + 1, date.bin = cut(date.num,seq(0,max(date.num),bin.size)), date.bin.num = as.numeric(date.bin))
metadata.filtered.interest <- metadata.filtered.interest %>% filter(!is.na(date.bin))

list_Re.name <- paste(date_w_space, "_list_Re_",a,".tsv",sep = "")
metadata.filtered.epi_set.list <- paste(metadata.filtered.interest$Accession.ID)
write.table(metadata.filtered.epi_set.list, list_Re.name, col.names=F, row.names=F, sep="\n", quote=F)

##count variants per day
metadata.filtered.interest.bin <- metadata.filtered.interest %>% group_by(date.bin.num, Pango_n_clade) %>% summarize(count = n()) %>% ungroup()

metadata.filtered.interest.bin.spread <- metadata.filtered.interest.bin %>% spread(key=Pango_n_clade,value = count)
metadata.filtered.interest.bin.spread[is.na(metadata.filtered.interest.bin.spread)] <- 0
metadata.filtered.interest.bin.spread <- metadata.filtered.interest.bin.spread

X <- as.matrix(data.frame(X0 = 1, X1 = metadata.filtered.interest.bin.spread$date.bin.num))

Y <- metadata.filtered.interest.bin.spread %>% select(- date.bin.num)
Y <- Y[,c(variant.ref,colnames(Y)[-which(colnames(Y)==variant.ref)])]

count.group <- apply(Y,2,sum)
count.total <- sum(count.group)
prop.group <- count.group / count.total

Y <- Y %>% as.matrix()
apply(Y,2,sum)

group.df <- data.frame(group_Id = 1:ncol(Y), group = colnames(Y))

Y_sum.v <- apply(Y,1,sum)

data.stan <- list(K = ncol(Y),
                  D = 2,
                  N = nrow(Y),
                  X = X,
                  Y = Y,
                  generation_time = generation_time,
                  bin_size = bin.size,
                  Y_sum = c(Y_sum.v))

fit.stan <- multi_nomial_model$sample(
  data=data.stan,
  iter_sampling=4000,
  iter_warmup=1000,
  seed=1234,
  parallel_chains = 4,
  #adapt_delta = 0.99,
  max_treedepth = 15,
  #pars=c('b_raw'),
  chains=4)

#growth rate
stat.info <- fit.stan$summary("growth_rate") %>% as.data.frame()
stat.info$Nextclade_pango <- colnames(Y)[2:ncol(Y)]

stat.info.q <- fit.stan$summary("growth_rate", ~quantile(.x, probs = c(0.025,0.975))) %>% as.data.frame() %>% rename(q2.5 = `2.5%`, q97.5 = `97.5%`)
stat.info <- stat.info %>% inner_join(stat.info.q,by="variable")

out.name <- paste('growth_rate.wo_strata.', a, '.txt', sep='')
write.table(stat.info, out.name, col.names=T, row.names=F, sep="\t", quote=F)


draw.df.growth_rate <- fit.stan$draws("growth_rate", format = "df") %>% as.data.frame() %>% select(! contains('.'))
draw.df.growth_rate.long <- draw.df.growth_rate %>% gather(key = class, value = value)

draw.df.growth_rate.long <- draw.df.growth_rate.long %>% mutate(group_Id = str_match(draw.df.growth_rate.long$class,'growth_rate\\[([0-9]+)\\]')[,2] %>% as.numeric() + 1)
draw.df.growth_rate.long <- merge(draw.df.growth_rate.long,group.df,by="group_Id") %>% select(value,group)
draw.df.growth_rate.long <- draw.df.growth_rate.long %>% group_by(group) %>% filter(value>=quantile(value,0.005),value<=quantile(value,0.995))
draw.df.growth_rate.long <- rbind(data.frame(group=variant.ref,value=1),draw.df.growth_rate.long)
draw.df.growth_rate.long <- draw.df.growth_rate.long %>% filter(group %in% lineage.interest.v)

draw.df.growth_rate.long <- draw.df.growth_rate.long %>% mutate(group = factor(group, levels=lineage.interest.v))

col.v <- brewer.pal(length(lineage.interest.v) + 1, "Set1")[c(1:8)]
g1 <- ggplot(draw.df.growth_rate.long,aes(x=group,y=value,color=group,fill=group))
g1 <- g1 + geom_hline(yintercept=1, linetype="dashed", alpha=0.5)
g1 <- g1 + geom_violin(alpha=0.6,scale="width")
g1 <- g1 + stat_summary(geom="pointrange",fun = median, fun.min = function(x) quantile(x,0.025), fun.max = function(x) quantile(x,0.975), size=0.5,fatten =1.5)
# g1 <- g1 + scale_fill_manual(valus = color_maping)
g1 <- g1 + scale_color_manual(values=col.v)
g1 <- g1 + scale_fill_manual(values=col.v)
g1 <- g1 + theme_classic()
g1 <- g1 + theme(panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(),
                 strip.text = element_text(size=8))
g1 <- g1 + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
g1 <- g1 + ggtitle(a)
g1 <- g1 + xlab('') + ylab('Relative Re')
g1 <- g1 + theme(legend.position = 'none')
g1 <- g1 + scale_y_continuous(limits=c(0.6,1.8),breaks=c(0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8))
#g1 <- g1 + scale_y_continuous(limits=c(0.5,2.0),breaks=c(0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0))
g1

pdf.name <- paste('Re_violin_plot.wo_strata.', a, '.pdf', sep='')
plot.l[[pdf.name]] <- g1
pdf(pdf.name,width=3.5,height=5.0)
plot(g1)
dev.off()

data_Id.df <- data.frame(data_Id = 1:length(X), date_Id = X[,2], Y_sum = Y_sum.v, date = as.Date(X[,2],origin=date.start)-1)


data.freq <- metadata.filtered.interest.bin %>% rename(group = Pango_n_clade) %>% group_by(date.bin.num) %>% mutate(freq = count / sum(count))
data.freq <- data.freq %>% mutate(date = as.Date(date.bin.num,origin=date.start)-1)



draw.df.theta <- fit.stan$draws("theta", format = "df") %>% as.data.frame() %>% select(! contains('.'))


#------- for USA (giant data) only

set.seed(123) 

draw.df.theta.sampled <- draw.df.theta %>%
  slice_sample(prop = 0.5)

rm(draw.df.theta)
gc()

draw.df.theta <- draw.df.theta.sampled

#-------



draw.df.theta.long <- draw.df.theta %>%
  pivot_longer(cols = everything(), names_to = "class", values_to = "value")

draw.df.theta.long <- draw.df.theta.long %>% mutate(data_Id = str_match(class,'theta\\[([0-9]+),[0-9]+\\]')[,2] %>% as.numeric(),
                                                    group_Id = str_match(class,'theta\\[[0-9]+,([0-9]+)\\]')[,2] %>% as.numeric())

draw.df.theta.long <- draw.df.theta.long %>% inner_join(data_Id.df %>% select(data_Id,date), by = "data_Id")

draw.df.theta.long.sum <- draw.df.theta.long %>% group_by(group_Id, date) %>% summarize(mean = mean(value),ymin = quantile(value,0.025),ymax = quantile(value,0.975))
draw.df.theta.long.sum <- draw.df.theta.long.sum %>% inner_join(group.df,by="group_Id")

draw.df.theta.long.sum.filtered <- draw.df.theta.long.sum %>% filter(group %in% lineage.interest.v) %>% mutate(group = factor(group,levels=lineage.interest.v))

g2 <- ggplot(draw.df.theta.long.sum.filtered,aes(x=date, y = mean, fill=group, color = group))
#g2 <- g2 + geom_point(aes(y=freq,size=count),alpha=0.4)
g2 <- g2 + geom_ribbon(aes(ymin=ymin,ymax=ymax), color=NA,alpha=0.2)
g2 <- g2 + geom_line(linewidth=0.3)
g2 <- g2 + scale_x_date(date_labels = "%y-%m", date_breaks = "1 months", date_minor_breaks = "1 month", limits = c(date.start, date.end))
g2 <- g2 + theme_classic()
g2 <- g2 + theme(panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(),
                 strip.text = element_text(size=8)
)
g2 <- g2 + ggtitle(a)
g2 <- g2 + scale_color_manual(values = col.v)
g2 <- g2 + scale_fill_manual(values = col.v)
g2 <- g2 + scale_y_continuous(limits=c(0,1.0),breaks=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0))
#g2 <- g2 + scale_y_continuous(limits=c(0,1.0),breaks=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0))
g2

pdf.name <- paste('lineage_dynamics.wo_strata.', a, '.pdf', sep='')
plot.l[[pdf.name]] <- g2
pdf(pdf.name,width=5,height=5)
plot(g2)
dev.off()
}



if(0){
########## Mutation frequency plot for tree reconstruction #########
### mut.info
### make sure to use python script that includes all needed period to make mut.info

lineage.interest.v <- c("XFG", "NB.1.8.1","XFJ","RE.2","RE.1")

mut.info <- fread(mut.info.name,header=T,sep="\t",quote="",check.names=T)

mut.info.merged <- mut.info %>% inner_join(metadata.filtered %>% select(Id=Accession.ID,Pango.lineage),by="Id")
mut.info.merged <- mut.info.merged %>% as.data.frame()
mut.info.merged <- mut.info.merged %>% filter(Pango.lineage %in% lineage.interest.v) #%>% slice_sample(n = 5000)
mut.info.merged <- mut.info.merged %>% mutate(mut = str_replace(mut, ":", "_")) %>% mutate(prot=str_split(mut, "_", simplify=T)[,1], mut.mod=gsub("[A-Z]", "", str_split(mut, "_", simplify=T)[,2], ignore.case=TRUE))

metadata.filtered.2 <- metadata.filtered %>% filter(Accession.ID %in% as.character(mut.info.merged$Id))
metadata.filtered.2 <- metadata.filtered.2 %>% filter(!(Pango.lineage=="XBB.1" & ((! grepl("Spike_F486S",AA.Substitutions)) | (!grepl("Spike_T478K",AA.Substitutions)))))
metadata.filtered.2 <- metadata.filtered.2 %>% distinct(Accession.ID,.keep_all=T)

list_mut.name <- paste(date_w_space,"_list_mut.tsv",sep = "")
metadata.filtered.epi_mut_set.list <- paste(metadata.filtered.2$Accession.ID)
write.table(metadata.filtered.epi_mut_set.list, list_mut.name, col.names=F, row.names=F, sep="\n", quote=F)

mut.info.merged <- mut.info.merged %>% filter(Id %in% as.character(metadata.filtered.2$Accession.ID))

count.pango.df <- metadata.filtered.2 %>% group_by(Pango.lineage) %>% summarize(count.pango = n())
count.pango_mut.df <- mut.info.merged %>% group_by(Pango.lineage,mut) %>% summarize(count.pango_mut = n())

count.pango_mut.df.merged <- count.pango_mut.df %>% inner_join(count.pango.df,by="Pango.lineage")
count.pango_mut.df.merged <- count.pango_mut.df.merged %>% mutate(mut.freq = count.pango_mut / count.pango)

mut.interest.v <- count.pango_mut.df.merged %>% filter(mut.freq > 0.5) %>% pull(mut) %>% unique()

count.pango_mut.df.merged.filtered <- count.pango_mut.df.merged %>% filter(mut %in% mut.interest.v)
count.pango_mut.df.merged.filtered <- count.pango_mut.df.merged.filtered %>% mutate(mut.freq.binary = ifelse(mut.freq > 0.2,1,0))

count.pango_mut.df.merged.filtered <- count.pango_mut.df.merged.filtered %>% mutate(gene = gsub("_[^_]+","",mut), AA_change = gsub("[^_]+_","",mut), pos = gsub("[A-Za-z]","",AA_change) %>% as.numeric())

mut.spread <- count.pango_mut.df.merged.filtered %>% select(Pango.lineage,mut,gene,AA_change,mut.freq.binary) %>% spread(key = Pango.lineage, value = mut.freq.binary)
mut.spread[is.na(mut.spread)] <- 0

mut.interest.v2 <- mut.spread[apply(mut.spread[,4:ncol(mut.spread)],1,sd) > 0,]$mut

count.pango_mut.df.merged.filtered <- count.pango_mut.df.merged.filtered %>% filter(mut %in% mut.interest.v2)

count.pango_mut.df.merged.filtered <- count.pango_mut.df.merged.filtered %>%
  mutate(mut_mod = ifelse(pos %in% 69:70,"HV69-70del",as.character(AA_change)))


count.pango_mut.df.merged.filtered <- count.pango_mut.df.merged.filtered %>% group_by(Pango.lineage,gene,mut_mod) %>% summarize(mut.freq = mean(mut.freq)) %>% ungroup() %>% mutate(mut = paste(gene,mut_mod,sep="_"))

mut.order.v <- count.pango_mut.df.merged.filtered %>% select(gene,mut_mod,mut) %>% mutate(pos = gsub("\\-.+","",gsub("[A-Za-z]","",mut_mod))) %>% distinct(mut_mod,.keep_all =T) %>% arrange(gene,as.numeric(pos)) %>% pull(mut)

count.pango_mut.df.merged.filtered <- count.pango_mut.df.merged.filtered %>% mutate(mut = factor(mut,levels=rev(mut.order.v)))
#count.pango_mut.df.merged.filtered <- count.pango_mut.df.merged.filtered %>% filter(Pango.lineage3 != "BA.2")


count.pango_mut.df.merged.filtered <- count.pango_mut.df.merged.filtered %>% mutate(Pango.lineage = factor(Pango.lineage,levels=lineage.interest.v))

g <- ggplot(count.pango_mut.df.merged.filtered, aes(x = Pango.lineage, y = mut, fill = mut.freq))
g <- g + geom_tile()
g <- g + scale_fill_gradientn(colours=brewer.pal(9, "BuPu"),limits=c(0,1))
g <- g + theme_set(theme_classic(base_size = 10, base_family = "Helvetica"))
g <- g + theme(panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.background = element_blank()
)
g <- g + theme(
  legend.key.size = unit(0.3, 'cm'), #change legend key size
  legend.key.height = unit(0.3, 'cm'), #change legend key height
  legend.key.width = unit(0.3, 'cm'), #change legend key width
  legend.title = element_text(size=6), #change legend title font size
  legend.text = element_text(size=6))
g <- g + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
g <- g + labs(fill = "Freq.")
g <- g + xlab("") + ylab("")
g <- g + theme(axis.ticks=element_line(colour = "black"),
               axis.text=element_text(colour = "black"))
g

pdf(mutation_figure.name,width=3,height=8)
#pdf(mutation_figure.name,width=2,height=5)
plot(g)
dev.off()
}



