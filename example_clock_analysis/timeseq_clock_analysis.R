#load packages
library("tidyverse")
library("data.table")
library('zoo')
library('ggpubr')

# Define working dir
wd <- "~/Desktop/TIME-Seq_Example/"

setwd(wd)

### Reading in coverage files and organizing data into matrixes

# Loading clock sites
clock_name <- 'TIME-Seq_Mouse_Blood_Clock.tsv'
clock <- read_tsv(paste('./clocks/',clock_name,sep=""))
cpgs_list <- clock$name

setwd(paste(wd,"coverage_files/",sep=""))

files <- str_subset(list.files(), '.+bismark.cov')

##Note: the names below for SE analysis. 
## For PE analysis, relace ".nonCG_filtered.bismark.cov" with "_bismark_merged.gz.bismark.cov"
sample <- str_match(files[1],pattern = "(.+).nonCG_filtered.bismark.cov")[,2]

#gathering methylation data and coverage data
reps_samples <- fread(files[1]) %>% 
  unite(col = "cpg",c(V1,V2) , sep="_") %>%
  filter(cpg %in% cpgs_list) %>%
  mutate(coverage=V5+V6) %>%
  select(cpg,V4,coverage)

mm <- reps_samples %>% select(cpg,V4) %>% rename(!!sample:=V4)
cm <- reps_samples %>% select(cpg,coverage) %>% rename(!!sample:=coverage)

n=1  

#gathering methylation data and coverage data for each file and making matrixes
for (f in files[2:length(files)]) {
  ##Note: the names below for SE analysis. 
  ## For PE analysis, relace ".nonCG_filtered.bismark.cov" with "_bismark_merged.gz.bismark.cov"
  sample <- str_match(f,pattern = "(.+).nonCG_filtered.bismark.cov")[,2]
  rep <- fread(f) %>%
    unite(col = "cpg",c(V1,V2) , sep="_") %>%
    filter(cpg %in% cpgs_list) %>%
    mutate(coverage=V5+V6) %>%
    select(cpg,V4,coverage) 
  
  mm <- rep %>% select(cpg,V4) %>% rename(!!sample:=V4) %>% full_join(mm,by="cpg")
  cm <- rep %>% select(cpg,coverage) %>% rename(!!sample:=coverage) %>% full_join(cm,by="cpg")
  
}

#tranposing methylation matrix
t.mm <- data.table::transpose(mm,keep.names = "sample_ID",make.names = T) 

# replacing NA coverage with 0
cm.f <- cm %>% mutate(across(where(is.numeric),
                             .fns = function(x) as.double(x))) %>%
  mutate(across(where(is.double),.fns = function(x) if_else(is.na(x)==T,0,x))) 

print(mm[1:5,1:5])
print(cm[1:5,1:5])


### Assessing low coverage at clock CpGs 
cutoff <- length(cm.f$cpg)*0.1

# Assessing how many clock CpGs samples have low coverage at.
low_coverages <- colSums(apply(cm.f[,-1], 2, function(x){as.numeric(x)<10})) %>%
  as.data.frame() %>%
  rownames_to_column(var = 'sample_ID') %>% 
  filter(sample_ID!="cpg") %>% 
  rename("lowcovs"=".") 

g <- low_coverages %>% 
  ggplot() +
  geom_histogram(aes(lowcovs)) +
  geom_vline(xintercept = cutoff, color="red", linetype="dashed") +
  theme_classic()
g

# Define samples to exclude
exclude <- filter(low_coverages, lowcovs>cutoff)$sample_ID
print(paste("Sample excluded:", exclude))

# Filter methylation matrix
t.mm.f <- filter(t.mm, sample_ID %in% exclude == F)

# Clock Analysis
setwd(wd)

# reading in sample info
sampleInfo <- read_csv(paste(wd,'/sampleInfo.csv',sep="")) %>% mutate(sample_ID=as.factor(sample_ID)) 

# reading in clock to analyze
clock_name <-  'TIME-Seq_Mouse_Blood_Clock.tsv'
clock <- read_tsv(paste('./clocks/',clock_name,sep="")) %>%
  rename("cpg"="name")

# coefficients
intercept <- filter(clock, cpg =='(Intercept)')$coefficient
a <- filter(clock, cpg =='a')$coefficient
c <- filter(clock, cpg =='c')$coefficient

# mean imputation in matrix to replace NAs.
t.mm.noNas <- cbind(tibble(sample_ID=t.mm.f$sample_ID), na.aggregate(t.mm.f[,-1]))

# tranposing back to CpG rowwise
mm.noNAs <- data.table::transpose(t.mm.noNas,keep.names = "cpg",make.names = T)

# reading in methylation matrix and melting
mm_clock <- mm.noNAs %>%
  melt(id.vars = "cpg") %>%
  rename("sample_ID"="variable")

# age predictions + deltaAge + sampleInfo
predictions <- clock  %>%
  filter(cpg %in% c('(Intercept)',"a","c")==F) %>% 
  full_join(mm_clock, by=c('cpg')) %>%
  mutate(weighted_mCpG = coefficient*(value)) %>%
  group_by(sample_ID) %>%
  summarise(sumWeighted_mCpG = sum(weighted_mCpG)) %>%
  mutate(b=sumWeighted_mCpG+intercept) %>%
  mutate(PredictedAge=a*b+c) %>%
  left_join(sampleInfo, by=c("sample_ID")) %>%
  mutate(deltaAge=PredictedAge-Age)
head(predictions)  


### Plot predictions
g <- predictions %>%
  ggplot() +
  geom_abline(slope=1,linetype="dashed",color="grey") +
  geom_point(aes(Age, PredictedAge),
             size=2,
             shape=21,
             alpha=0.5,
             fill="red") +
  theme_classic() +
  xlab('Chronological age (months)') +
  ylab('Predicted age (months)') +
  stat_cor(aes(Age, PredictedAge)) +
  xlim(0,NA) + ylim (0,NA)
g
