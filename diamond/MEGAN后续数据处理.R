
setwd("~/Desktop/Shared_Folder/metagenome_owe/Megagenome_learing/diamond/")

#----Taxonomy物种注释表格整理#-----
# --导入数据
library(tidyverse)
library(phyloseq)
# BiocManager::install("phyloseq")
path = "./"

fl.0 =dir(path, pattern = c("Taxonomy.txt"), full.names = TRUE, ignore.case = TRUE)

for (i in 1:length(fl.0)) {
  file = fl.0[i]
  samp.name <- sapply(strsplit(basename(file), "Taxonomy.txt"), `[`, 1)
  
  
  otut <- read.delim(file,header = F)
  
  #-去除MEGAN给出的不同分类等级的统计数据
  otut.1 <- otut[!str_detect(otut$V1, "root"),]
  
  head(otut.1)
  dim(otut.1)
  # 去除重复的基因，重复的注释，去掉,注意这里去掉的基因和注释信息都是一样的
  otut.1 %>% distinct( V1,V2, .keep_all = TRUE) %>%
    dim()
  otut.2 <- otut.1 %>% distinct( V1,V2, .keep_all = TRUE) 
  # 统计count
  
  otut.3 <- as.data.frame(table(otut.2$V2))
  
  head(otut.3)
  colnames(otut.3) <- c("ASV",samp.name)
  
  write.csv(otut.3,paste(samp.name,"_count.csv",sep = ""),row.names = F)
  
}

# 不同样本物种注释结果合并
path <- "./"
fl.1 =dir(path, pattern = c("_count.csv"), full.names = TRUE, ignore.case = TRUE)


if (length(fl.1) > 2) {
  data <- read.csv(fl.1[1]) %>% full_join(read.csv(fl.1[2 ]),by = "ASV")
  for (i in 3:length(fl.1)) {
    data <- data %>%
      full_join(read.csv(fl.1[i ]),by = "ASV")
  } 
  
  }else if(length(fl.1) == 2 ) {
  data <- read.csv(fl.1[1]) %>% full_join(read.csv(fl.1[2 ]),by = "ASV")
}

head(data)

write.csv(data,"./otu_count.csv",row.names = F)

#--注释文件

otutax = separate(data, col = ASV,
         into = c("subKingdom","sub2Kingdom","Kingdom","Phylum","subPhylum","Class","Order","Family","Genus","Species"),
         sep = ";",
         remove = F)
head(otutax)
otutax$OTUname = paste("ASV_",1:length(data$ASV),sep = "")


otu = otutax[,c(colnames(data)[-1])]
row.names(otu) = otutax$OTUname

tax = otutax[,c(4,5,7,8,9,10,11)]
row.names(tax) = otutax$OTUname
head(tax)


ps <- phyloseq(otu_table(as.matrix(otu),taxa_are_rows = T),
         tax_table(as.matrix(tax)))

saveRDS(ps,"ps_tax.rds")



#----eggnog功能表格整理#-----


fl.0 =dir(path, pattern = c("eggnog.txt"), full.names = TRUE, ignore.case = TRUE)
i= 1
for (i in 1:length(fl.0)) {
  file = fl.0[i]
  samp.name <- sapply(strsplit(basename(file), "eggnog.txt"), `[`, 1)
  
  
  otut <- read.delim(file,header = F)

  # 去除重复的基因，重复的注释，去掉,注意这里去掉的基因和注释信息都是一样的
  otut %>% distinct( V1,V2, .keep_all = TRUE) %>%
    dim()
  otut.2 <- otut %>% distinct( V1,V2, .keep_all = TRUE) 
  # 统计count
  
  otut.3 <- as.data.frame(table(otut.2$V2))
  
  head(otut.3)
  colnames(otut.3) <- c("eggnog",samp.name)
  
  write.csv(otut.3,paste(samp.name,"_eggnogcount.csv",sep = ""),row.names = F)
  
}





# 不同样本物种注释结果合并
path <- "./"
fl.1 =dir(path, pattern = c("_eggnogcount.csv"), full.names = TRUE, ignore.case = TRUE)


if (length(fl.1) > 2) {
  data <- read.csv(fl.1[1]) %>% full_join(read.csv(fl.1[2 ]),by = "eggnog")
  for (i in 3:length(fl.1)) {
    data <- data %>%
      full_join(read.csv(fl.1[i ]),by = "eggnog")
  } 
  
}else if(length(fl.1) == 2 ) {
  data <- read.csv(fl.1[1]) %>% full_join(read.csv(fl.1[2 ]),by = "eggnog")
}

head(data)

write.csv(data,"./eggnog_count.csv",row.names = F)


#--eggnog注释文件整理

otutax = separate(data, col = eggnog,
                  into = c("eggNOG.id","eggNOG.group","eggNOG.descrip","COG1192"),
                  sep = ";",
                  remove = F)
head(otutax)
otutax$OTUname = paste("eggnog_",1:length(data$eggnog),sep = "")


otu = otutax[,c(colnames(data)[-1])]
row.names(otu) = otutax$OTUname

tax = otutax[,c(2,3,4,5)]
row.names(tax) = otutax$OTUname

head(tax)


ps <- phyloseq(otu_table(as.matrix(otu),taxa_are_rows = T),
               tax_table(as.matrix(tax)))
ps
saveRDS(ps,"ps_eggnog.rds")

#--INTERPRO2GO功能注释数据整理#-----------


fl.0 =dir(path, pattern = c("INTERPRO2GO.txt"), full.names = TRUE, ignore.case = TRUE)
i= 1
for (i in 1:length(fl.0)) {
  file = fl.0[i]
  samp.name <- sapply(strsplit(basename(file), "INTERPRO2GO.txt"), `[`, 1)
  
  
  otut <- read.delim(file,header = F)
  
  # 去除重复的基因，重复的注释，去掉,注意这里去掉的基因和注释信息都是一样的
  otut %>% distinct( V1,V2, .keep_all = TRUE) %>%
    dim()
  otut.2 <- otut %>% distinct( V1,V2, .keep_all = TRUE) 
  # 统计count
  
  otut.3 <- as.data.frame(table(otut.2$V2))
  
  head(otut.3)
  colnames(otut.3) <- c("INTERPRO2GO",samp.name)
  
  write.csv(otut.3,paste(samp.name,"_INTERPRO2GOcount.csv",sep = ""),row.names = F)
  
}





# 不同样本物种注释结果合并
path <- "./"
fl.1 =dir(path, pattern = c("_INTERPRO2GOcount.csv"), full.names = TRUE, ignore.case = TRUE)


if (length(fl.1) > 2) {
  data <- read.csv(fl.1[1]) %>% full_join(read.csv(fl.1[2 ]),by = "INTERPRO2GO")
  for (i in 3:length(fl.1)) {
    data <- data %>%
      full_join(read.csv(fl.1[i ]),by = "INTERPRO2GO")
  } 
  
}else if(length(fl.1) == 2 ) {
  data <- read.csv(fl.1[1]) %>% full_join(read.csv(fl.1[2 ]),by = "INTERPRO2GO")
}

head(data)

write.csv(data,"./INTERPRO2GO_count.csv",row.names = F)


#--INTERPRO2GO注释文件整理

otutax = separate(data, col = INTERPRO2GO,
                  into = c("INTERPRO2GO","GO.group","GO.sub","IPR"),
                  sep = ";",
                  remove = F)
head(otutax)
otutax$OTUname = paste("INTERPRO2GO_",1:length(data$INTERPRO2GO),sep = "")


otu = otutax[,c(colnames(data)[-1])]
row.names(otu) = otutax$OTUname

tax = otutax[,c(2,3,4,5)]
row.names(tax) = otutax$OTUname

head(tax)


ps <- phyloseq(otu_table(as.matrix(otu),taxa_are_rows = T),
               tax_table(as.matrix(tax)))
ps
saveRDS(ps,"ps_INTERPRO2GO.rds")

#--SEED功能注释整理#----------



fl.0 =dir(path, pattern = c("SEED.txt"), full.names = TRUE, ignore.case = TRUE)
i= 1
for (i in 1:length(fl.0)) {
  file = fl.0[i]
  samp.name <- sapply(strsplit(basename(file), "SEED.txt"), `[`, 1)
  
  
  otut <- read.delim(file,header = F)
  
  # 去除重复的基因，重复的注释，去掉,注意这里去掉的基因和注释信息都是一样的
  otut %>% distinct( V1,V2, .keep_all = TRUE) %>%
    dim()
  otut.2 <- otut %>% distinct( V1,V2, .keep_all = TRUE) 
  # 统计count
  
  otut.3 <- as.data.frame(table(otut.2$V2))
  
  head(otut.3)
  colnames(otut.3) <- c("SEED",samp.name)
  
  write.csv(otut.3,paste(samp.name,"_SEEDcount.csv",sep = ""),row.names = F)
  
}





# 不同样本物种注释结果合并
path <- "./"
fl.1 =dir(path, pattern = c("_SEEDcount.csv"), full.names = TRUE, ignore.case = TRUE)


if (length(fl.1) > 2) {
  data <- read.csv(fl.1[1]) %>% full_join(read.csv(fl.1[2 ]),by = "SEED")
  for (i in 3:length(fl.1)) {
    data <- data %>%
      full_join(read.csv(fl.1[i ]),by = "SEED")
  } 
  
}else if(length(fl.1) == 2 ) {
  data <- read.csv(fl.1[1]) %>% full_join(read.csv(fl.1[2 ]),by = "SEED")
}

head(data)

write.csv(data,"./SEED_count.csv",row.names = F)


#--SEED注释文件整理

otutax = separate(data, col = SEED,
                  into = c("SEED","SEED.group","gene.id","SEED.discrip"),
                  sep = ";",
                  remove = F)
head(otutax)
otutax$OTUname = paste("SEED_",1:length(data$SEED),sep = "")


otu = otutax[,c(colnames(data)[-1])]
row.names(otu) = otutax$OTUname

tax = otutax[,c(2,3,4,5)]
row.names(tax) = otutax$OTUname

head(tax)


ps <- phyloseq(otu_table(as.matrix(otu),taxa_are_rows = T),
               tax_table(as.matrix(tax)))
ps
saveRDS(ps,"ps_SEED.rds")


