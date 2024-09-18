

# info <- readRDS(runonce::download_file(
#   "https://figshare.com/ndownloader/files/36360900",
#   dir = "tmp-data", fname = "map_hm3.rds"))

map <- readRDS("./data/hapmap3/map_hm3.rds")
str(map)
typeof(map)

writeLines(map$rsid, "./data/hapmap3/euro_rsid.txt")

length(map$rsid)
chr_nums <- 1:22

for(chr in chr_nums){
  chr_map <- map[which(map$chr == chr),]
  writeLines(unique(chr_map$rsid), paste0("./data/hapmap3/euro_rsid_chr",chr,".txt"))
}

snp_count_by_chr <- map %>% select(chr, rsid) %>% distinct(.keep_all = TRUE) %>% group_by(chr) %>% summarise(count=n())


# chr count
# <int> <int>
#   1     1 87145
# 2     2 87521
# 3     3 73414
# 4     4 65490
# 5     5 66330
# 6     6 71046
# 7     7 58036
# 8     8 57084
# 9     9 48475
# 10    10 56549

# map$rsid[is.na(map$rsid)] # check there's no missing
