library(pgenlibr)
library(data.table)

# File paths
psam_path <- "/dcs05/legacy-dcl01-arking/data/UK_Biobank/static/Genetic/GenotypingArray/Imputed/plink/chr22.psam"
# Read in the sample data
sample_data <- fread(psam_path)

# uniqueN(variant_data$ID) # 7379615

pgen_path <- "/dcs05/legacy-dcl01-arking/data/UK_Biobank/static/Genetic/GenotypingArray/Imputed/plink/chr22.pgen"
pvar_path <- "/dcs05/legacy-dcl01-arking/data/UK_Biobank/static/Genetic/GenotypingArray/Imputed/plink/chr22.pvar"

# Initialize the .pgen file
pvar <- NewPvar(pvar_path)
pgen <- NewPgen(pgen_path, pvar=pvar)

# read a single variant
buf <- pgenlibr::Buf(pgen)
pgenlibr::Read(pgen, buf, 1)

# read a list of variants
## var.idxs (list of variants you would like to read)
geno_mat <- ReadList(pgen, 1:3, meanimpute=F)
# 1:3 correspond to the variant ids with index 1:3 in the pvar table

GetVariantCt(pvar)

dim(geno_mat) # 487409      3

ClosePgen(pgen)
ClosePvar(pvar)
