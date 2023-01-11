geno = "~/data/geno/GENO"
pheno = "~/data/feno/All_rbSTATs.txt"
ids = "~/data/ids/ids.txt"
process GWAS {
	"""
	plink2 --bfile $geno --pheno $pheno --keep $ids --glm allow-no-covars --out RESULTS_GWAS
	"""
}
