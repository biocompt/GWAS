input = "~/data/geno/REQUITE"
pheno = "~/data/feno/All_rbSTATs.txt"
ids = "~/data/feno/ids.txt"

process GWAS {
	"""
	plink2 --bfile $input --pheno $pheno --keep $ids --glm allow-no-covars --out RESULTS_GWAS
	awk '{ print \$3}' RESULTS_GWAS.rbAll.glm.linear > ids_SIGNIFICATIVE.txt
	plink2 --bfile $input --extract ids_SIGNIFICATIVE.txt --freq --out FREQUENCY
	"""
}
