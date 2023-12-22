#bash prs_step123_microbiome.bash <target file location> <folder prefix> <basename list>
#bash prs_step123_microbiome.bash /mnt/data_schen_1/jingchun/dbgap_AI_project/dbGaP_dataset/AD_microglia_dbgap/phs000168/AD_merged/AD_All.QC1 ad_from_microbiome_ /mnt/data_schen_1/jingchun/gwas_data/Microbiome_GWAS_data/AD_microbiome_scores/microbiome1.list


targetfile=$1
folderprefix=$2
basenamelist=$3

while read celltype
do
echo $celltype start
/mnt/data_schen_1/jingchun/scz_plink_scores/PRSICE2.3.5/PRSice_linux --fastscore --base /mnt/data_schen_1/jingchun/gwas_data/Microbiome_GWAS_data/data.bris.ac.uk/datasets/22bqn399f9i432q56gt3wfhzlc/mGWAS/gwas/method_score/fgfp/addsnp/${celltype}_RNT_allchr.addsnp.txt.gz --chr CHR --snp SNP1 --bp BP --A1 A1 --A2 A2 --pvalue P --stat BETA --clump-kb 250 --clump-p 1 --clump-r2 0.1 --bar-levels 0.00000005,0.00001,0.001,0.1 --model add --out ${folderprefix}${celltype}/${folderprefix}${celltype} --all-score ${folderprefix}${celltype}/${folderprefix}${celltype} --binary-target T --target ${targetfile} --thread 1 --upper 1

/mnt/data_schen_1/jingchun/scz_plink_scores/PRSICE2.3.5/PRSice_linux --fastscore --base /mnt/data_schen_1/jingchun/gwas_data/Microbiome_GWAS_data/data.bris.ac.uk/datasets/22bqn399f9i432q56gt3wfhzlc/mGWAS/gwas/method_score/fgfp/addsnp/${celltype}_RNT_allchr.addsnp.txt.gz --chr CHR --snp SNP1 --bp BP --A1 A1 --A2 A2 --pvalue P --stat BETA --clump-kb 250 --clump-p 1 --clump-r2 0.1 --bar-levels 0.00000005,0.00001,0.001,0.1 --model add --out ${folderprefix}${celltype}/${folderprefix}${celltype} --all-score ${folderprefix}${celltype}/${folderprefix}${celltype} --binary-target T --target ${targetfile} --thread 1 --upper 1 --extract ${folderprefix}${celltype}/${folderprefix}${celltype}.valid

/mnt/data_schen_1/jingchun/scz_plink_scores/PRSICE2.3.5/PRSice_linux --base /mnt/data_schen_1/jingchun/gwas_data/Microbiome_GWAS_data/data.bris.ac.uk/datasets/22bqn399f9i432q56gt3wfhzlc/mGWAS/gwas/method_score/fgfp/addsnp/${celltype}_RNT_allchr.addsnp.txt.gz --chr CHR --snp SNP1 --bp BP --A1 A1 --A2 A2 --pvalue P --stat BETA --clump-kb 250 --clump-p 1 --clump-r2 0.1 --interval 0.00005 --lower 0.0001 --model add --out ${folderprefix}${celltype}/${folderprefix}${celltype}.best --binary-target T --target ${targetfile} --thread 1 --extract ${folderprefix}${celltype}/${folderprefix}${celltype}.valid --upper 0.5

echo $celltype done

done < ${basenamelist}

#plot
#Rscript /mnt/data_schen_1/jingchun/scz_plink_scores/PRSICE2.3.5/PRSice.R --prsice /mnt/data_schen_1/jingchun/scz_plink_scores/PRSICE2.3.5/PRSice_linux --base /mnt/data_schen_1/jingchun/gwas_data/Microbiome_GWAS_data/data.bris.ac.uk/datasets/22bqn399f9i432q56gt3wfhzlc/mGWAS/gwas/method_score/fgfp/addsnp/${p}_RNT_allchr.addsnp.txt.gz --chr CHR --snp SNP1 --bp BP --A1 A1 --A2 A2 --pvalue P --stat BETA --clump-kb 250 --clump-p 1 --clump-r2 0.1 --interval 0.00005 --lower 0.0001 --model add --out ad_from_microbiome_${p}/ad_from_microbiome_${p}.plot --binary-target T --target /mnt/data_schen_1/jingchun/dbgap_AI_project/dbGaP_dataset/AD_microglia_dbgap/phs000168/AD_merged/AD_All.QC1 --thread 1 --extract ad_from_microbiome_${p}/ad_from_microbiome_${p}.valid --upper 0.5
