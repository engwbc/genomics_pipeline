# GWAS script #
#!/usr/bin/bash
INPUTBIN=$1 #autosome plink files
KeepFile=$2 #keep.node#.txt
OUTFILE=$3 #custom name
SCRIPT_PATH="./PipelineScripts"

keep_name=${KeepFile%.txt} #remove trailing .txt
keep_name=${keep_name#*.} #remove leading keep.
tmp_name=${INPUTBIN}_${keep_name}

#Exclude outliers from PCA
plink --allow-no-sex --bfile ${INPUTBIN} --keep $KeepFile --make-bed --out ${tmp_name}

echo "Pruning high LD SNPs to calculate IBD"
plink --allow-no-sex --bfile ${tmp_name} --indep-pairwise 50 5 0.2 --out ${tmp_name}
plink --allow-no-sex --bfile ${tmp_name} --extract ${tmp_name}.prune.in --make-bed --out ${tmp_name}_pruned

#Simplify calling repeat plink commands by containing them into functions; change parameters here once
postQC() {
    plink --allow-no-sex --bfile $1 --mind 0.1 --geno 0.05 --hwe 1e-06 --maf 0.01 --make-bed --out $2
}

#Remove family members (PI_HAT > 0.185)
plink --noweb --allow-no-sex --bfile ${tmp_name}_pruned --missing --genome --min 0.185 --out ${tmp_name}_families
if [ $(cat "${tmp_name}_families.genome" | wc -l) -ge 2 ]; then
    printf "\nFound `tail -n +2 ${tmp_name}_families.genome | wc -l` samples with IBD >= 0.185\n"
    R --vanilla --slave --args ${tmp_name}_families.imiss ${tmp_name}_families.genome < $SCRIPT_PATH/removeRelated.R > remove.${tmp_name}families.txt
    echo "`head remove.${tmp_name}families.txt`" 
    plink --bfile ${tmp_name} --remove "remove.${tmp_name}families.txt" --make-bed --out ${OUTFILE}_nofam
    postQC ${OUTFILE}_nofam ${OUTFILE}_${keep_name} 
else
printf "\n\nNo additional related individuals found.\nProceeding with post-QC...\n"
postQC ${tmp_name} ${OUTFILE}_${keep_name} 
fi

printf "\nGenerating FasTLMM phenotype file\n"
awk '{print $1,$2,$6}' ${OUTFILE}_${keep_name}.fam > ${OUTFILE}_${keep_name}.phenotype.txt

rm ${tmp_name}* ${OUTFILE}_nofam*

echo "Post-QC complete"

#echo "Compressing output files..."
#tar -czvf ${OUTFILE}_${keep_name}.tar.gz ${OUTFILE}_${keep_name}.{bed,bim,fam}

exit 0
#Unmerged sample QC
#N93 - 188 samples (from 202 - removed 14; 3 from high/low heterozygosity and 11 from sex mismatch)
#N86 - 835 samples (from 863 - removed 28; 3 from IBD>=0.185, 3 from IBD>=0.9, 8 from heterozygosity and 14 from sex mismatch)

#Merged sample QC
#1021 samples (from 1026 - removed 5 from cluster group)

#Is it better to clear IBD >0.185 before merging datasets?