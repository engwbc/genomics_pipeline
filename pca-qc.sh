#!/usr/bin/bash

#Future plans: wait for IPCAPS to finish then continue this script based on argument value
INPUTBIN=$1
removeFiles=$2

if [ -z "${INPUTBIN}" ]; then
    printf "Error! Plink filename not given (e.g., dataset1 or NCD1)"
    exit 1
fi

#Principal Component Analysis
#5. LD pruning, remove high LD SNPs
plink --bfile ${INPUTBIN} --indep-pairwise 50 5 0.2 --out ${INPUTBIN}
echo "`wc -l ${INPUTBIN}.prune.in ${INPUTBIN}.prune.out`"

#5a. Exclude pruned SNPs
plink --bfile ${INPUTBIN} --extract ${INPUTBIN}.prune.in --make-bed --out ${INPUTBIN}_pruned

#5b. IPCAPS additional SNP filtering
plink --bfile ${INPUTBIN}_pruned --hwe 0.001 --mind 0.05 --geno 0.02 --maf 0.05 --make-bed --out ${INPUTBIN}_pca

echo "Generating sample list for IPCAPS"
awk '{print $1,$2,$6}' ${INPUTBIN}_pca.fam > list.group.$1.txt
echo "Substituting numeric phenotype category to Case/Control/Other"
awk '{gsub("1","Control",$3); gsub("2","Case",$3); gsub("-9","Other",$3)}1' list.group.$1.txt > list.group.tmp \
&& mv list.group.tmp list.group.$1.txt

if [ $removeFiles == "del" ]; then
    echo "Removing intermediate files"
    rm *.prune.{in,out} ${INPUTBIN}_pruned*
fi
echo "PCA QC completed"

##Extract selected samples
#plink --bfile pm4.call0.98.gender.autosome --keep gwas-sample-list.txt --make-bed --out data4gwas
