#Run this script from terminal using bash as interpreter
#Does not need file extension in filename.
#EXAMPLE: bash dataQC.sh filename1
#!/usr/bin/bash
BINARY_NAME=$1
PCAqc=$2
SCRIPT_PATH="./PipelineScripts"
RSCRIPTS=("CallRateIBD_QC.R" "check_heterozygosity_rate.R" "heterozygosity_outliers_list.R" "removeRelated.R")

#Check if PLINK binary files are found, if not abort script
if [ ! -e ${BINARY_NAME}.bed -o ! -e ${BINARY_NAME}.bim -o ! -e ${BINARY_NAME}.fam ]; then
    echo "Error: input files not found!"
    exit 1
fi
printf "Found binary files with prefix: $BINARY_NAME\n"
printf "`cat ${BINARY_NAME}.fam | wc -l` samples in ${BINARY_NAME}.fam\n\n"

#Check if script folder exists, if not create one and move R files there
if [ ! -d $SCRIPT_PATH ]; then
    echo "PipelineScripts folder not found. Creating new folder."
    mkdir $SCRIPT_PATH
    for s in "${RSCRIPTS[@]}"; do
        echo "Moving '$s' to $SCRIPT_PATH" 
        mv "$s" $SCRIPT_PATH
    done
else
echo "Found '$SCRIPT_PATH' folder, checking script files." 
fi

#Check for missing R scripts, then warn user which one(s) is/are missing.
if find $SCRIPT_PATH -maxdepth 0 -empty | read v; then echo "'$SCRIPT_PATH' is empty! Move R scripts to this directory.'"; fi

cd $SCRIPT_PATH
miss_count=0
for s in "${RSCRIPTS[@]}"; do
    if [[ ! -e "$s" ]]; then
        echo "$s - missing"
        ((miss_count++))
    fi
done
if [ $miss_count -gt 0 ]; then
    printf "Missing $miss_count R scripts. Check if filename has been altered.\nAborting...\n"
    exit 1
else
printf "All R scripts found in $SCRIPT_PATH\n\n" 
cd "$OLDPWD"
fi

##GWAS Pre-QC
#1. Sample QC - filter call rate <0.98
printf "Running pre-QC step 1: call rate filtering\n"
plink --allow-no-sex --bfile "$BINARY_NAME" --missing --make-bed --out ${BINARY_NAME}_missing

#R script to remove call rate (arguments: script function [1] = filter <98% call rate, .imiss file)
R --vanilla --slave --args 1 ${BINARY_NAME}_missing.imiss < $SCRIPT_PATH/CallRateIBD_QC.R
#Using -ge 2 to workaround header line, to check if call rate outliers exist 
#!/usr/bin/bash
if [ $(cat "list.calloutliers0.98.txt" | wc -l) -ge 2 ]; then
    echo " `tail -n +2 list.calloutliers0.98.txt | wc -l` Sample with call rate <0.98 found. Running PLINK to exclude." 
    plink --noweb --allow-no-sex --bfile "$BINARY_NAME" --remove list.call.outliers0.98.txt --make-bed --out ${BINARY_NAME}_Call98
else
printf "No samples below call rate threshold (<0.98) found.\n\n"
plink --noweb --allow-no-sex --bfile "$BINARY_NAME" --make-bed --out ${BINARY_NAME}_Call98
fi

#2. Sample QC 2 - sex mismatch (--check-sex)
printf "\nChecking for sample sex discrepancies.\n"
plink --noweb --allow-no-sex --bfile ${BINARY_NAME}_Call98 --check-sex --out ${BINARY_NAME}_Call98
#If PROBLEM found, run plink to remove samples with sex mismatches
if [ $(cat "${BINARY_NAME}_Call98.sexcheck" | grep -o PROBLEM | wc -l) -ge 1 ]; then
    echo "Found `grep -c PROBLEM ${BINARY_NAME}_Call98.sexcheck` samples with mismatched sex."
    grep PROBLEM ${BINARY_NAME}_Call98.sexcheck | tee list.sample.sexprobs.txt
    printf "Filtering these samples using PLINK --remove\n\n"
    plink --noweb --allow-no-sex --bfile ${BINARY_NAME}_Call98 --remove list.sample.sexprobs.txt --make-bed --out ${BINARY_NAME}_Call98.gender
else
printf "\nNo samples with conflicting sex.\n\n"
plink --noweb --allow-no-sex --bfile ${BINARY_NAME}_Call98 --make-bed --out ${BINARY_NAME}_Call98.gender
fi

#2a. Check sample heterozygosity
printf "\nChecking for sample heterozygosity outliers (>3 s.d.)\n"
plink --bfile ${BINARY_NAME}_Call98.gender --het --out R_check
Rscript --no-save $SCRIPT_PATH/check_heterozygosity_rate.R
Rscript --no-save $SCRIPT_PATH/heterozygosity_outliers_list.R
sed 's/"//g' "fail-het-qc.txt" | awk '{print $1, $2}' > "het_fail_ind_${BINARY_NAME}.txt"
if [ $(cat "het_fail_ind_${BINARY_NAME}.txt" | wc -l) -gt 1 ]; then
    printf "Found `tail -n +2 het_fail_ind_${BINARY_NAME}.txt | wc -l` heterozygosity outliers.\n\n"
    plink --allow-no-sex --bfile ${BINARY_NAME}_Call98.gender --remove het_fail_ind_${BINARY_NAME}.txt --make-bed --out ${BINARY_NAME}_noHet
else
printf "\nNo heterozygous outliers found.\n\n"
plink --noweb --allow-no-sex --bfile ${BINARY_NAME}_Call98.gender --make-bed --out ${BINARY_NAME}_noHet
fi
#3. Check sample IBD > 0.9
findIBD() {
    plink --noweb --allow-no-sex --bfile $1 --missing --genome --min $3 --make-bed --out $2
}

echo "Pruning high LD SNPs to calculate IBD"
plink --allow-no-sex --bfile ${BINARY_NAME}_noHet --indep-pairwise 50 5 0.2 --out ${BINARY_NAME}_noHet
plink --allow-no-sex --bfile ${BINARY_NAME}_noHet --extract ${BINARY_NAME}_noHet.prune.in --make-bed --out ${BINARY_NAME}_pruned

findIBD "${BINARY_NAME}_pruned" "${BINARY_NAME}_pihat0.9" 0.9
if [ $(cat "${BINARY_NAME}_pihat0.9.genome" | wc -l) -ge 2 ]; then
    printf "\nFound `tail -n +2 ${BINARY_NAME}_pihat0.9.genome | wc -l` samples with IBD >= 0.9\n"
    R --vanilla --slave --args ${BINARY_NAME}_pihat0.9.imiss ${BINARY_NAME}_pihat0.9.genome < $SCRIPT_PATH/removeRelated.R > remove.${BINARY_NAME}related.txt 
    plink --allow-no-sex --bfile ${BINARY_NAME}_noHet --remove "remove.${BINARY_NAME}related.txt" --make-bed --out ${BINARY_NAME}_ibd0.9
else
printf "\nNo highly related (PIHAT >= 0.9) samples found.\n"
plink --allow-no-sex --bfile ${BINARY_NAME}_noHet --make-bed --out ${BINARY_NAME}_ibd0.9
fi

#3a. Check sample IBD > 0.185
findIBD "${BINARY_NAME}_pruned" "${BINARY_NAME}_pihat0.185" 0.185
if [ $(cat "${BINARY_NAME}_pihat0.185.genome" | wc -l) -ge 2 ]; then
    printf "\nFound `tail -n +2 ${BINARY_NAME}_pihat0.185.genome | wc -l` samples with IBD >= 0.185\n"
    R --vanilla --slave --args ${BINARY_NAME}_pihat0.185.imiss ${BINARY_NAME}_pihat0.185.genome < $SCRIPT_PATH/removeRelated.R > remove.${BINARY_NAME}fam.txt 
    plink --allow-no-sex --bfile ${BINARY_NAME}_ibd0.9 --remove "remove.${BINARY_NAME}fam.txt" --make-bed --out ${BINARY_NAME}_nofam
else
printf "\nNo related (PIHAT >= 0.185) samples found.\n"
plink --allow-no-sex --bfile ${BINARY_NAME}_ibd0.9 --make-bed --out ${BINARY_NAME}_nofam
fi

#4. SNP QC - select for autosomal SNPs
plink --noweb --allow-no-sex --bfile ${BINARY_NAME}_nofam --autosome --make-bed --out ${BINARY_NAME}_autosome

#Intermediate file cleanup 
rm ${BINARY_NAME}{_nofam,_noHet,_ibd0.9,_Call98.gender,_Call98,_missing,_pihat*}.{bed,bim,fam} *.hh *.{imiss,lmiss} *.prune.{in,out} *_pruned.{bed,bim,fam,log}

#Merge - QC on both -> call rate, IBD 0.9, autosome
# plink --allow-no-sex --bfile "$datasetA" --bmerge "$datasetB" --make-bed --out AB_merged

printf "\nPre-QC completed.\n\n"

if [ ! -z ${PCAqc} ]; then
    if [ $PCAqc == "pca" ]; then
        echo "Running QC for PCA"
        bash pca-qc.sh ${BINARY_NAME}_autosome del
    fi
fi
