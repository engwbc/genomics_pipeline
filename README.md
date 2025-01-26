# Main files:
1. `dataQC.sh` - Primarily sample QC including: relatedness (IBD) > 0.9 and >0.185, heterogeneity, sex discrepancies and call rate.
2. `pca-qh.sh` - Runs pairwise LD pruning and SNP QC based on PCA tool's recommendations, refer to [IPCAPS](https://scfbm.biomedcentral.com/articles/10.1186/s13029-019-0072-6)
3. `gwas-qc.sh` - Selects samples based on PCA results then runs additional QC steps required for association tests.

## PipelineScripts
These R scripts are sourced within the shell scripts. <br>
1. `CallRateIBD.R` - outputs a .txt file listing samples with call rate > 0.98
2. `check_heterozygosity_rate.R` - plots a histogram showing the heterozygosity rates obtained from PLINK `--het`
3. `heterozygosity_outliers_list.R` - creates a list of samples who have heterozygosity rates >3 standard deviations from the mean.
4. `removeRelated.R` - creates a list of samples with an IBD rate greater than a pre-determined threshold (PLINK ``--genome --min x``) and lowest call rate from the paired comparisons.
