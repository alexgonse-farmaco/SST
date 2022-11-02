#Calculate SST PRS from CMC genotypes, using imputed data

## Load modules

module load bio/PLINK
module load R

### Set the path for your scripts and input files. Create name for your output file
## Paths
SCRIPTPATH='/external/rprshnas01/kcni/asegura/SST_Fernanda'
TARGETPATH='/external/rprshnas01/netdata_kcni/stlab/Alex_PRS/CMC_genotypes/imputed_annotated'
SUMSTATPATH='/external/rprshnas01/kcni/asegura/SST_Fernanda'
RESULTPATH='/external/rprshnas01/kcni/asegura/SST_Fernanda/results'
SUMSTATFILE='CommonMind_ge_DLPFC_naive.permuted.sstprs.summarystats_prs1_prs2_v2'
OUT='SSTPRS_cmcimp'

export fields='$16,$18,$17,$22,$12'
THRESHOLD1='0.1'
THRESHOLD2='1'
cd $RESULTPATH

for i in {1..22}; do
cat $TARGETPATH/cmcrod_chr$i.bim > $RESULTPATH/target_chr$i.bim
done

awk '{print '$fields'}'  $SUMSTATPATH/$SUMSTATFILE'.txt' > $RESULTPATH/sumstat.txt

Rscript $SCRIPTPATH/calculate_SST_PRS_Rscript1.R #execute script from Fernanda

for i in {1..22}; do #delete problematic SNPs by using a prevoulsy created file for CMC imputed?annotated data and clump
plink --bfile /external/rprshnas01/netdata_kcni/stlab/Alex_PRS/CMC_genotypes/imputed_annotated/cmcrod_chr$i --exclude $TARGETPATH/cmcrod-merge.missnp --make-bed --out $RESULTPATH/cmcsst_chr$i
plink --bfile cmcsst_chr$i --clump-p1 1 --clump-p2 1 --clump-r2 0.10 --clump-kb 250 --clump $RESULTPATH/SNPs_P_forclumping.txt --out $RESULTPATH/target_clumped_1_chr$i
done

Rscript  $SCRIPTPATH/calculate_SST_PRS_Rscript2_alex.R #execute scrpit from Fernanda, modified to work for 22 chrom in parallel

echo "qv"$THRESHOLD1" 0.00 "$THRESHOLD1 >q.ranges.txt
echo "qv"$THRESHOLD2" 0.00 "$THRESHOLD2 >>q.ranges.txt

for i in {1..22}; do
plink --allow-no-sex --bfile $TARGETPATH/cmcrod_chr$i --extract targetsnps.bim --exclude $TARGETPATH/cmcrod-merge.missnp --score $RESULTPATH/target_clumped_1_nothreshold.raw 1 2 3 header sum  --write-snplist --q-score-range q.ranges.txt target_clumped_1_nothreshold.snppvalues --out $RESULTPATH/$OUT'_chr'$i
done

Rscript $SCRIPTPATH/sumscores.R #execute script to get a merged file with both thresholds
