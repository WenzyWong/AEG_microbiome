#####################################################
# AEG microbiome
# Yunzhe Wang
#
# Step 00. Extracting bacterial reads from RNA-seq
# using bacNeo --bacc
#
#####################################################
DIR_CLEAN="/data/yzwang/project/AEG_seiri/CleanData/"
INDEX_HISAT2="/data/yzwang/reference/hisat/hg38"
DIR_OUT="/data/yzwang/project/AEG_seiri/bacNeo_out/"

mkdir -p "${DIR_OUT}"

ls "${DIR_CLEAN}" | while read sample
do
    fq1="${DIR_CLEAN}/${sample}/${sample}_1.fq.gz"
    fq2="${DIR_CLEAN}/${sample}/${sample}_2.fq.gz"
    
    # Run bacc
    out_bacc="${DIR_OUT}/bacc/${sample}"
    mkdir -p "${DIR_OUT}"
    bacNeo --bacc -1 "${fq1}" -2 "${fq2}" -m RNA -r "${INDEX_HISAT2}" -o "${out_bacc}" -l s -l g -t 64
done

# Extract bacc normalized matrix
bacNeo --extract-matrix -d "${DIR_OUT}/bacc" -l s -l g -n abundance