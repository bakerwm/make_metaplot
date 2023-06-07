# Date: 2023-06-07
# Author: Wang Ming
# version: 1.0

# Purpose
# extract AD/BD pairs from high through sequencing data
# the following graph explain the library structure
# ---(AD/BD)----[vec-1][attL][vec-2]----(BD/AD)---
#
# The sequences (5'->3'):
# 1. vec-1: TAGAACCCAGCTTTCTTGTACAAAGTGGTGAGCTTGGGCCCGTTTAAAC
# 2. vec-2: GATTATAAGGATGACGACGATAAAGGGCACTCGAGATATCTAGACCCAGCTTTCTTGTACAAAGTGGTGAGCTC (rev-comp)
# 3. attL: GTGCCAGGGCGTGCCCTTGAGTTCTCTCAGTTG (sens)
# 4. attL_rc: CAACTGAGAGAACTCAAGGGCACGCCCTGGCAC (rev-comp)

# How-To
# 1. fetch attL: extract all reads contain full-length 33-bp-attL or its rev-comp sequence
# 2. remove attL: remove attL and downstream sequence from each 3' end of read (minlen=69)
# 3. remove vector: remove AD/BD vector sequence from 3' end of each read (minlen=20)
# 4. pairing read12: output format, ID,read1_seq,read2_seq


CPU=16
################################################################################
## modules
## SE mode
function trim_attL() {
    local out_dir=$1
    local fq_in=$2
    local attL=$3 # sens, anti
    # check output 
    [[ ${fq_in} = *_2.fq.gz ]] && rd=2 || rd=1
    if [[ ${attL} = anti ]]
    then
        as="a"
        ad="GTGCCAGGGCGTGCCCTTGAGTTCTCTCAGTTG"
    else
        as="s"
        ad="CAACTGAGAGAACTCAAGGGCACGCCCTGGCAC"
    fi # anti/sens
    local suffix="_${rd}${as}"
    local fname=$(basename ${fq_in/.fq.gz})
    local fq_out="${out_dir}/${fname}${suffix}.fq.gz"
    local log="${out_dir}/${fname}${suffix}.cutadapt.log"
    [[ -f ${fq_out} ]] && echo "... file exists: ${fq_out}" && return 1
    [[ ! -d ${out_dir} ]] && mkdir -p ${out_dir}
    cutadapt -j ${CPU} --action trim --discard-untrimmed -O 33 -e 0.1 -m 69 -n 1 -a ${ad} -o ${fq_out} ${fq_in} > ${log}
    # echo ${fq_out} # file name
}
export -f trim_attL


function trim_vec() {
    local out_dir=$1
    local fq_in=$2
    local fq_out="${out_dir}/$(basename ${fq_in})"
    local log=${fq_out/.fq.gz/.cutadapt.log}
    local vec_1="TAGAACCCAGCTTTCTTGTACAAAGTGGTGAGCTTGGGCCCGTTTAAAC" # AD
    local vec_2="GATTATAAGGATGACGACGATAAAGGGCACTCGAGATATCTAGACCCAGCTTTCTTGTACAAAGTGGTGAGCTC" # BD
    [[ -f ${fq_out} ]] && echo "... file exists: ${fq_out}" && return 1
    [[ ! -d ${out_dir} ]] && mkdir -p ${out_dir}
    cutadapt -j ${CPU} --discard-untrimmed -O 12 -e 0.1 -m 20 -n 2 -a ${vec_1} -a ${vec_2} -o ${fq_out} ${fq_in} > ${log}
    # echo ${fq_out} # file name
}
export -f trim_vec




# Annotate fq, by gene_name/gene_id/NULL ...
# human library
function align_hg38() {
    local out_dir=$1
    local fq=$2 # fasta
    local idx="/data/yulab/hiseq001/data/genome/hg38/hisat2_index/hg38"
    local fname=$(basename ${fq/.fq.gz})
    local bam="${out_dir}/${fname}.bam"
    local bed="${bam/.bam/.bed}"
    local log="${out_dir}/${fname}.hisat2.log"
    # [[ -f ${bam} ]] && echo "file exists: ${bam}" && return 1
    [[ ! -d ${out_dir} ]] && mkdir -p ${out_dir}
    [[ ! -f ${bam} ]] && \
        hisat2 -p ${CPU} --very-sensitive --add-chrname -x ${idx} -U ${fq} 2> ${log} | \
        samtools view -Sub -F 0x4 -F 2048 -q 10 - | \
        samtools sort -@ ${CPU} -o ${bam} - && \
        samtools index ${bam}
    # convert bam to bed
    [[ ! -f ${bed} ]] && \
        bedtools bamtobed -i ${bam} > ${bed} #
    # echo ${bed}
}
export -f align_hg38


function anno_bed() {
    local out_dir=$1
    local in_bed=$2 # query
    local bname=$(basename ${in_bed/.bed/.anno.bed})
    local out_bed="${out_dir}/${bname}" # output 
    [[ ! -d ${out_dir} ]] && mkdir -p ${out_dir}
    # GRCh38, all genes
    local gene_bed="/data/yulab/hiseq001/user/wangming/0606_dnaseq/data/db/Homo_sapiens.GRCh38.106.gene.bed.gz"
    local anno_py="/data/yulab/hiseq001/user/wangming/0606_dnaseq/scripts/anno_bed.py"
    [[ ! -f ${out_bed} ]] && python ${anno_py} ${gene_bed} ${in_bed} ${out_bed}
}
export -f anno_bed


# annotate fastq, from 
function anno_fq() {
    local out_dir=$1
    local fq_in=$2
    local fname="$(basename ${fq_in/.fq.gz})"
    local bed="${out_dir}/${fname}.bed"
    [[ ! -d ${out_dir} ]] && mkdir -p ${out_dir}
    # 1. map to hg38 genome
    [[ ! -f ${bed} ]] && align_hg38 ${out_dir} ${fq_in}
    # 2. get annotation (gene_name)
    anno_bed ${out_dir} ${bed} 
    # 3. log
    echo "... save to file: ${bed}"
}
export -f anno_fq


# output common reads (by_id) from two fastq files
# save to new directory
function merge_fq_by_id() {
    local out_dir=$1
    local fq1=$2
    local fq2=$3
    [[ ! -d ${out_dir} ]] && mkdir -p ${out_dir}
    # file-1:
    out1="${out_dir}/$(basename ${fq1/.fq/.fa})"
    [[ -f ${out1} ]] && echo "... file exists: ${out1}" || \
        bioawk -cfastx 'FNR==NR {a[$name]=$seq; next} {if($name in a){print ">"$name"\n"a[$name]"\n"}}' ${fq1} ${fq2} | gzip > ${out1}
    # file-2:
    out2="${out_dir}/$(basename ${fq2/.fq/.fa})"
    [[ -f ${out2} ]] && echo "... file exists: ${out2}" || \
        bioawk -cfastx 'FNR==NR {a[$name]=$seq; next} {if($name in a){print ">"$name"\n"a[$name]"\n"}}' ${fq2} ${fq1} | gzip > ${out2}
}
export -f merge_fq_by_id


# input: bed6+pos+gene_name
# ouput: read_id,gene_name,gene_name
function merge_bed_by_id() {
    local out_dir=$1
    local bed1=$2
    local bed2=$3
    [[ ! -d ${out_dir} ]] && mkdir -p ${out_dir}
    # file-1:
    out1="${out_dir}/$(basename ${bed1/.bed/.read12.bed})"
    [[ ! -f ${out1} ]] && \
        awk 'FNR==NR {a[$4]=$8; next} {if($4 in a){print $4,a[$4],$8}}' ${bed1} ${bed2} > ${out1}
    # # file-2:
    # out1="${out_dir}/$(basename ${bed1/.bed/.read12.bed})"
    # [[ ! -f ${out2} ]] && \
    #     awk 'FNR==NR {a[$4]=$8; next} {if($4 in a){print $4,$8,a[$4]}}' ${bed2} ${bed1} > ${out2}
}
export -f merge_bed_by_id


# wrap all steps
function main() {
    local out_dir=$1
    local fq1=$2
    local fq2=$3
    # 1. extract/trim attL
    echo "[1/4] trimming attL"
    out_dir1="${out_dir}/1.trim_attL"
    trim_attL ${out_dir1} ${fq1} sens
    trim_attL ${out_dir1} ${fq1} anti
    trim_attL ${out_dir1} ${fq2} sens
    trim_attL ${out_dir1} ${fq2} anti

    # 2. trim vector
    echo "[2/4] trimming vector"
    out_dir2="${out_dir}/2.trim_vector"
    for fq in ${out_dir1}/*gz 
    do 
        trim_vec ${out_dir2} ${fq}
    done

    # 3. annotate reads
    echo "[3/4] annotating reads"
    out_dir3="${out_dir}/3.annotate_reads"
    for c in ${out_dir2}/*fq.gz 
    do 
        anno_fq ${out_dir3} ${c}
    done

    # 4. paring read12
    echo "[4/4] find read1,2 paires"
    out_dir4="${out_dir}/4.pairing_reads"
    r1_name=$(basename ${fq1/.fq.gz})
    r2_name=$(basename ${fq2/.fq.gz})
    bed1s="${out_dir3}/${r1_name}_1s.anno.bed"
    bed1a="${out_dir3}/${r1_name}_1a.anno.bed"
    bed2s="${out_dir3}/${r2_name}_2s.anno.bed"
    bed2a="${out_dir3}/${r2_name}_2a.anno.bed"
    ## bed format:
    merge_bed_by_id ${out_dir4} ${bed1s} ${bed2a}
    merge_bed_by_id ${out_dir4} ${bed1a} ${bed2s}

    # 5. finish
    echo "Finish!"


    # # 3. pairing read12: output as fasta format
    # out_dir3="${out_dir}/3.pairing_read12"
    # r1_name=$(basename ${fq1/.fq.gz})
    # r2_name=$(basename ${fq2/.fq.gz})
    # fq1s="${out_dir2}/${r1_name}_1s.fq.gz"
    # fq1a="${out_dir2}/${r1_name}_1a.fq.gz"
    # fq2s="${out_dir2}/${r2_name}_2s.fq.gz"
    # fq2a="${out_dir2}/${r2_name}_2a.fq.gz"
    # ## 1s + 2a
    # merge_fq_by_id ${out_dir3} ${fq1s} ${fq2a}
    # ## 1a + 2s
    # merge_fq_by_id ${out_dir3} ${fq1a} ${fq2s}

}
export -f main

[[ $# -lt 3 ]] && echo "Usage: find_ad_bd.sh <out_dir> <fq1> <fq2>" && exit 1
main $1 $2 $3


# ## 1. fetch attL
# ## 2. remove attL
# out_dir1="results3/1.trim_attL"
# sens="GTGCCAGGGCGTGCCCTTGAGTTCTCTCAGTTG"
# anti="CAACTGAGAGAACTCAAGGGCACGCCCTGGCAC"
# r1="data/raw_data/ATAC_100mixlib_3hIP_rep1_1.fq.gz"
# r2="data/raw_data/ATAC_100mixlib_3hIP_rep1_2.fq.gz"
# mkdir -p ${out_dir1}

# ## 1s,read1,sens
# out_1s="$out_dir1/$(basename ${r1/.fq.gz/_1s.fq.gz})"
# log1s="$out_dir1/$(basename ${r1/.fq.gz/_1s.cutadapt.log})"
# [[ ! -f ${out_1s} ]] && \
#     cutadapt -j 8 --action trim --discard-untrimmed -O 33 -e 0.1 -m 69 -n 1 -a ${sens} -o ${out_1s} ${r1} > ${log1s}
# ## 2a,read2,anti
# out_2a="$out_dir1/$(basename ${r2/.fq.gz/_2a.fq.gz})"
# log2a="$out_dir1/$(basename ${r1/.fq.gz/_2a.cutadapt.log})"
# [[ ! -f ${out_2a} ]] && \
#     cutadapt -j 8 --action trim --discard-untrimmed -O 33 -e 0.1 -m 69 -n 1 -a ${anti} -o ${out_2a} ${r2} > ${log2a}

# ## 1a,read1,anti
# out_1a="$out_dir1/$(basename ${r1/.fq.gz/_1a.fq.gz})"
# log1a="$out_dir1/$(basename ${r1/.fq.gz/_1a.cutadapt.log})"
# [[ ! -f ${out_1a} ]] && \
#     cutadapt -j 8 --action trim --discard-untrimmed -O 33 -e 0.1 -m 69 -n 1 -a ${anti} -o ${out_1a} ${r1} > ${log1a}
# ## 2s,read2,sens
# out_2s="$out_dir1/$(basename ${r2/.fq.gz/_2s.fq.gz})"
# log2s="$out_dir1/$(basename ${r1/.fq.gz/_2s.cutadapt.log})"
# [[ ! -f ${out_2s} ]] && \
#     cutadapt -j 8 --action trim --discard-untrimmed -O 33 -e 0.1 -m 69 -n 1 -a ${sens} -o ${out_2s} ${r2} > ${log2s}


# ################################################################################
# ## 3. remove vector
# out_dir2="results3/2.trim_vector"
# vec_1="TAGAACCCAGCTTTCTTGTACAAAGTGGTGAGCTTGGGCCCGTTTAAAC" # AD
# vec_2="GATTATAAGGATGACGACGATAAAGGGCACTCGAGATATCTAGACCCAGCTTTCTTGTACAAAGTGGTGAGCTC" # BD
# mkdir -p ${out_dir2}

# function trim_vec() {
#     local fq=$1
#     local out_dir=$2
#     local out="${out_dir}/$(basename ${fq})"
#     local log=${out/.fq.gz/.cutadapt.log}
#     [[ -f ${out} ]] && echo "file exists: ${out}" && return 1
#     mkdir -p ${out_dir}
#     cutadapt -j 8 --discard-untrimmed -O 12 -e 0.1 -m 20 -n 2 -a ${vec_1} -a ${vec_2} -o ${out} ${fq} > ${log}
# }
# export -f trim_vec

# for a in ${out_dir1}/*gz
# do
#     trim_vec $a ${out_dir2}	
# done


# ################################################################################
# ## 4. pairing read12
# out_dir3="results3/3.pairing_read12"
# mkdir -p ${out_dir3}
# ## 1s+2a
# bash merge_fq.sh ${out_dir3} ${out_dir2}/*1s.fq.gz ${out_dir2}/*2a.fq.gz
# ## 1a+2s
# bash merge_fq.sh ${out_dir3} ${out_dir2}/*1a.fq.gz ${out_dir2}/*2s.fq.gz


# # 
# # 
# # cutadapt -j 8 --rc --action trim --discard-untrimmed -O 33 -e 0.1 -m 69 -n 2 -a GTGCCAGGGCGTGCCCTTGAGTTCTCTCAGTTG -o result2/ATAC_100mixlib_3hIP_rep1_3_1.fq.gz data/raw_data/ATAC_100mixlib_3hIP_rep1_1.fq.gz > ATAC_100mixlib_3hIP_rep1_3_1.cutadapt.log 
# # 
# # ## 2a,read2,anti
# # 
# # 
# # ## 1a,read1,anti
# # ## 2s,read2,sens
# # 
# # cutadapt -j 8 --rc --action trim --discard-untrimmed -O 33 -e 0.1 -m 69 -n 2 -a CAACTGAGAGAACTCAAGGGCACGCCCTGGCAC -o result2/ATAC_100mixlib_3hIP_rep1_5_1.fq.gz data/raw_data/ATAC_100mixlib_3hIP_rep1_1.fq.gz > ATAC_100mixlib_3hIP_rep1_5_1.cutadapt.log 
# # 
# # cutadapt -j 8 --rc --action trim --discard-untrimmed -O 33 -e 0.1 -m 69 -n 2 -a GTGCCAGGGCGTGCCCTTGAGTTCTCTCAGTTG -o result2/ATAC_100mixlib_3hIP_rep1_3_2.fq.gz data/raw_data/ATAC_100mixlib_3hIP_rep1_2.fq.gz > ATAC_100mixlib_3hIP_rep1_3_2.cutadapt.log
# # 
# # cutadapt -j 8 --rc --action trim --discard-untrimmed -O 33 -e 0.1 -m 69 -n 2 -a CAACTGAGAGAACTCAAGGGCACGCCCTGGCAC -o result2/ATAC_100mixlib_3hIP_rep1_5_2.fq.gz ./data/raw_data/ATAC_100mixlib_3hIP_rep1_2.fq.gz > ATAC_100mixlib_3hIP_rep1_5_2.cutadapt.log 
