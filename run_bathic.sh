#!/usr/bin/env bash 


## variables
# fq_dir="/data1/yuyang/data/seq_data/YY353_20230629/results/"
# clean_dir="/data1/yuyang/yeziling/20230715_batHiC/00hictrim"
# hicpro_dir="/data1/yuyang/yeziling/20230715_batHiC/03hicresult"
# hic_dir="/data1/yuyang/yeziling/20230715_batHiC/06hic"



################################################################################
# The following scripts were designed for BAT Hi-C data processing             #
# see: https://doi.org/10.1093/nar/gkw809 for more details                     #
#                                                                              #
# There are 11 steps described in the original publication,                    #
#                                                                              #
# in brief                                                                     #
# 1. scan and trim bridge linker: trimLinker                                   #
# 2. run HiC-Pro                                                               #
# 3. Filter for valid interactions.                                            #
# 4. Parse the uniquely valid contacts pairs                                   #
# 5. Use ICE to normalize the contact matrix                                   #
# 6. Convert *allValidPairs files into Juice file .hic for visualization       #
################################################################################

# Project directory structure 
#
# Project
# ├── 01.raw_data
# │   ├── HiC_mHap_DMSO_2hr_rep1_1.fq.gz
# │   └── HiC_mHap_DMSO_2hr_rep1_2.fq.gz
# ├── 02.clean_data
# │   └── HiC_mHap_DMSO_2hr_rep1
# │       ├── HiC_mHap_DMSO_2hr_rep1_1.valid.fastq
# │       ├── HiC_mHap_DMSO_2hr_rep1_2.valid.fastq
# │       └── HiC_mHap_DMSO_2hr_rep1.trim.stat
# ├── 03.hicpro_results
# │   ├── bowtie_results
# │   ├── config-hicpro.txt
# │   ├── hic_results
# │   ├── logs
# │   ├── rawdata -> ...path-to-clean_data
# │   └── tmp
# ├── 04.hic_files
# │   └── HiC_mHap_DMSO_2hr_rep1.hic
# ├── 05.TAD_files
# │   └── HiC_mHap_DMSO_2hr_rep1
# ├── 06.juicer_files
# │   ├── 01.contact
# │   ├── 02.tad
# │   ├── 03.loop
# │   └── 04.compartments
# └── 07.report

################################################################################
# Change the following variables according to your system
## global variables
N_CORES=12
## software
HICPRO_PATH="/data1/yuyang/wangm/biosoft/HiC-Pro_3.1.0/" # HiC-Pro: PATH/bin/HiC-Pro
ChIAPET2_PATH="/data1/yuyang/wangm/biosoft/ChIA-PET2/" # trimLinker: PATH/bin/trimLinker
JUICER_PATH="/data1/yuyang/wangm/biosoft/juicer/" # juicer_tools: PATH/juicer_tools.jar
# hicpro_sif="/data1/yuyang/data/docker/hicpro_slurm.sif"
## customized variables for mouse ## 
CHROM_SIZES="/data1/yuyang/data/genome/mm10/bigZips/mm10.chrom.sizes"
RES_FRAG="/data1/yuyang/wangm/biodata/hic/restriction_sites/mm10_Alu1.bed"
################################################################################


## 1. trimlinker
function trim_linker() {
    # see ChIA-PET2/bin/trimLinker
    # scan and trim bridge
    local fq1=$1
    local fq2=$2
    local out_dir=$3
    local trimlinker="${ChIAPET2_PATH}/bin/trimLinker"
    [[ ! -f ${trimlinker} ]] && echo "    file not exists, see ChIA-PET2, ${trimlinker}" && return 1
    # [[ ! ${fq1} = *_1.fq.gz ]] && echo "not fq1: ${fq1}" && return 1
    [[ ! -f ${fq1} ]] && echo "fq1 not exists: ${fq1}" && return 1
    [[ ! -f ${fq2} ]] && echo "fq2 not exists: ${fq2}" && return 1
    local fq_name="$(basename "${fq1/_1.fq.gz}")"
    local sub_dir="${out_dir}/${fq_name}"
    [[ ! -d ${sub_dir} ]] && mkdir -p "${sub_dir}"
    local fq1_clean="${sub_dir}/${fq_name}_1.valid.fastq"
    local fq2_clean="${sub_dir}/${fq_name}_2.valid.fastq"
    if [[ -f ${fq1_clean} && -f ${fq2_clean} ]]
    then
        echo "    trim_linker skipped, file exists: ${sub_dir}" 
    else
        ${trimlinker} \
            -t ${N_CORES} -m 1 -k 1 -l 16 -e 0 \
            -o "${sub_dir}" -n "${fq_name}" \
            -A ACGCGATATCTTATC -B AGTCAGATAAGATAT \
            "${fq1}" "${fq2}"
    fi
}
export -f trim_linker 


## 2. run HiC-Pro, singularity
function run_hicpro() {
    local config=$1
    local data_dir=$2
    local out_dir=$3
    [[ ! -f ${config} ]] && echo "HiC-Pro config file not exists: ${config}"  && return 1
    [[ ! -d ${data_dir} ]] && echo "fastq data directory not exists: ${data_dir}" && return 1
    # [[ ! -d ${out_dir} ]] && mkdir -p ${out_dir} # !! do not need
    # check fastq files in data
    ## fixed pattern here: _1.valid.fastq, _2.valid.fastq
    n_fq=$(ls ${data_dir}/*/* | grep -c _1.valid.fastq) # number of fastq files, HiC-Pro skipped
    [[ ${n_fq} == 0 ]] && echo "error: no fastq files found: ${data_dir}/*/*gz" && return 1
    # check if HiC-Pro finished or not, allValidPairs
    echo "    >>>Check if HiC-Pro finished or not<<<"
    icheck=0
    for i in ${data_dir}/*/*_1.valid.fastq 
    do 
        iname=$(basename ${i/_1.valid.fastq})
        vp=${out_dir}/hic_results/data/${iname}/${iname}.allValidPairs
        [[ ! -f ${vp} ]] && icheck=1
        [[ -f ${vp} ]] && itag="ok" || itag="failed"
        printf "%10s - %s\n" ${itag} ${iname}
    done
    # using singularity for HiC-Pro
    local hicpro_cmd="${HICPRO_PATH}/bin/HiC-Pro"
    if [[ -f ${hicpro_sif} ]] 
    then 
        apptainer run --bind /var:/var,/data1:/data1 ${hicpro_sif} \
            HiC-Pro -c ${config} -i ${data_dir} -o ${out_dir}
    elif [[ -f ${hicpro_cmd} ]] 
    then
        # switch to hicpro env
        source ${HOME}/miniconda3/etc/profile.d/conda.sh
        conda activate hicpro
        [[ ${icheck} -gt 0 ]] && \
            ${hicpro_cmd} -c ${config} -i ${data_dir} -o ${out_dir}
        # echo $CONDA_PREFIX
        # echo $(conda env list)
        conda deactivate # exit env
    else
        echo "HiC-Pro not found, the singularity image not exists"
        echo "Please install the HiC-Pro first, see:"
        echo "https://github.com/nservant/HiC-Pro"
        return 1
    fi
}
export -f run_hicpro


## !!3. filter HiC-Pro, validpairs
## Finished in HiC-Pro: finished in HiC-Pro pipeline
## HiC-Pro_output/hic_results/data/smp_name/smp_name.allValidPairs
function filt_hicpro() {
    ## remove dangline pairs and self-circle
    echo 1
}
export -f filt_hicpro 


## !!4. convert to matrix
## Finished in HiC-Pro: finished in HiC-Pro pipeline
## HiC-Pro_output/hic_results/
function build_matrix() {
    # from validPairs to matrix
    local valid_pairs=$1
    local out_dir=$2
    local bin_size=$3
    [[ ! -f ${valid_pairs} ]] && echo "file not exists: ${valid_pairs}" && return 1
    [[ -z $3 ]] && bin_size=25000 # 25kb
    # build matrix
    local build_cmd="${HICPRO_PATH}/scripts/build_matrix"
    local ice_cmd="${HICPRO_PATH}/scripts/ice"
    [[ ! -f ${build_cmd} ]] && echo "`build_matrix`command not found: ${build_cmd}" && return 1
    [[ ! -f ${ice_cmd} ]] && echo "`ice` command not found: ${ice_cmd}" && return 1
    local smp_name="$(basename ${valid_pairs/.allValidPairs})"
    local prefix="${smp_name}_${bin_size}"
    local raw_dir="${out_dir}/${smp_name}/raw/${bin_size}"
    local matrix_raw="${raw_dir}/${prefix}.matrix"
    [[ ! -d ${raw_dir} ]] && mkdir -p ${raw_dir}
    [[ ! -f ${matrix_raw} ]] && \
        ${build_cmd} \
        -–binsize ${bin_size} -–chrsizes ${CHROM_SIZES} \
        -–ifile ${valid_pairs} \
        --matrix-format upper \
        –-oprefix ${out_dir}/${prefix}
    # ICE, nomalized
    local iced_idr="${out_dir}/${smp_name}/iced/${bin_size}"
    local matrix_iced="${iced_dir}/${prefix}_iced.matrix"
    [[ ! -f ${matrix_iced} ]] && \
        python ${ice_cmd} \
        --results_filename ${matrix_iced} \
        –-filter_low_counts_perc 0.05 \
        –-filter_high_counts_perc 0.02 \
        –-max_iter 100 –-eps 0.1 \
        –-remove-all-zeros-loci \
        -–verbose 1 ${matrix_raw}
}
export -f build_matrix


## 5. convert to .hic
## allValidPairs: HiC-Pro_output/hic_results/data/smp_name/smp_name*.allValidPairs
function hicpro_to_juicer() {
    # from validPairs to .hic
    local valid_pairs=$1
    local out_dir=$2
    [[ ! -f ${valid_pairs} ]] && echo "    file not exists: ${valid_pairs}" && return 1
    [[ ! -d ${out_dir} ]] && mkdir -p ${out_dir}
    local hic_file="${out_dir}/$(basename ${valid_pairs}).hic"
    # convert HiC-Pro validpairs to .hic format, supported by Juicer
    # for further analysis using Juicer_tools
    ## run
    local hic2juice="${HICPRO_PATH}/bin/utils/hicpro2juicebox.sh"
    local juicer_tools="${JUICER_PATH}/juicer_tools.jar"
    [[ ! -f ${hic2juice} ]] && echo "command not found: ${hic2juice}" && return 1
    [[ ! -f ${juicer_tools} ]] && echo "error: juicer_tools not exists: ${juicer_tools}" && return 1
    if [[ -f ${hic_file} ]] 
    then 
        echo "    file exists: ${hic_file}" 
    else
        ${hic2juice} -i ${valid_pairs} \
            -o ${out_dir} -j ${juicer_tools} \
            -g ${CHROM_SIZES} -r ${RES_FRAG}
    fi
}
export -f hicpro_to_juicer 


## 6. call TAD - HiCExplorer tools
function find_TAD() {
    local valid_pairs=$1
    local out_dir=$2
    local bin_size=$3
    [[ ! -f ${valid_pairs} ]] && echo "    file not exists: ${valid_pairs}" && return 1
    [[ -z $3 ]] && bin_size=10000 # 10kb
    local smp_name="$(basename ${valid_pairs/.allValidPairs})"
    local sub_dir="${out_dir}/${smp_name}/${bin_size}"
    local prefix="${out_dir}/${smp_name}.${bin_size}"
    local cool_file="${prefix}.cool"
    [[ ! -d ${sub_dir} ]] && mkdir -p ${sub_dir}
    if [[ -f ${cool_file} ]] 
    then
        echo "    file eixists: ${cool_file}"
    else

        conda activate hicexplorer
        cooler cload pairs \
            ${CHROM_SIZES}:${bin_size} "${valid_pairs}"  "${cool_file}" \
            -c1 2 -p1 3 -c2 5 -p2 6
        conda deactivate
            # -c1, chr1 field number
            # -p1, pos1 field number
            # -c2, chr2 field number
            # -p2, pos2 field number
    fi
    # 1. tad
    # local tad=${sub_dir}/${smp_name}.${bin_size}.tad_domains.bed
    local tad="${prefix}.tad_domains.bed"
    local tad_prefix="${prefix}.tad"
    if [[ -f ${tad} ]] 
    then 
        echo "    file exists: ${tad}"
    else
        source ${HOME}/miniconda3/etc/profile.d/conda.sh
        conda activate hicexplorer
        hicFindTADs --matrix ${cool_file} --outPrefix ${tad_prefix} \
            -p ${N_CORES} \
            --thresholdComparisons 0.05 \
            --correctForMultipleTesting fdr
        conda deactivate
    fi
}
export -f find_TAD


## 7. juicer tools
function juicer_arrowhead() {
    # finding contact domains, TADs
    local hic=$1
    local out_dir=$2
    [[ ! -f ${hic} ]] && echo "file not exists: ${hic}" && return 1
    [[ ! -d ${out_dir} ]] && mkdir -p ${out_dir}
    local smp_name=$(basename ${hic/.allValidPairs.hic})
    local juicer_tools="${JUICER_PATH}/juicer_tools.jar"
    local tad_file=""
    ## midium res: -c -m 2000 -r 10000 -k KR
    ## high res:   -c -m 2000 -r 5000 -k KR
    if [[ -f ${tad_file} ]] 
    then 
        echo "file exists: ${tad_file}" 
    else
        # switch to juicer env
        source ${HOME}/miniconda3/etc/profile.d/conda.sh
        conda activate juicer
        java -Xmx20g -jar ${juicer_tools} \
            arrowhead -m 2000 -r 10000 -k KR \
            --threads ${N_CORES} \
            ${hic} ${out_dir}
        conda deactivate
    fi
}
export -f juicer_arrowhead


function juicer_hiccups() {
    # finding chromatin loops
    local hic=$1
    local out_dir=$2
    [[ ! -f ${hic} ]] && echo "file not exists: ${hic}" && return 1
    [[ ! -d ${out_dir} ]] && mkdir -p ${out_dir}
    local smp_name=$(basename ${hic/.allValidPairs.hic})
    local juicer_tools="${JUICER_PATH}/juicer_tools.jar"
    local loop_file=""
    if [[ -f ${loop_file} ]] 
    then 
        echo "file exists: ${loop_file}" 
    else
        # switch to juicer env
        source ${HOME}/miniconda3/etc/profile.d/conda.sh
        conda activate juicer
        java -Xmx20g -jar ${juicer_tools} \
            hiccups --cpu -k KR --threads ${N_CORES} --ignore_sparsity \
            ${hic} ${out_dir} "${smp_name}_loop"
        conda deactivate
    fi
}
export -f juicer_hiccups

## 8. report 
## figures

function run_bathic() {
    # pipeline for BAT Hi-C data
    # see: doi:10.1016/j.ymeth.2019.08.004 for more details
    # in brief
    # 1. trim linker using ChIA-PET2, trimLinker
    # 2. run HiC-Pro 
    # 3. convert to .hic files
    # 4. call TADs
    # 5. call loop, (deeploop?) 
    # 6. report
    local fq1=$1
    local fq2=$2
    local hicpro_config=$3 # required
    local out_dir=$4
    [[ ! -f ${fq1} ]] && echo "fq1 not found: ${fq1}" && return 1
    [[ ! -f ${fq2} ]] && echo "fq2 not found: ${fq2}" && return 1
    [[ ! ${fq1} = *1.fq.gz ]] && "fq1 should be in [*_1.fq.gz] format" && return 1
    [[ ! ${fq2} = *2.fq.gz ]] && "fq2 should be in [*_2.fq.gz] format" && return 1
    local smp_name=$(basename ${fq1/_1.fq.gz}) # !!!
    fq1=$(realpath -s ${fq1})
    fq2=$(realpath -s ${fq2})
    out_dir=$(realpath -s ${out_dir}) # absolute path

    ## 1. copy/link raw_data
    echo "[1/7] - copy raw data"
    local raw_dir="${out_dir}/01.raw_data/${smp_name}"
    [[ ! -d ${raw_dir} ]] && mkdir -p ${raw_dir}
    ln -fs ${fq1} ${raw_dir}/$(basename ${fq1})
    ln -fs ${fq2} ${raw_dir}/$(basename ${fq2})

    ## 2. trim linker
    echo "[2/7] - scan and trim linker"
    local clean_dir="${out_dir}/02.clean_data"
    trim_linker ${fq1} ${fq2} ${clean_dir}

    ## 3. HiC-Pro
    echo "[3/7] - run HiC-Pro"
    local hicpro_dir="${out_dir}/03.hicpro_results"
    run_hicpro ${hicpro_config} ${clean_dir} ${hicpro_dir}
    local valid_pairs="${hicpro_dir}/hic_results/data/${smp_name}/${smp_name}.allValidPairs"

    ## 4. convert to .hic
    echo "[4/7] - convert to .hic"
    local hic_dir="${out_dir}/04.hic_files"
    hicpro_to_juicer ${valid_pairs} ${hic_dir}
    local hic_file="${hic_dir}/${smp_name}.allValidPairs.hic"

    ## 5. call TADs, HiCExplorer
    echo "[5/7] - call TADs, using HiCExplorer"
    local tad_dir="${out_dir}/05.TAD_files"
    for bs in 10000 25000 50000 
    do 
        find_TAD ${valid_pairs} ${tad_dir} ${bs}
    done

    ## 6. call loop, Juicer
    echo "[6/7] - call loop, using Juicer"
    local loop_dir="${out_dir}/06.juicer_files/03.loop"
    juicer_hiccups ${hic_file} ${loop_dir}

    # ## 7. report
    echo "finished!"
}
export -f run_bathic


function usage() {
    cat << EOF
usage : run_bathic -1 demo_1.fq.gz -2 demo_2.fq.gz -o results -c config [-h] [-v]
Use option -h for more information, -v for documentation
EOF
}

function help {
    usage;
    cat << EOF
    -1 FQ1     : read1 of Paired-End file, named as _1.fq.gz 
    -2 FQ2     : read2 of Paired-End file, named as _2.fq.gz 
    -o OUT_DIR : output folder
    -c CONFIG  : config for HiC-Pro, see HiC-Pro/config-hicpro.txt 
    [-h]       : help
    [-v]       : show Documentation
EOF
}
export -f help


function verbose() {
    cat << EOF
Analysis BAT Hi-C data (Ji Xiong Lab, PKU)
see doi:10.1016/j.ymeth.2019.08.004 for more details

## 1. Description:

This pipeline is consist of the following steps
1. trim linker using ChIA-PET2, trimLinker
2. run HiC-Pro 
3. convert to .hic files
4. call TADs
5. call loop, (deeploop?) 
6. report

## 2. Requirements:

For the BAT Hi-C pipeline, the following software and conda env are required

Software:
  - HiC-Pro: v3.1.0, https://github.com/nservant/HiC-Pro 
  - ChIA-PET2: v0.9.3, https://github.com/GuipengLi/ChIA-PET2
  - HiCExplorer: v3.6, https://github.com/deeptools/HiCExplorer
  - Juicer: v1.9.9, jcuda 0.8, https://github.com/aidenlab/juicer
Conda env:
  - hicpro: for HiC-Pro 
  - ChIA-PET2: for ChIA-PET2, trimLinker 
  - hicexplorer: for hicFindTADs, and a bunch of tools for downstream analysis
  - juicer: for juicer tools

## 3. Installation  

### 3.1 Install HiC-Pro

$ git clone https://github.com/nservant/HiC-Pro 
$ cd HiC-Pro
$ sed -i 's/^name:.*/name: hicpro/' environment.yml # rename the env to hicpro
$ conda env create -f environment.yml 
$ conda activate hicpro
$ path-to-hicpro

Generate a GENOME_FRAGMENT file by:
$ HiC-Pro/bin/utils/digest_genome.py -r AC^GT -o mm10_Alul.bed mm10.fa

### 3.2 Install ChIA-PET2

$ git clone https://github.com/GuipengLi/ChIA-PET2
$ cd ChIA-PET2
$ chmod +x bin/ChIA-PET2
$ make
$ realpath bin/trimLinker
path/bin/trimLinker
# record the above path for trimLinker

### 3.3 Install HiCExplorer

$ conda create -n hicexplorer hicexplorer=3.6 python=3.8 -c bioconda -c conda-forge
$ hicFindTADs -h 
usage: hicFindTADs --matrix MATRIX --outPrefix OUTPREFIX --correctForMultipleTesting {fdr,bonferroni,None}
                   [--minDepth INT bp] [--maxDepth INT bp] [--step INT bp]
                   [--TAD_sep_score_prefix TAD_SEP_SCORE_PREFIX] [--thresholdComparisons THRESHOLDCOMPARISONS]
                   [--delta DELTA] [--minBoundaryDistance MINBOUNDARYDISTANCE]
                   [--chromosomes CHROMOSOMES [CHROMOSOMES ...]] [--numberOfProcessors NUMBEROFPROCESSORS] [--help]
                   [--version]

### 3.4 Install Juicer

CPU version:

$ conda create -n juicer python=3.9 openjdk -c bioconda -c conda-forge
$ conda activate juicer

$ git clone https://github.com/aidenlab/juicer
$ mkdir -p <path-to-my-install-dir>
$ cd <path-to-my-install-dir>
$ ln -fs <path-to-juicer>/CPU scripts
$ cd scripts/common
$ wget https://hicfiles.tc4ga.com/public/juicer/juicer_tools.1.9.9_jcuda.0.8.jar
$ ln -s juicer_tools.1.9.9_jcuda.0.8.jar  juicer_tools.jar
$ cd ../..
# record the path to juicer_tools.jar

## 4 Getting Started

+ Prepare fastq files, demo_1.fq.gz, demo_2.fq.gz 
+ Prepare hicpro_config file

$ mkdir work_dir
$ cd work_dir
$ cp <path>/run_bathic.sh .
$ cp <path-to-HiC-Pro>/config-hicpro.txt .

## file: config-hicpro.txt ##
PAIR1_EXT =_1.valid
PAIR2_EXT =_2.valid
BOWTIE2_IDX_PATH = <path to index> 
GENOME_SIZE = <path to chrome.sizes>
GENOME_FRAGMENT = <path to Alu1.bed>
## file: Alu1.bed
$ HiC-Pro/bin/utils/digest_genome.py -r AC^GT -o mm10_Alul.bed mm10.fa

$ conda env list 
ChIA-PET2                /home/wangm/miniconda3/envs/ChIA-PET2
hicexplorer              /home/wangm/miniconda3/envs/hicexplorer
hicpro                *  /home/wangm/miniconda3/envs/hicpro
juicer                   /home/wangm/miniconda3/envs/juicer
$ bash run_bathic.sh demo_1.fq.gz demo_2.fq.gz results config-hicpro.txt
EOF
}
export -f verbose 

## INPUT
[[ $# -lt 1 ]] && usage && exit 1

# Parse command-line options
while getopts ":1:2:c:o:vh" opt
do
    case $opt in
        1) fq1=$OPTARG ;;
        2) fq2=$OPTARG ;;
        c) config=$OPTARG ;;
        o) out_dir=$OPTARG ;;
        v) verbose ;;
        h) help ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            usage
            exit 1
            ;;
        :)
            echo "Option -$OPTARG requires an argument." >&2
            usage
            exit 1
            ;;
    esac
done

if [[ -z ${fq1} || -z ${fq2} || -z ${config} ]]
then
    usage
    exit
fi

# echo $fq1 ${fq2} ${config} ${out_dir}
run_bathic ${fq1} ${fq2} ${config} ${out_dir} 

