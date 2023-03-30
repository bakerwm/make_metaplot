#!/usr/bin/bash 

# Date: 2023-03-21

function example() {
    cat <<EOF
# 1. run STAR
## key: --outFilterMultimapNmax 100 --winAnchorMultimapNmax 200
$ STAR --runThreadN \${threads} --genomeDir Index/genome \\
    --readFilesIn \${fq1} \${fq2} --readFilesCommand zcat \\
    --sjdbGTFfile \${gtf} --outFileNamePrefix \${prefix} \\
    --sjdbOverhang 100 --alignSJDBoverhangMin 1 \\
    --outFilterMismatchNmax 999 --runRNGseed 777 \\
    --outFilterMultimapNmax 100 --winAnchorMultimapNmax 200 \\
    --outMultimapperOrder Random --outSAMmultNmax 1 \\
    --outSAMtype BAM --outFilterType BySJout --alignSJoverhangMin 8 \\
    --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000

# 2. run TEtranscripts
## required: treatment, control, 2-reps for each
$ TEtranscripts --format BAM --stranded reverse \\
    -t piwiKD_ Aligned.out.bam \\
    -c control_Aligned.out.bam \\
    --GTF \${gene.gtf} --TE \${te.gtf} \\
    --mode multi --minread 1 -i 10 --padj 0.05 \\
    --project \${project_name} --outdir \${outdir}

# 3. run TEcount
## required: for single BAM
$ TEcount --format BAM --stranded reverse \\
    -b \${bam} --GTF \${gene.gtf} --TE \${te.gtf} \\
    --mode multi --minread 1 -i 10 \\
    --project \${project_name} --outdir \${outdir}

# see: https://github.com/mhammell-laboratory/TEtranscripts

EOF
}
export -f example


# # Description
# function usage() {
#     cat <<'EOF'
# Usage: bash run_rnaseq_TE.sh <1> <2> <3> ...

# Quantification for Gene + Transposon Element (TE) using TEtranscripts/TEcoun.

# EOF
# }
# export -f usage


# for aligner index
# species: human, mouse, fruitfly
# source: Ensembl, GENCODE, USSC
# build: latest (GRCh38, GRCm38, dm6, hg38, mm10) # mm39 not yet
# aligner: STAR (v2.7.2a), bowtie2 (v2.4.5), bwa (v0.7.17-r1188), hisat2 (v2.1.0)
function db_table() {
    cat << EOF
human Ensembl GRCh38
human GENCODE GRCh38
human UCSC hg38
mouse Ensembl GRCm38
mouse GENCODE GRCm38
mouse UCSC mm10
fruitfly Ensembl dm6
fruitfly UCSC dm6
EOF
}
export -f db_table


function get_species() {
    case $1 in 
        hg38|human|Human)
            echo human
            ;;
        mm10|mouse|Mouse)
            echo mouse
            ;;
        dm6|fruitfly|Fruitfly)
            echo fruitfly
            ;;
        *)
            echo null
    esac
}
export -f get_species


function get_source() {
    so=${1,,} # to lowercase
    case ${so} in 
        ensembl) 
            echo Ensembl
            ;;
        gencode)
            echo GENCODE
            ;;
        ucsc)
            echo UCSC
            ;;
        *)
            echo null
            ;;
    esac
}
export -f get_source 


function get_build() {
    # arguments: species, source
    species=$(get_species $1)
    source=$(get_source $2)
    # echo "!!!C ${species} ${source}"
    db_table | awk -v sp=${species} -v so=${source} '$1==sp && $2 == so {print $3}'
}
export -f get_build


# 0=yes, 1=no
function is_valid_aligner() {
    tag=1 # 1=no, 0=yes
    for a in STAR bowtie2 bwa hisat2; do
        [[ $1 == ${a} ]] && tag=0 && break
    done
    echo ${tag}
}
export -f is_valid_aligner


# 0=yes, 1=no
function is_STAR_idx() {
    # arguments: idx
    # 0=yes, 1=no
    tag=0
    for i in Genome SA SAindex; do
        [[ ! -f $1/${i} ]] && tag=1 && break
    done
    echo ${tag}
}
export -f is_STAR_idx


# 0=yes, 1=no
function is_bowtie2_idx() {
    # arguments: idx
    # 0=yes, 1=no
    tag=0
    for i in $(seq 1 4) rev.1 rev.2; do
        [[ ! -f $1/${i}.bt2 ]] && tag=1 && break
    done
    echo ${tag}
}
export -f is_bowtie2_idx


# 0=yes, 1=no
function is_hisat2_idx() {
    # arguments: idx
    # 0=yes, 1=no
    tag=0
    for i in $(seq 1 8) ; do
        [[ ! -f $1/${i}.ht2 ]] && tag=1 && break
    done
    echo ${tag}
}
export -f is_hisat2_idx


# 0=yes, 1=no
function is_bwa_idx() {
    # arguments: idx
    # 0=yes, 1=no
    tag=0
    for i in amb ann bwt pac sa; do
        [[ ! -f $1/${i}.bt2 ]] && tag=1 && break
    done
    echo ${tag}
}
export -f is_bowtie2_idx


# 0=yes, 1=no
function is_valid_idx() {
    # arguments: aligner idx
    # 0=yes, 1=no
    case $1 in 
        STAR)
            is_STAR_idx $2
            ;;
        bowtie2)
            is_bowtie2_idx $2
            ;;
        hisat2)
            is_hisat2_idx $2
            ;;
        bwa)
            is_bwa_idx $2
            ;;
        *)
            echo 1 # no
    esac
}
export -f is_valid_idx


function get_index() {
    # argument: species, source, aligner
    genome_dir="${HOME}/data/genome_db" # custome_db
    species=$(get_species $1) # human mouse fruitfly
    source=$(get_source $2)   # Ensembl GENCODE UCSC
    build=$(get_build $1 $2)  # GRCh38 hg38
    [[ $(is_valid_aligner $3) -gt 0 ]] && aligner=null || aligner=$3
    # construct index
    idx="${genome_dir}/${build}/${source}/${aligner}_index/genome"
    [[ $(is_valid_idx ${aligner} ${idx}) -gt 0 ]] && echo "No index found: $1 $2 $3 ${idx}" && return 1
    echo ${idx}
}
export -f get_index


function get_gtf() {
    # argument: species, source
    genome_dir="${HOME}/data/genome_db" # custome_db
    species=$(get_species $1) # human mouse fruitfly
    source=$(get_source $2)   # Ensembl GENCODE UCSC
    build=$(get_build $1 $2)  # GRCh38 hg38
    gene_gtf="${genome_dir}/${build}/${source}/gtf/genome.gtf"
    te_gtf="${genome_dir}/${build}/${source}/TE_GTF/TE.gtf"
    # check file exists
    tag=0 # 
    [[ ! -f ${gene_gtf} ]] && gene_gtf="-" && tag=1
    [[ ! -f ${te_gtf} ]] && te_gtf="-" && tag=1
    # echo "GENE_GTF: ${gene_gtf}"
    # echo "TE_GTF: ${te_gtf}"
    [[ ${tag} -gt 0 ]] && echo "GTF files not exists" && return 1
    echo ${gene_gtf} ${te_gtf}
}
export -f get_gtf


function fx_prefix() {
    # extract prefix of fastq name
    echo $@ | xargs -n 1 basename | sed -Ee 's/\.\w+$//' -Ee 's/(_[0-9])?.f(ast)?[aq]$//i' | xargs 
    # basename $1 |  sed -Ee 's/\.\w+$//' -Ee 's/(_[0-9])?.f(ast)?[aq]$//i'
    # basename $1 | sed -E 's/(_[0-9])?.f(ast)?[aq](.gz)?//i'
}
export -f fx_prefix


function file_abspath() {
    # convert bam to absolute path
    # require "" in arguments
    echo $@ | xargs -n 1 realpath -s | xargs
}
export -f file_abspath


function file_exists_log() {
    for b in $@ ; do 
        [[ -f ${b} ]] && t="ok" || t="failed"
        printf "%8s : ${b}\n" ${t}
    done
}
export -f file_exists_log 


function file_exists() {
    # file exists: 0=yes, 1=not
    tag=0
    for b in $@ ; do 
        [[ -f ${b} ]] && tag=1 && break
    done
    return ${tag}
}
export -f file_exists


# STAR, for TEtranscripts
function run_STAR() {
    # arguments: genome, out_dir, fq1, fq2
    local threads=24 # default
    local genome=$(get_species $1) #
    local out_dir=$2 #
    [[ ${genome} == "fruitfly" ]] && source="Ensembl" || source="GENCODE" # !!! default
    local fq1=$(file_abspath $3) # read1 of PE
    local fq2=$(file_abspath $4) # read2 of PE
    local STAR_idx=$(get_index ${genome} ${source} STAR)
    local gtf=$(get_gtf ${genome} ${source})
    [[ ${gtf} == 1 ]] && echo "GTF not found: ${genome} ${source}" && return 1
    local gene_gtf=$(echo ${gtf} | awk '{print $1}')
    # local te_gtf=$(echo ${gtf} | awk '{print $2}')
    # check arguments
    [[ ${genome} == null ]] && echo "unknown genome: $1" && return 1
    [[ ! -f ${fq1} ]] && echo "fq1 not exists: $3" && return 1
    [[ ! -f ${fq2} ]] && fq2="" # skip read2
    # output
    local prefix=$(fx_prefix ${fq1})
    out_dir=${out_dir}/${prefix} # !!! update output
    [[ ! -d ${out_dir} ]] && mkdir -p ${out_dir}
    out_dir=$(file_abspath ${out_dir}) # absolute path
    local STAR_prefix="${out_dir}/${prefix}." # prefix.Aligned.out.bam
    local STAR_bam=${STAR_prefix}Aligned.out.bam
    local STAR_stderr=${STAR_prefix}stderr
    local bam=${out_dir}/${prefix}.bam # final output
    local cmd="${out_dir}/cmd.sh"
    # run STAR
    [[ -f ${bam} ]] && echo "${bam}" && return 1
    # command line
    local cmd_txt=$(cat <<EOF
STAR --runThreadN ${threads} --genomeDir ${STAR_idx} \
    --readFilesIn ${fq1} ${fq2} --readFilesCommand zcat \
    --sjdbGTFfile ${gene_gtf} --outFileNamePrefix ${STAR_prefix} \
    --sjdbOverhang 100 --alignSJDBoverhangMin 1 \
    --outFilterMismatchNmax 999 --runRNGseed 777 \
    --outFilterMultimapNmax 100 --winAnchorMultimapNmax 200 \
    --outMultimapperOrder Random --outSAMmultNmax 1 \
    --outSAMtype BAM Unsorted --outFilterType BySJout --alignSJoverhangMin 8 \
    --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 \
    2> ${STAR_stderr} && \
    ln -fs $(basename ${STAR_bam}) ${bam} # create symlink
EOF
)
    echo ${cmd_txt} | tee ${cmd} | bash
    # output    
    echo ${bam}
}
export -f run_STAR


function run_TEcount() {
    # arguments: genome, out_dir, bam
    # default: GENCODE for human,mouse; Ensembl for fruitfly
    local genome=$(get_species $1) #
    local out_dir=$2
    [[ ${genome} == "Fruitfly" ]] && source="Ensembl" || source="GENCODE" # !!! default
    local gtf=$(get_gtf ${genome} ${source})
    [[ ${gtf} == 1 ]] && echo "GTF not found: ${genome} ${source}" && return 1
    local gene_gtf=$(echo ${gtf} | awk '{print $1}')
    local te_gtf=$(echo ${gtf} | awk '{print $2}')
    # check bam files
    local bam=$(file_abspath "$3") # absolute path    
    file_exists_log ${bam} # show log
    [[ $(file_exists ${bam}) -gt 0 ]] && echo "BAM file failed" && return 1
    local prefix=$(fx_prefix ${bam}) #
    out_dir=${out_dir}/${prefix} # !!! update out_dir
    [[ ! -d ${out_dir} ]] && mkdir -p ${out_dir}
    out_dir=$(file_abspath ${out_dir}) # absolute path
    local cmd=${out_dir}/cmd.sh
    # command line
    local cmd_txt=$(cat <<EOF
TEcount -b ${bam} --GTF ${gene_gtf} --TE ${te_gtf} \
    --format BAM --stranded reverse --mode multi -i 10 \
    --project ${prefix} --outdir ${out_dir}
EOF
)
    # prepare command file
    echo ${cmd_txt} | tee ${cmd} | bash
    echo ${cmd}
}
export -f run_TEcount


function run_TEtranscripts() {
    # arguments: genome, out_dir, "bam_t", "bam_c"
    # default: GENCODE for human,mouse; Ensembl for fruitfly
    # support multiple bam files
    local genome=$(get_species $1) #
    local out_dir=$2
    [[ ${genome} == "fruitfly" ]] && source="Ensembl" || source="GENCODE"
    local gtf=$(get_gtf ${genome} ${source})
    [[ ${gtf} == 1 ]] && echo "GTF not found: ${genome} ${source}" && return 1
    local gene_gtf=$(echo ${gtf} | awk '{print $1}')
    local te_gtf=$(echo ${gtf} | awk '{print $2}')
    # check bam files
    local bam_t=$(file_abspath "$3")
    local bam_c=$(file_abspath "$4")
    local bname_t=($(fx_prefix "${bam_t}"))
    bname_c=($(fx_prefix "${bam_c}"))
    bname_t1=$(echo ${bname_t[0]} | sed -Ee 's/_r(ep)?[0-9]+$//i') # treat name
    bname_c1=$(echo ${bname_c[0]} | sed -Ee 's/_r(ep)?[0-9]+$//i') # ctrl name
    # bam exists
    file_exists_log ${bam_t} ${bam_c} # show log
    [[ $(file_exists ${bam_t} ${bam_c}) -gt 0 ]] && echo "BAM file failed" && return 1
    local prefix="${bname_c1}_vs_${bname_t1}"
    out_dir=${out_dir}/${prefix} # !!! update out_dir
    [[ ! -d ${out_dir} ]] && mkdir -p ${out_dir}
    out_dir=$(file_abspath ${out_dir})
    local cnt_table=${out_dir}/${prefix}.cntTable
    [[ -f ${cnt_table} ]] && echo ${cnt_table} && return 1
    local cmd=${out_dir}/cmd.sh
    # command line
    local cmd_txt=$(cat <<EOF
TEtranscripts --format BAM --stranded reverse \
    -t ${bam_t} -c ${bam_c} \
    --GTF ${gene_gtf} --TE ${te_gtf} \
    --mode multi --minread 1 -i 10 --padj 0.05 \
    --project ${prefix} --outdir ${out_dir}
EOF
)
    # prepare command file
    echo ${cmd_txt} | tee ${cmd} | bash
    echo ${cnt_table}
}
export -f run_TEtranscripts


function run_rnaseq() {
    # arguments: genome, out_dir, fq1, fq2
    local genome=$(get_species $1) #
    local out_dir=$2
    local ctrl_fq1=$(file_abspath $3)  # read1 of wildtype
    local ctrl_fq2=$(file_abspath $4)  # read2 of wildtype
    local treat_fq1=$(file_abspath $5) # read1 of mutant
    local treat_fq2=$(file_abspath $6) # read2 of mutant
    [[ ! -d ${out_dir} ]] && mkdir -p ${out_dir}
    out_dir=$(file_abspath ${out_dir})
    # 1. run STAR
    local STAR_dir=${out_dir}/01.align_STAR
    local TEtool_dir=${out_dir}/02.quant
    local bam_ctrl=()
    for fq1 in ${ctrl_fq1}; do 
        fq2=$(echo ${fq1} | sed -E 's/_1.(f(ast)?q)/_2.\1/g')
        [[ ! -f ${fq2} ]] && fq2=""
        bam=$(run_STAR ${genome} ${STAR_dir} ${fq1} ${fq2})
        # run_TEcount ${genome} ${TEtool_dir} ${bam}
        bam_ctrl+=("${bam}")
    done
    # echo "!AAA" ${STAR_dir}
    bam_treat=()
    for fq1 in ${treat_fq1}; do
        fq2=$(echo ${fq1} | sed -E 's/_1.(f(ast)?q)/_2.\1/g')
        [[ ! -f ${fq2} ]] && fq2=""
        bam=$(run_STAR ${genome} ${STAR_dir} ${fq1} ${fq2})
        # run_TEcount ${genome} ${TEtool_dir} ${bam}
        bam_treat+=("${bam}")
    done
    # 2. run TEtranscripts
    de_dir=${out_dir}/03.deseq
    local bam_t=$(echo ${bam_treat[@]})
    local bam_c=$(echo ${bam_ctrl[@]})
    cnt_tab=$(run_TEtranscripts ${genome} ${de_dir} "${bam_t}" "${bam_c}")
    echo ${cnt_tab}
}
export -f run_rnaseq


function usage() {
  echo "Usage: $0 [-g <genome>] [-o <out_dir>] [-w <ctrl_fq1>] [-W <ctrl_fq2>] [-m <treat_fq1>] [-M <treat_fq2>]"
  echo "    -g, --genome    Specify genome (default: fruitfly)"
  echo "    -o, --out_dir   Specify output diretory (default: ./)"
  echo "    -c, --ctrl_fq1    Specify read1 of control sample, rep1,rep2,... for multiple files (default: null)"
  # echo "    -C, --ctrl_fq2    Specify read2 of control sample, rep1,rep2,... for multiple files (default: null)"
  echo "    -t, --treat_fq1   Specify read1 of treatment sample, rep1,rep2,... for multiple files (default: null)"
  # echo "    -T, --treat_fq2   Specify read2 of treatment sample, rep1,rep2,... for multiple files (default: null)"
  echo "Note: "
  echo "    Guess fq2 by fq1 name, _1.fq to _2.fq"
  exit 1
}


function main() {
    # Default values
    local genome="fruitfly"
    local out_dir="./"
    local ctrl_fq1="null"
    local ctrl_fq2="null"
    local treat_fq1="null"
    local treat_fq2="null"
    local threads=4
    # Parse command-line arguments
    opts=$(getopt -o g:o:c:C:t:T:p: --long genome:,out_dir:,ctrl_fq1:,ctrl_fq2:,treat_fq1:,treat_fq2:threads: -- "$@")
    eval set -- "$opts"
    while true ; do 
        case "$1" in 
            -g|--genome)
                genome="$2"
                shift 2
                ;;
            -o|--out_dir)
                out_dir="$2"
                shift 2
                ;;
            -c|--ctrl_fq1)
                ctrl_fq1="$2"
                shift 2
                ;;
            -C|--ctrl_fq2)
                ctrl_fq2="$2"
                shift 2
                ;;
            -t|--treat_fq1)
                treat_fq1="$2"
                shift 2
                ;;
            -T|--treat_fq2)
                treat_fq2="$2"
                shift 2
                ;;
            -p|--threads)
                threads="$2"
                shift 2
                ;;
            --)
                shift
                break
                ;;
        esac
    done
    # usage    
    # message
    if [ "${ctrl_fq1}" == "null" ] || [ "${treat_fq1}" == "null" ]; then
        usage
    fi
    # prepare fq1
    ctrl_fq1=$(echo ${ctrl_fq1} | sed 's/,/ /g')
    treat_fq1=$(echo ${treat_fq1} | sed 's/,/ /g')
    # guess fq2 (by fq1)
    ctrl_fq2=$(echo ${ctrl_fq1} | sed -E 's/_1.(f(ast)?q)/_2.\1/g')
    treat_fq2=$(echo ${treat_fq1} | sed -E 's/_1.(f(ast)?q)/_2.\1/g')
    # check file exists
    echo "control fastq files:"
    file_exists_log ${ctrl_fq1} ${ctrl_fq2} # show log
    echo "treatment fastq files:"
    file_exists_log ${treat_fq1} ${treat_fq2} # show log
    [[ $(file_exists ${ctrl_fq1} ${treat_fq1}) -gt 0 ]] && echo "fastq files failed ..." && return 1 
    [[ $(file_exists ${ctrl_fq2}) -gt 0 ]] && ctrl_fq2=null
    [[ $(file_exists ${treat_fq2}) -gt 0 ]] && treat_fq2=null
    run_rnaseq ${genome} ${out_dir} "${ctrl_fq1}" "${ctrl_fq2}" "${treat_fq1}" "${treat_fq2}"
}


main $@
