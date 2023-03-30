# Run Salmon for gene-level quantification
# return
# 1. quant.sf (TPM)
# 2. quant.gene.tsv (TPM)
#
# Date: 2023-03-18


function get_index() {
    gdir="/data/yulab/wangming/data/genome/"
    case $1 in 
        dm6|mm10|hg38)
            echo ${gdir}/${1}/salmon_index/${1}
            ;;
        *)
            echo 1
            ;;
    esac
}
export -f get_index


function tx2gene_R() {
    # convert transcript level to gene level
    # using scripts
    cat << EOF
## quant gene-level TPM
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1){
  print("Usage: Rscript run_salmon.R <quant.sf>")
  print("")
  print("Option:")
  print("  quant.sf     The output of salmon")
  stop("arguments failed")
}
sf <- args[1]

# locate the info.json
info <- file.path(dirname(sf), "cmd_info.json")
da   <- jsonlite::read_json(info)
idx  <- da\$index
t2g_f <- file.path(idx, "tx2gene.csv")

if(! file.exists(t2g_f)) {
  warning(glue::glue("t2g file not exists: {t2g_f}"))
} else {
  suppressPackageStartupMessages(library(tximport))
  tx2gene <- readr::read_csv(t2g_f, show_col_types = FALSE)
  txi <- tximport::tximport(sf, type = "salmon", tx2gene = tx2gene,
                            abundanceCol = "TPM")
  # gene, length, count, tpm
  df1 <- cbind(as.data.frame(txi\$length),
               as.data.frame(txi\$abundance),
               as.data.frame(txi\$counts))
  df1 <- round(df1, 4)
  colnames(df1) <- c("length", "TPM", "count")
  df1 <- tibble::rownames_to_column(df1, "id")
  # output
  gene_tpm <- file.path(dirname(sf), "quant.gene.tsv")
  message(glue::glue("save gene counts to file: {gene_tpm}"))
  readr::write_tsv(df1, gene_tpm, col_names = TRUE)
}
EOF
}
export -f tx2gene_R

# prepare working directory
# genome, out_dir, fq1, fq2
function run_salmon() {
    genome=$1
    out_dir=$2
    fq1=$3
    fq2=$4
    # for single-end
    [[ ${fq1} == ${fq2} ]] && fq2="" #
    # check index
    idx=$(get_index ${genome})
    [[ ! -d ${idx} ]] && echo "index not exists: ${idx}" && return 1
    prefix="$(basename ${fq1/.gz})"
    prefix=$(echo ${prefix} | sed -Ee 's/(_1)?.f(ast)?q(.gz)?//')
    out_dir="${out_dir}/${prefix}"
    [[ ! -d ${out_dir} ]] && mkdir -p ${out_dir}
    out_dir=$(realpath -s ${out_dir})
    # check exists
    sf="${out_dir}/quant.sf"
    log="${out_dir}/log.stdout"
    # 1. run salmon
    if [[ -f ${sf} ]] ; then
        echo "file exists: ${sf}"
    else
        if [[ -f ${fq2} ]] ; then
            # Paired-end mode
            salmon quant --gcBias --validateMappings -l A -p 8 -i ${idx} -o ${out_dir} -1 ${fq1} -2 ${fq2} 2> ${log}
        else
            # Single-end mode
            salmon quant --gcBias --validateMappings -l A -p 8 -i ${idx} -o ${out_dir} -r ${fq1} 2> ${log}
        fi
    fi
    # 2. run tx2gene 
    run_t2g="${out_dir}/run_tx2gene.sh"
    t2g_R="${out_dir}/tx2gene.R"
    gg="${out_dir}/quant.gene.tsv"
    if [[ -f ${gg} ]] ; then
        echo "file exists: ${gg}"
    else
        tx2gene_R > ${t2g_R}
        # echo "Rscript $(basename ${t2g_R}) $(basename ${sf})" > ${run_t2g}
        echo "Rscript $(basename ${t2g_R}) $(basename ${sf})" > ${run_t2g}
        pre_dir=$(pwd)
        cd ${out_dir}
        bash ${run_t2g}
        cd ${pre_dir}
    fi
}
export -f run_salmon


function run_salmon2() {
    fq1=$1
    out_dir=$2
    fq2=${fq1/_1.fq/_2.fq}
    [[ ${fq1} == ${fq2} ]] && fq2="" #
    run_salmon hg38 ${out_dir} ${fq1} ${fq2}
}
export -f run_salmon2 

# parallel -j 4 run_salmon2 {} results/RNAseq_salmon ::: data/clean_data/*/*1.fq.gz 

[[ $# -lt 3 ]] && echo "Usage: run_salmon.sh <genome> <out_dir> <fq1> <fq2>" && exit 1
run_salmon $1 $2 $3 $4
