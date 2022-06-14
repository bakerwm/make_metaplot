# convert bam to fastq (PE)
#
# exclude -F 2048 supplementary, -F 256 secondary
# include -f 2 mapped in proper pair

function bam2fq(){
    local bam=$1
    local outdir=$2
    bname=$(basename ${bam/.bam})
    # outdir=$(dirname ${x})
    fq1="${outdir}/${bname}_1.fq"
    fq2="${outdir}/${bname}_2.fq"
    if [[ -f ${bam} ]] 
    then
        [[ ! -d ${outdir} ]] && mkdir -p ${outdir}
        echo ${bname}
        samtools fastq -@ 8 -1 ${fq1} -2 ${fq2} -F 2304 -f 2 ${bam}
    fi
}
export -f bam2fq

[[ $# -lt 1 ]] && echo "bam2fq.sh <outdir> bam1 [bam2 ...]" && exit 1
outdir=$1
bam=${@}
bam=(${bam[@]/${outdir}}) # update
# echo ${bam[@]}
parallel --jobs 4 bam2fq {} ${outdir} ::: ${bam[@]}
