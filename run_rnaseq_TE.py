#!/usr/bin/env python
#-*- encoding:utf-8 -*-
"""
Quantification of TE using TEtranscripts/TEcount
Home: https://hammelllab.labsites.cshl.edu/software/#TEtranscripts
code: https://github.com/mhammell-laboratory/TEtranscripts

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
"""

# to-do
# 1. add .lock to make sure each command run only once
# 2. add .


import os
import sys
import re
import argparse
import logging
import subprocess
from multiprocessing import Pool


logging.basicConfig(
    format='[%(asctime)s %(levelname)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    stream=sys.stdout)
log = logging.getLogger(__name__)
log.setLevel('INFO')



def is_valid_aligner(x):
    return x in ['STAR', 'bowtie2', 'bwa', 'hisat2']


def is_STAR_index(x):
    fx = ['Genome', 'SA', 'SAindex']
    x2 = [os.path.join(x, i) for i in fx]
    return all(file_exists(x2))


def is_bowtie2_index(x):
    fx = [1, 2, 3, 4, 'rev.1', 'rev.2']
    x2 = [f'{x}.{i}.bt2' for i in fx]
    return all(file_exists(x2))


def is_hisat2_index(x):
    fx = list(range(1, 9)) # 1 to 8
    x2 = [f'{x}.{i}.ht2' for i in fx]
    return all(file_exists(x2))


def is_bwa_index(x):
    fx = ['amb', 'ann', 'bwt', 'pac', 'sa']
    x2 = [f'{x}.{i}' for i in fx]
    return all(file_exists(x2))


def is_valid_index(x, aligner='STAR'):
    if is_valid_aligner(aligner):
        fx = f'is_{aligner}_index'
        out = eval(fx)(x)
    else:
        out = False
    return out


def get_index(x, aligner='STAR', genome_dir='~/data/genome_db'):
    genome_dir = os.path.expanduser(genome_dir) # default
    genome = get_species(x)
    out = None
    for source in ['GENCODE', 'Ensembl', 'UCSC']:
        # source = get_db_source(x)
        build = get_build(genome, source)
        index = os.path.join(
            genome_dir, build, source, f'{aligner}_index', 'genome'
        )
        if is_valid_index(index, aligner):
            out = index
            break # stop
    return out


def get_gtf(x, genome_dir='~/data/genome_db'):
    genome_dir = os.path.expanduser(genome_dir) # default
    genome = get_species(x)
    # source = get_db_source(genome) # GENCODE/Ensembl/UCSC
    # build = get_build(genome, source) # hg38/GRCh38
    out = None
    for source in ['GENCODE', 'Ensembl', 'UCSC']:
        build = get_build(genome, source)
        f1 = os.path.join(genome_dir, build, source, 'gtf', 'genome.gtf')
        f2 = os.path.join(genome_dir, build, source, 'TE_GTF', 'TE.gtf')
        if all(file_exists([f1, f2])):
            out = [f1, f2]
            break
    return out


def get_species(x):
    sp = {
        'hg38': 'human',
        'human': 'human',
        'mm10': 'mouse',
        'mouse': 'mouse',
        'dm6': 'fruitfly',
        'fruitfly': 'fruitfly'
    }
    # x is str
    return sp.get(x.lower(), None)


def get_build(x, source):
    """
    default:
    human Ensembl GRCh38
    human GENCODE GRCh38
    human UCSC hg38
    mouse Ensembl GRCm38
    mouse GENCODE GRCm38
    mouse UCSC mm10
    fruitfly Ensembl dm6
    fruitfly UCSC dm6
    """
    if get_species(x) == 'human':
        out = 'hg38' if source == 'UCSC' else 'GRCh38'
    elif get_species(x) == 'mouse':
        out = 'mm10' if source == 'UCSC' else 'GRCm38'
    elif get_species(x) == 'fruitfly':
        out = 'dm6'
    else:
        out = None
    return out


def get_db_source(x, db=None):
    # fruitfly: Ensembl > UCSC
    # human/mouse: GENCODE > Ensembl > UCSC
    if get_species(x) == 'fruitfly':
        out = 'Ensembl'
    else:
        out = 'GENCODE' 
    return out


def file_abspath(x):
    if x is None or x == 'None': # in case toml format?!
        out = None
    elif isinstance(x, str):
        out = os.path.abspath(os.path.expanduser(x))
    elif isinstance(x, list):
        out = [file_abspath(i) for i in x]
    else:
        log.warning('x, expect str,list, got {}'.format(type(x).__name__))
        out = x
    return out


def file_prefix(x, fix_pe=False, fix_rep=False):
    if isinstance(x, str):
        # out = file_prefix(x)
        out, ext1 = os.path.splitext(x) # remove extension
        if ext1 in ['.gz', '.bz2']:
            out = os.path.splitext(out)[0] # remove extension
        if fix_pe:
            out = re.sub('[._](r)?[12]$', '', out, flags=re.IGNORECASE)
        if fix_rep:
            out = re.sub('[._](rep|r)[0-9]+$', '', out, flags=re.IGNORECASE)
    elif isinstance(x, list):
        out = [file_prefix(i, fix_pe, fix_rep) for i in x]
    else:
        out = None
    return out


def file_exists(x):
    if x is None:
        out = False
    elif isinstance(x, str):
        out = os.path.exists(x) # file/dir/link
    elif isinstance(x, list):
        out = [file_exists(i) for i in x]
    else:
        log.warning('x, expect str,list, got {}'.format(type(file).__name__))
        out = False
    return out


def check_fx_paired(x, y):
    out = False
    if isinstance(x, str) and isinstance(y, str):
        xname = file_prefix(x, fix_pe=True)
        yname = file_prefix(y, fix_pe=True)
        out = xname == yname
    elif isinstance(x, list) and isinstance(y, list):
        if len(x) == len(y) and len(x) > 0:
            out = all([check_fx_paired(i, j) for i,j in zip(x, y)])
    else:
        pass
    return out


def write_file(s, file, overwrite=False):
    """
    Write str to file
    """
    fdir = os.path.dirname(file)
    if not file_exists(fdir):
        os.makedirs(fdir)
    if file_exists(file) and not overwrite:
        log.info('file exists')
    else:
        with open(file, 'wt') as w:
            w.write(s+'\n')


def check_file(x, **kwargs):
    show_error = kwargs.get('show_error', False)
    show_log = kwargs.get('show_log', False)
    check_empty = kwargs.get('check_empty', False)
    if isinstance(x, str):
        if file_exists(x):
            x_size = os.stat(x).st_size
            # empty gzipped file, size=20
            q_size = 20 if x.endswith('.gz') else 0
            out = x_size > q_size if check_empty else True
            if show_log:
                flag = 'ok' if out else 'failed'
                log.info('{:<6s} : {}'.format(flag, x))
        else:
            if show_error:
                log.error('file not exists: {}'.format(x))
            out = False # failed
    elif isinstance(x, list):
        out = all([check_file(i, **kwargs) for i in x])
    else:
        if show_error:
            log.error('x expect str or list, got {}'.format(type(x).__name__))
        out = False
    return out


def run_shell_cmd():
    p = subprocess.Popen(['/bin/bash','-o','pipefail'], # to catch error in pipe
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
        preexec_fn=os.setsid) # to make a new process with a new PGID
    pid = p.pid
    pgid = os.getpgid(pid)
    cmd_name = '{} ...'.format(os.path.basename(cmd.split()[0]))
    log.info('run_shell_cmd: PID={}, PGID={}, CMD={}'.format(pid, pgid, cmd_name))
    stdout, stderr = p.communicate(cmd)
    rc = p.returncode
    err_str = 'PID={}, PGID={}, RC={}\nSTDERR={}\nSTDOUT={}'.format(
        pid,
        pgid,
        rc,
        stderr.strip(),
        stdout.strip())
    if rc:
        # kill all child processes
        try:
            os.killpg(pgid, signal.SIGKILL)
        except:
            log.error(err_str)
    return (rc, stdout.strip('\n'), stderr.strip('\n'))


def get_bam(x, out_dir):
    if isinstance(x, str):
        xname = file_prefix(x, fix_pe=True)
        out = os.path.join(out_dir, xname, xname+'.bam')
    elif isinstance(x, list):
        out = [get_bam(i, out_dir) for i in x]
    else:
        out = None
    return out


def run_STAR(genome, fq1, fq2=None, out_dir='./', **kwargs):
    """
    Run STAR to mapping reads to reference genome for TEtranscripts
    optimize arguments:
    --outFilterMultimapNmax 100 --winAnchorMultimapNmax 200

    Args:
        genome (str): Name of the genome assembly.
        fq1 (list or str): Path(s) to the BAM file(s) for the treatment condition.
        fq2 (list or str): Path(s) to the BAM file(s) for the control condition.
        out_dir (str): Path to the output directory. Default is current directory.
        overwrite (bool): if True, overwrite existing output files. Default is False.
        stranded (str): Library strandedness. Can be 'no', 'forward' or 'reverse', Default is 'reverse'.

    Returns:
        str: Path to the output count table.

    Raises:
        ValueError: If any input is invalid.    
    """
    threads = kwargs.get('threads', 4)
    overwrite = kwargs.get('overwrite', False)
    # out_dir = kwargs.get('out_dir', './')
    # parse arguments 
    genome = get_species(genome) 
    fq1 = file_abspath(fq1)
    fq2 = file_abspath(fq2)
    out_dir = file_abspath(out_dir)
    star_index = get_index(genome, 'STAR')
    gene_gtf, _ = get_gtf(genome)
    # check_file(fq1, fq2) # paired, file exists # skipped
    prefix = file_prefix(fq1)
    fx_reader = 'zcat' if fq1.endwith('.gz') else '-' # default: '-'
    # output
    star_prefix = os.path.join(out_dir, prefix, f'{prefix}.')
    star_bam = f'{star_prefix}Aligned.out.bam'
    star_stderr = f'{star_prefix}stderr'
    bam = f'{star_prefix}bam'
    cmd_txt = os.path.join(out_dir, prefix, 'cmd.sh')

    # build command line
    cmd = ' '.join([
        f'STAR --runThreadN {threads}',
        f'--genomeDir {star_index}',
        f'--readFilesIn {fq1} {fq2}',
        f'--readFilesCommand {fx_reader}',
        f'--sjdbGTFfile {gene_gtf}',
        f'--outFileNamePrefix {star_prefix}',
        '--outFilterMultimapNmax 100 --winAnchorMultimapNmax 200',
        '--sjdbOverhang 100 --alignSJDBoverhangMin 1',
        '--outFilterMismatchNmax 999 --runRNGseed 777',
        '--outMultimapperOrder Random --outSAMmultNmax 1',
        '--outSAMtype BAM Unsorted --outFilterType BySJout',
        '--alignSJoverhangMin 8 --alignIntronMin 20',
        '--alignIntronMax 1000000 --alignMatesGapMax 1000000',
        f'2>{star_stderr} &&',
        f'ln -fs {os.path.basename(star_bam)} {bam}'
    ])
    write_file(cmd, cmd_txt)
    # run command
    if file_exists(bam) and not overwrite:
        log.info(f'BAM file already exists, skipping run_STAR for {prefix}')
        return(bam)
    else:
        # run_shell_cmd(cmd_txt)
        print('!AAA')
    # output
    return bam 


def run_TEtranscripts(genome, bam_t, bam_c, **kwargs):
    """
    Run TEtranscripts to quantify transcript expression and differential
    expression between two conditions

    Args:
        genome (str): Name of the genome assembly.
        bam_t (list or str): Path(s) to the BAM file(s) for the treatment condition.
        bam_c (list or str): Path(s) to the BAM file(s) for the control condition.
        out_dir (str): Path to the output directory. Default is current directory.
        overwrite (bool): if True, overwrite existing output files. Default is False.
        stranded (str): Library strandedness. Can be 'no', 'forward' or 'reverse', Default is 'reverse'.

    Returns:
        str: Path to the output count table.

    Raises:
        ValueError: If any input is invalid.    
    """
    out_dir = kwargs.get('out_dir', './')
    overwrite = kwargs.get('overwrite', False)
    stranded = kwargs.get('stranded', 'reverse')
    # Check input parameters
    if isinstance(bam_t, str):
        bam_t = [bam_t]
    if isinstance(bam_c, str):
        bam_c = [bam_c]
    for bam_file in bam_t + bam_c:
        if not file_exists(bam_file):
            raise ValueError(f'BAM file not fount: {bam_file}')
        try:
            pysam.AlignmentFile(bam_file, 'rb')
        except ValueError:
            raise ValueError(f'Invalid BAM file: {bam_file}')
    # parse arguments 
    genome = get_species(genome)
    out_dir = file_abspath(out_dir)
    gene_gtf, te_gtf = get_gtf(genome)
    if not all(file_exists(gene_gtf, te_gtf)):
        log.error('[error]: gtf files error')
        return None
    # Set up output files
    name_t = file_prefix(bam_t, rm_rep=True)
    name_c = file_prefix(bam_c, rm_rep=True)
    prefix = f'{name_c[0]}_vs_{name_t[0]}'
    cnt_prefix = os.path.join(out_dir, prefix, f'{prefix}.')
    cnt_table = f'{cnt_prefix}cntTable'
    cnt_stderr = f'{cnt_prefix}stderr'
    cmd_txt = os.path.join(out_dir, prefix, 'cmd.sh')
    bam_t_str = ' '.join(bam_t)
    bam_c_str = ' '.join(bam_c)
    # Build command
    cmd = ' '.join([
        'TEtranscripts --format BAM',
        f'--stranded {stranded}',
        f'-t {bam_t_str} -c {bam_c_str}',
        f'--GTF {gene_gtf}',
        f'--TE {te_gtf}',
        f'--project {prefix} --outdir {out_dir}',
        '--mode multi --minread 1 -i 10 --padj 0.05',
        f'2> {cnt_stderr}'
    ])
    write_file(cmd, cmd_txt)

    # Run command
    if file_exists(cnt_table) and not overwrite:
        log.info(f'File exists, TEtranscripts skipped {cnt_table}')
    else:
        # run_shell_cmd(cmd_txt)
        print('!BBB')
    # Output
    if not file_exists(cnt_table):
        raise ValueError('TEtranscripts output file does not exists.')
    return cnt_table


def run_TEcount(genome, bam, **kwargs):
    # arguments
    out_dir = kwargs.get('out_dir', './') # default
    overwrite = kwargs.get('overwrite', False)
    stranded = kwargs.get('stranded', 'reverse')
    if not isinstance(bam, str) or not file_exists(bam):
        log.error('[error]: bam file error')
        return None
    # parse arguments 
    genome = get_species(genome)
    bam = file_abspath(bam)
    out_dir = file_abspath(out_dir)
    source = get_db_source(genome) # source of genome: Ensembl, GENCODE, UCSC
    gene_gtf, te_gtf = get_gtf(genome, source)
    prefix = file_prefix(bam)
    # output files
    cnt_prefix = os.path.join(out_dir, prefix, f'{prefix}.')
    cnt_table = f'{cnt_prefix}cntTable'
    cnt_stderr = f'{cnt_prefix}stderr'
    cmd_txt = os.path.join(out_dir, prefix, 'cmd.sh')
    # build command
    cmd = ' '.join([
        'TEcount --format BAM',
        f'--stranded {stranded}',
        f'-b {bam}',
        f'--GTF {gene_gtf}',
        f'--TE {te_gtf}',
        f'--project {prefix} --outdir {out_dir}',
        '--mode multi',
        f'2> {cnt_stderr}'
    ])
    write_file(cmd, cmd_txt)
    # run command
    if file_exists(cnt_table) and not overwrite:
        log.info(f'file exists, TEcount skipped ...')
    else:
        # run_shell_cmd(cmd_txt)
        print('!CCC')
    # output
    return cnt_table


def run_rnaseq(genome, t_fq1, t_fq2, c_fq1, c_fq2, **kwargs):
    out_dir = kwargs.get('out_dir', './')
    overwrite = kwargs.get('overwrite', False)
    stranded = kwargs.get('stranded', 'reverse')
    threads = kwargs.get('threads', 4)
    parallel = kwargs.get('parallel', 2)
    # arguments, required, PE?
    if isinstance(t_fq2, list):
        if not check_fx_paired(t_fq1, t_fq2):
            log.error('fq1, fq2 not paired')
            return None
        args_t_fq2 = t_fq2
    else:
        args_t_fq2 = [None] * len(t_fq1)
    if isinstance(c_fq2, list):
        if not check_fx_paired(c_fq1, c_fq2):
            log.error('fq1, fq2 not paired')
            return None
        args_c_fq2 = c_fq2
    else:
        args_c_fq2 = [None] * len(c_fq1)
    # parse arguments 
    genome = get_species(genome)
    out_dir = file_abspath(out_dir)
    t_fq1, t_fq2 = file_abspath([t_fq1, t_fq2])
    c_fq1, c_fq2 = file_abspath([c_fq1, c_fq2])
    gene_gtf, te_gtf = get_gtf(genome)
    if not all(file_exists([gene_gtf, te_gtf])):
        log.error('[error]: gtf files error')
        return None
    # 1. run STAR (multi)
    star_dir = os.path.join(out_dir, '01.align')
    args_genome = [genome] * len(t_fq1 + c_fq1)
    args_fq1 = t_fq1 + c_fq1
    args_fq2 = t_fq2 + c_fq2
    args_out_dir = [star_dir] * len(t_fq1 + c_fq1)
    print(args_genome, args_fq1, args_fq2, args_out_dir)
    if parallel > 1:
        # run in parallel:
        with Pool(processes=parallel) as pool:
            args = zip(args_genome, args_fq1, args_fq2, args_out_dir)
            # pool.map(run_STAR, args)

    # 2. run TEtranscripts
    te_dir = os.path.join(out_dir, '03.de')
    bam_t = get_bam(t_fq1, star_dir)
    bam_c = get_bam(c_fq1, star_dir)
    run_TEtranscripts(genome, bam_t, bam_c, out_dir=te_dir)

    # 3. report
    # to-do: R script


def get_args():
    example = '\n'.join([
        'Example:',
        '1. run program',
        '$ python run_rnaseq_TE.py -g mouse -o results -c1 c_1.fq.gz -c2 c_2.fq.gz -t1 t_1.fq.gz -t2 t_2.fq.gz'
    ])
    parser = argparse.ArgumentParser(
        prog='run_rnaseq_TE', description='run_rnaseq_TE', epilog=example,
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-g', '--genome', default='mouse', required=True,
        help='genome name, eg: mouse, mm10, ... , default [mouse]')
    parser.add_argument('-o', '--out-dir', dest='out_dir', default='./',
        help='directory to save results, default: [./]')
    parser.add_argument('-c1', '--c-fq1', dest='c_fq1', default=None, 
        required=True,
        nargs='+', help='read1 files of control, support multiple files')
    parser.add_argument('-c2', '--c-fq2', dest='c_fq2', default=None, 
        nargs='+', help='read2 files of control, support multiple files')
    parser.add_argument('-t1', '--t-fq1', dest='t_fq1', default=None,
        nargs='+', help='read1 files of treatment, support multiple files')
    parser.add_argument('-t2', '--t-fq2', dest='t_fq2', default=None, 
        nargs='+', help='read2 files of treatment, support multiple files')
    parser.add_argument('-O', '--overwrite', action='store_true',
        help='Overwrite files')
    parser.add_argument('-p', '--threads', default=4, type=int,
        help='Number of threads for each STAR program, default: [4]')
    parser.add_argument('-j', '--parallel', default=2, type=int,
        help='Number of jobs to run in parallel, default: [2]')
    parser.add_argument('-s', '--stranded', choices=['no', 'reverse', 'forward'],
        default='reverse', help=' '.join([
            'stranded library, forward for "first-strand"',
            'library, eg: TruSeq; reverse for "second-strand" library, eg:',
            'QIAseq and dUTP library; no for non-stranded library; default: [no]']))
    return parser


def main():
    args = get_args().parse_args() # input
    kwargs = vars(args).copy() # copy
    # update args
    k_del = ['genome', 'c_fq1', 'c_fq2', 't_fq1', 't_fq2']
    tmp = [kwargs.pop(k) for k in k_del] # remove positional keys
    # run pipeline
    run_rnaseq(
        args.genome, args.c_fq1, args.c_fq2, args.t_fq1, args.t_fq2, **kwargs
    )


if __name__ == '__main__':
    main()
