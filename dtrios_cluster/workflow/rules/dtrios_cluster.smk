shell.prefix("source /data/antwerpen/grp/asvardal/share/hscon5_setup.sh; ")
import os, math, shutil
import pandas as pd



wildcard_constraints:
    job = "\d+",
    split = "\d+",
    filebase = "[^/]+",
    outtype = "(BBAA)|(Dmin)|"

localrules:
    aggregate_per_job_per_split,
    check_final_output,
    concatenate,
    create_trio_splits

try:
    vcfs = [config['vcf']]
    assert 'vcfs' not in config.keys(), "config can only contain 'vcf' or 'vcfs', not both"
    logger.warning('Mutltiple jobs are only implemented across vcfs.'
                   ' For optimal perfomance provide as many vcfs as you want jobs to run (e.g., one per chromosome). '
                   'Each job will be parallelised across config["n_threads"]. ')
except KeyError:
    vcfs = config['vcfs']

dsuite_path = config['dsuite_path']
remove_intermediate = config['remove_intermediate']

sets = pd.read_csv(f"../config/{analysis}/sets.txt", sep='\t', header=None)
n_populations = len(set([s for s in sets.iloc[:,1].values if s not in ['xxx', 'Outgroup']]))
assert 'Outgroup' in sets.iloc[:,1].values, "Outgroup not in sets: ".format(sets)


n_splits = config['n_splits']
combinations = math.comb(n_populations, 3)
lines_per_split = math.ceil(combinations/n_splits)
lines_per_split = math.ceil(combinations/n_splits)
n_splits = math.ceil(combinations/lines_per_split)
suffix_length = math.ceil(math.log10(n_splits))

logger.info(f"There are {n_populations} populations present in set.txt. Hence there are {combinations} trios.")
logger.info(f"Trio output is going to be split into {n_splits} chunks for intermediate processing.")


def get_vcf(wildcards):
    return vcfs[int(wildcards.job)]

checkpoint dtrios_parallel_no_combine:
    input:
         vcf = get_vcf,
         sets = "../config/{analysis}/sets.txt"
    output:
        out_dir = directory("../output/{analysis}/DTparallel/job{job}")
    params:
        dtrios_parallel_path = os.path.abspath(os.path.join(config['dsuite_path'],
                                                        "../utils/DtriosParallel")),
        dsuite_path = config['dsuite_path']
    resources:
        walltime = 72,
        mem_mb = 120000
        #walltime = 336:00:00
    threads: config['cores_per_node']
    shell:
        """
         module load GCC; \
         ../../scripts/DtriosParallel --JKnum 1 --no-combine --dsuite-path {params.dsuite_path} \
         --intermediate-prefix {output.out_dir}/DT \
         --keep-intermediate --cores {threads} {input.sets} {input.vcf};
        """

rule create_trio_splits:
    input:
        combine = "../output/{analysis}/DTparallel/job{job}/{filebase}_combine.txt",
        stderr = "../output/{analysis}/DTparallel/job{job}/{filebase}_combine_stderr.txt"
    params:
        lines_per_split = lines_per_split,
        suffix_length = suffix_length,
        outdir = "../output/{analysis}/DTparallel/split_job{job}",
    output:
        combines = [("../output/{{analysis}}/DTparallel/split_job{{job}}/{{filebase}}_{:0" + str(suffix_length)
                    + "d}_combine.txt").format(i) for i in range(n_splits)],
        stderrs = [("../output/{{analysis}}/DTparallel/split_job{{job}}/{{filebase}}_{:0" + str(suffix_length)
                    + "d}_combine_stderr.txt").format(i) for i in range(n_splits)],
    resources:
        walltime=72,
        mem_mb=5000
    shell:
        """
        mkdir -p {params.outdir}
        split --additional-suffix=_combine.txt -a {params.suffix_length} \
          -d --lines {params.lines_per_split} {input.combine} {params.outdir}/{wildcards.filebase}_ ; \
        split --additional-suffix=_combine_stderr.txt -a {params.suffix_length} \
          -d --lines {params.lines_per_split} {input.stderr} {params.outdir}/{wildcards.filebase}_ ;
        """


def get_combine_input(wildcards, stderr=False):
    job = wildcards.job
    if stderr:
        sfx = "_stderr"
    else:
        sfx = ""

    dtp_files = []
    split = wildcards.split
    analysis = wildcards.analysis

    dirname = checkpoints.dtrios_parallel_no_combine.get(analysis=wildcards.analysis, job=job).output[0]

    combine_filebases = glob_wildcards(os.path.join(dirname,"{filebase,[^/]+}_combine.txt")).filebase

    error_str = \
    """
    Number of combine files should match number of cores per node
    Please check/ delete output folder of 
    dtrios_parallel_no_combine and rerun. There are {} cores_per_node, but filebases are: {}."
    """
    assert len(combine_filebases) == config['cores_per_node'], error_str.format(config['cores_per_node'],
                                                                                combine_filebases)
    fns = []
    for fb in combine_filebases:
        fn = os.path.join(f"../output/{analysis}/DTparallel/split_job{job}", ("{}_{:0" + str(suffix_length) + "d}_combine{}.txt").format(fb,
                            int(split), sfx))
        fns.append(fn)

    return fns


rule aggregate_per_job_per_split:
    input:
        combines = get_combine_input,
        stderrs = lambda wildcards: get_combine_input(wildcards,stderr=True)
    output:
        agg="../output/{analysis}/DTparallel/{job}_{split}.aggregated"
    run:
        with open(output.agg,'w') as f:
            for fn in input.combines:
                base = fn.rsplit('_',1)[0]
                f.write(base + '\n')


rule combine:
    input:
        aggregated = expand("../output/{{analysis}}/DTparallel/{job}_{{split}}.aggregated", job=range(len(vcfs)))
    output:
        bbaa = "../output/{analysis}/combined/{split}/{analysis}_combined_BBAA.txt",
        dmin = "../output/{analysis}/combined/{split}/{analysis}_combined_Dmin.txt",
        out_dir = directory("../output/{analysis}/DTparallel/combined/{split}/")
    params:
        dsuite_path = config['dsuite_path'],
        #combines_bases = lambda wildcards, input: [fn.rsplit('_',1)[0] for fn in input.combines],
        out_prefix = lambda wildcards, output: output.bbaa.rsplit('_',2)[0],
    resources:
        walltime=72,
        mem_mb=5000
    run:
        out_prefix = params.out_prefix
        dsuite_path = params.dsuite_path
        combines_bases = []

        for fn in input.aggregated:
            with open(fn) as f:
                for line in f.readlines():
                    combines_bases.append(line.strip())
        combines_bases = ' '.join(combines_bases)
        os.makedirs(output.out_dir, exist_ok=True)
        shell("module load GCC; {dsuite_path}/Dsuite DtriosCombine --out-prefix {out_prefix} {combines_bases}")


rule concatenate:
    input:
        bbaas = expand("../output/{{analysis}}/combined/{split}/{{analysis}}_combined_{{outtype}}.txt",
                    split=[("{:0" + str(suffix_length) + "d}").format(s) for s in range(n_splits)])
    output:
        bbaa = "../output/{analysis}/{analysis}_{outtype}.txt",
    shell:
        """
        head -n1 {input.bbaas[0]} > {output.bbaa}; \
        awk 'FNR!=1' {input.bbaas} >> {output.bbaa}; \
        """
import subprocess

rule check_final_output:
    input:
        bbaa = "../output/{analysis}/{analysis}_BBAA.txt",
        dmin= "../output/{analysis}/{analysis}_Dmin.txt",
    output:
        check = "../output/{analysis}/BBAA_check_passed.txt"
    params:
        combinations = combinations
    run:
        analysis = wildcards.analysis
        bbaa = input.bbaa
        p = subprocess.Popen(f'wc -l {bbaa}', shell=True, stdout=subprocess.PIPE)
        o, e = p.communicate()
        lines = int(o.decode().split(' ')[0]) - 1
        if lines != params.combinations:
            logger.error(f"There are {n_populations} populations so we expect {combinations} trios. "
                         f"However, there are {lines} trios present in {bbaa}")
        else:
            outstr = (f"There are {n_populations} populations so we expect {combinations} trios.\n "
                         f"Check successful! {lines} trios present in {bbaa}")
            logger.info(outstr)
            if remove_intermediate:
                shutil.rmtree(os.path.abspath(f'../output/{analysis}/DTparallel'))
                shutil.rmtree(os.path.abspath(f'../output/{analysis}/combined'))
            with open(output.check, 'w') as f:
                f.write(f"Analysis {analysis}: \n" + outstr + '\n')
