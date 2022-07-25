import glob, os

# This snakefile runs lpWGS samples through the CN caller WisecondorX, and improves the output through shifting and purity estimation
# The workflow was made by Kristena Daley on Dec 9, 2020 (Started on Dec 9, 2020, finished on Aug 31, 2021)

# Declare a variable for where to find your bams
bam_path_all = "/projects/rmorin/projects/NHL_ctDNA_analysis/wisecondorX/lpWGS_MarkDuplicates_test/"
bam_path_normals = "/projects/rmorin/projects/NHL_ctDNA_analysis/wisecondorX/lpWGS_all_500kb_bincorrection/00_reference_bam/"
base_path = "/projects/rmorin/projects/NHL_ctDNA_analysis/wisecondorX/lpWGS_all_500kb_bincorrection/"
base_path_unadjusted = "/projects/rmorin/projects/NHL_ctDNA_analysis/wisecondorX/lpWGS_all_500kb_bincorrection/03_modified_code_output_all/01_unadjusted_output/"
base_path_adjusted = "/projects/rmorin/projects/NHL_ctDNA_analysis/wisecondorX/lpWGS_all_500kb_bincorrection/03_modified_code_output_all/02_adjusted_output/""

samples = []
for bamfile in glob.glob(bam_path_all + "*.bam"):
    base = os.path.basename(bamfile)
    sampleid = base.split(".")[0]
    samples.append(sampleid)

rule all:
    input:
        expand("/projects/rmorin/projects/NHL_ctDNA_analysis/wisecondorX/lpWGS_all_500kb_bincorrection/03_modified_code_output_all/02_adjusted_output/{sample}_output/new_offset_mean_and_purity.txt",sample=samples)

##
## Section 1: 
##    Step 1) Run WisecondorX original parameters (500kb bins)
##              a) bam_to_npz_samples - Converts bam files to npz file format for all samples and reference files
##              b) npz_ref_creation - Creates a reference npz file from the individual reference npz files. Note these should all be in a separate directory from the sample.npz files. The docs for WisecondorX mention that having at least 50 samples for reference creation is ideal for optimal results (normalization). 
##              c) CN_prediction - Predicts copy number alterations and outputs multiple files and plots. These may require shifting due to the way wisecondorX calculates aberration cutoffs and centers log ratios of bins. 
##    Step 2) Plot_CN_density_unshifted - Plots unshifted CN plot with their corresponding density plot from the mixture model (This is calculated using one cluster/classification in mclust).
##    Step 3) Determine_offset_mean_and_purity - Calculates the offset mean and purity estimation for each sample, through a gaussian mixture model.
## 

rule bam_to_npz_samples:
    input:
        f"{bam_path_all}" + "{sample}.LOWPASS.markDupl.bam"
    output:
        f"{base_path}" + "01_samples_npz/{sample}.npz"
    params:
        wisecondorX="/home/krdaley/anaconda3/envs/wisecondorX/bin/WisecondorX"
    conda:
        "/projects/rmorin/projects/NHL_ctDNA_analysis/snakemake/wisecondorX.yaml"
    shell:
        "{params.wisecondorX} convert {input} {output} --binsize 500000"

rule bam_to_npz_normals:
    input:
        f"{bam_path_normals}" + "{sample}.bam"
    output:
        f"{base_path}" + "01_reference_npz/{sample}.npz"
    params:
        wisecondorX="/home/krdaley/anaconda3/envs/wisecondorX/bin/WisecondorX"
    conda:
        "/projects/rmorin/projects/NHL_ctDNA_analysis/snakemake/wisecondorX.yaml"
    shell:
        "{params.wisecondorX} convert {input} {output} --binsize 500000"

#Chris Rushton made a reference and it is located here: 
#/projects/rmorin_scratch/NHL_ctDNA_analysis_scratch/wisecondorX/01wisecondorx_ref_fromchris.500kb.2020-11-24.86.npz

rule npz_ref_creation:
    input:
        f"{base_path}" + "01_reference_npz/"
    output:
        f"{base_path}" + "01_reference_npz/WisecondorX_reference_500kb.npz"
    params:
        wisecondorX="/home/krdaley/anaconda3/envs/wisecondorX/bin/WisecondorX"
    conda:
        "/projects/rmorin/projects/NHL_ctDNA_analysis/snakemake/wisecondorX.yaml"
    shell:
        "{params.wisecondorX} newref {input}/*.npz {output} --binsize 500000"


rule CN_prediction:
    input:
        npz=f"{base_path}" + "01_samples_npz/{sample}.npz",
        ref="/projects/rmorin_scratch/NHL_ctDNA_analysis_scratch/wisecondorX/01wisecondorx_ref_fromchris.500kb.2020-11-24.86.npz"
    output:
        bed=f"{base_path_unadjusted}" + "{sample}_output/{sample}_bins.bed",
        png=f"{base_path_unadjusted}" + "{sample}_output/{sample}.plots/genome_wide.png"
    params:
        wisecondorX="/home/krdaley/anaconda3/envs/wisecondorX/bin/WisecondorX"
    conda:
        "/projects/rmorin/projects/NHL_ctDNA_analysis/snakemake/wisecondorX.yaml"
    log: 
        f"{base_path}" + "01_samples_npz/log/{sample}.wisecondor.log"
    shell:
        "{params.wisecondorX} predict {input.npz} {input.ref} $(dirname {output.bed})/{wildcards.sample} --bed --plot --ylim [-2,2] > {log} 2>&1 "


rule Plot_CN_and_density_unshifted:
    input: 
        bed=f"{base_path_unadjusted}" + "{sample}_output/{sample}_bins.bed",
        png=f"{base_path_unadjusted}" + "{sample}_output/{sample}.plots/genome_wide.png" 
    output:
        pdf=f"{base_path_unadjusted}" + "{sample}_output/genome_wide_density.pdf"
    conda:
        "/projects/rmorin/projects/NHL_ctDNA_analysis/snakemake/R.yaml"
    log:
        f"{base_path}" + "01_samples_npz/log_merged_plots/{sample}.wisecondor.log"
    shell:
        "Rscript /home/krdaley/anaconda3/envs/wisecondorX/lib/python3.8/site-packages/wisecondorX/include/cn_dist_plot_merge.R {input.bed} {input.png} {output.pdf} > {log} 2>&1 "

rule Determine_offset_mean_and_purity:
    input:
        bed=f"{base_path_unadjusted}" + "{sample}_output/{sample}_bins.bed",
        pdf=f"{base_path_unadjusted}" + "{sample}_output/genome_wide_density.pdf"
    output:
        txt=f"{base_path_unadjusted}" + "{sample}_output/new_offset_mean_and_purity.txt"
    conda:
        "/projects/rmorin/projects/NHL_ctDNA_analysis/snakemake/R.yaml"
    log:
        f"{base_path}" + "01_samples_npz/log_purity_est/{sample}.wisecondor.log"
    shell:
        "Rscript /home/krdaley/anaconda3/envs/wisecondorX/lib/python3.8/site-packages/wisecondorX/include/offset_mean_and_purity_calc.R {input.bed} {output.txt} > {log} 2>&1 "

##
## Section 2: 
##    Step 1) Shift_mean_add_beta - Reruns WisecondorX using the offset mean and beta (beta = purity if it is between 0 and 1, otherwise beta is not used and no amplifications will be seen in the output plots).
##    Step 2) Plot_CN_and_density_postshift - Plots shifted CN plots with corresponding density plots. Side note; the "offset" line in the density plot is not relevant after shifting is complete. 
##    Step 3) Recalculate_shifted_mean_and_purity - New purity estimation (equivalent to ichorCNAs tumour fraction output value) for each sample. Offset mean is calculated again, but is not relevant as shifting has already occurred, so the focus here is on the new purity estimation. 
##

#This rule can be finicky due to the bash shell portion, so note that the do_rerun step works fine on a mac, but has not been used on a windows computer, so there might be an error about return carriages or something? I havent had it but I was told it might occur in the future, if so, the do_rerun step might need to be tweaked..just a thought. 
rule Shift_mean_add_beta:
    input:
        npz=f"{base_path}" + "01_samples_npz/{sample}.npz",
        ref="/projects/rmorin_scratch/NHL_ctDNA_analysis_scratch/wisecondorX/01wisecondorx_ref_fromchris.500kb.2020-11-24.86.npz",
        mean_purity_file=f"{base_path_unadjusted}" + "{sample}_output/new_offset_mean_and_purity.txt"
    output:
        bed=f"{base_path_adjusted}" + "{sample}_output/{sample}_bins.bed",
        png=f"{base_path_adjusted}" + "{sample}_output/{sample}.plots/genome_wide.png"
        #done=f"{base_path_adjusted}" + "{sample}_output/{sample}.done"
    params:
        wisecondorX="/home/krdaley/anaconda3/envs/wisecondorX/bin/WisecondorX"
    conda:
        "/projects/rmorin/projects/NHL_ctDNA_analysis/snakemake/wisecondorX.yaml"
    log: 
        f"{base_path}" + "01_samples_npz/log/{sample}.wisecondor.log"
    shell:
        """
        offset_value=$(head -1 {input.mean_purity_file} | grep -Eo [0-9]+\.[0-9]+.+)
        purity_value=$(tail -2 {input.mean_purity_file} | grep -Eo [0-9]+\.[0-9]+)
        do_rerun=$(tail -1 {input.mean_purity_file} | cut -d'=' -f2 | cut -d' ' -f2)
        if [[ $do_rerun == "TRUE" ]]
        then
            {params.wisecondorX} predict {input.npz} {input.ref} $(dirname {output.bed})/{wildcards.sample} --bed --plot --offset_neut_peak $offset_value --ylim [-2,2] > {log} 2>&1 
        else
            {params.wisecondorX} predict {input.npz} {input.ref} $(dirname {output.bed})/{wildcards.sample} --bed --plot --offset_neut_peak $offset_value --ylim [-2,2] --beta $purity_value > {log} 2>&1
        fi
        """


rule Plot_CN_and_density_postshift:
    input: 
        bed=f"{base_path_adjusted}" + "{sample}_output/{sample}_bins.bed",
        png=f"{base_path_adjusted}" + "{sample}_output/{sample}.plots/genome_wide.png"
    output:
        pdf=f"{base_path_adjusted}" + "{sample}_output/genome_wide_density.pdf"
    conda:
        "/projects/rmorin/projects/NHL_ctDNA_analysis/snakemake/R.yaml"
    log:
        f"{base_path}" + "01_samples_npz/log_merged_plots/{sample}.wisecondor.log"
    shell:
        "Rscript /home/krdaley/anaconda3/envs/wisecondorX/lib/python3.8/site-packages/wisecondorX/include/cn_dist_plot_merge.R {input.bed} {input.png} {output.pdf} > {log} 2>&1 "


rule Recalculate_shifted_mean_and_purity:
    input:
        bed=rules.Shift_mean_add_beta.output.bed,
        pdf=rules.Plot_CN_and_density_postshift.output.pdf
    output:
        txt=f"{base_path_adjusted}" + "{sample}_output/new_offset_mean_and_purity.txt"
    conda:
        "/projects/rmorin/projects/NHL_ctDNA_analysis/snakemake/R.yaml"
    log:
        f"{base_path}" + "01_samples_npz/log_purity_est/{sample}.wisecondor.shifted.log"
    shell:
        "Rscript /home/krdaley/anaconda3/envs/wisecondorX/lib/python3.8/site-packages/wisecondorX/include/offset_mean_and_purity_calc.R {input.bed} {output.txt} > {log} 2>&1 "
