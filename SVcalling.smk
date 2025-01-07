import os
import gzip

configfile: "SVcalling_config.yaml"

# 提取文件名的基部分（去除路径和扩展名）
ref_basename=os.path.splitext(os.path.basename(config["ref"]))[0]
fastq_suffix=config.get("fastq_suffix")

rule all:
    input:
        expand("vcf/sv/germline/{sample}.bcf", sample=config["sample"]),
        "vcf/sv/final.vcf.gz",
        "mapping/merged_depth_stats.txt"

rule bwa_index:
    input:
        reference_genome=config["ref"]
    output:
        "genome_index/{ref_basename}.amb",
        "genome_index/{ref_basename}.ann",
        "genome_index/{ref_basename}.bwt",
        "genome_index/{ref_basename}.pac",
        "genome_index/{ref_basename}.sa"
    log:
        "logs/index/bwa_index_{ref_basename}.log"
    shell:
        """
        bwa index \
        -p genome_index/{wildcards.ref_basename} \
        {input.reference_genome} \
        2> {log}
        """

rule samtools_fai_index:
    input:
        reference_genome=config["ref"]
    output:
        "genome_index/{ref_basename}.fai"
    log:
        "logs/index/samtools_index_{ref_basename}.log"
    shell:
        """
        samtools faidx {input.reference_genome} 2> {log}
        cp {input.reference_genome}.fai {output}
        """

rule gatk_dict_index:
    input:
        reference_genome=config["ref"]
    output:
        "genome_index/{ref_basename}.dict"
    log:
        "logs/index/gatk_index_{ref_basename}.log"
    shell:
        """
        gatk CreateSequenceDictionary \
        -R {input.reference_genome} \
        -O {output} \
        2> {log}
        """

rule QualityControlfastp:
    input:
        f"raw_data/{{sample}}_clean_1{fastq_suffix}",
        f"raw_data/{{sample}}_clean_2{fastq_suffix}"
    output:
        "clean_data/{sample}_1_clean.fq.gz",
        "clean_data/{sample}_2_clean.fq.gz",
        "logs/fastp/fastp_report/{sample}.fastp.html"
    threads: 2
    params:
        qualified_quality_phred=config["qualified_quality_phred"],
        unqualified_percent_limit=config["unqualified_percent_limit"],
        trim_front=config["trim_front"]
    log:
        "logs/fastp/{sample}.log"
    shell:
        """
        fastp \
        --thread {threads} \
        -i {input[0]} \
        -I {input[1]} \
        -o {output[0]} \
        -O {output[1]} \
        -h {output[2]} \
        -j /dev/null \
        -q {params.qualified_quality_phred} \
        -u {params.unqualified_percent_limit} \
        -f {params.trim_front} \
        2> {log}
        """

rule BWA_map:
    input:
        r1="clean_data/{sample}_1_clean.fq.gz",
        r2="clean_data/{sample}_2_clean.fq.gz",
        bwa_index=expand("genome_index/{ref_basename}.{ext}", ref_basename=ref_basename, ext=["amb", "ann", "bwt", "pac", "sa"])
    output:
        temp("mapping/{sample}.sorted.bam")
    threads: 8
    params:
        rg=r"@RG\tID:{sample}\tSM:{sample}\tLB:{sample}\tPL:ILLUMINA",
        prefix=f"genome_index/{ref_basename}"
    log:
        "logs/bwa/bwa_map_{sample}.log"
    shell:
        """
        bwa mem \
        -R '{params.rg}' \
        -t {threads} \
        {params.prefix} {input.r1} {input.r2} \
        | samtools sort -@ {threads} -o {output} \
        2> {log}
        """

rule RemoveDuplicates:
    input:
        "mapping/{sample}.sorted.bam"
    output:
        "mapping/{sample}.sorted.markdup.bam",
        "mapping/markdup_metrics/{sample}.sorted.markdup_metrics.txt"
    log:
        "logs/bwa/RemoveDuplicates_{sample}.log"
    shell:
        """
        gatk MarkDuplicates \
        -I {input} \
        -O {output[0]} \
        -M {output[1]} \
        --REMOVE_DUPLICATES true \
        --CREATE_INDEX true \
        2> {log}
        """

rule BAMDepthStat:
    input:
        bam="mapping/{sample}.sorted.markdup.bam"
    output:
        temp("mapping/{sample}.chr.stat.gz")
    log:
        "logs/depth/{sample}_pandepth.log"
    threads: 2
    shell:
        """
        pandepth \
        -i {input.bam} \
        -o mapping/{wildcards.sample} \
        -t {threads} \
        2> {log}
        """

rule MergeDepthStats:
    input:
        expand("mapping/{sample}.chr.stat.gz", sample=config["sample"])
    output:
        "mapping/merged_depth_stats.txt"
    log:
        "logs/depth/merge_depth_stats.log"
    run:
        with open(output[0], 'w') as out_file:
            for stat_file in input:
                sample=stat_file.split("/")[-1].split(".")[0]
                with gzip.open(stat_file, 'rt') as f: 
                    last_line=f.readlines()[-1].strip()
                    out_file.write(f"{sample}\t{last_line}\n")

rule DellyIndividualSVcalling:
    input:
        bam="mapping/{sample}.sorted.markdup.bam",
        ref=config["ref"]
    output:
        bcf="vcf/sv/germline/{sample}.bcf"
    conda:
        config["conda_env"]["delly_conda"]
    log:
        "logs/sv/{sample}.delly_individual_call.log"
    threads: 2
    shell:
        """
        delly call \
        -g {input.ref} \
        -o {output.bcf} \
        {input.bam} \
        2> {log}
        """
        
rule DellyIndividualSVmerge:
    input:
        expand("vcf/sv/germline/{sample}.bcf", sample=config["sample"])
    output:
        sitebcf="vcf/sv/sites.bcf"
    conda:
        config["conda_env"]["delly_conda"]
    log:
        "logs/sv/merged_site_sv.log"
    shell:
        """
        delly merge \
        -o {output.sitebcf} \
        {input} \
        2> {log}
        """
        
rule DellySVgenotype:
    input:
        sitebcf="vcf/sv/sites.bcf",
        ref=config["ref"],
        bam="mapping/{sample}.sorted.markdup.bam"
    output:
        genobcf="vcf/sv/germline/{sample}.geno.bcf"
    conda:
        config["conda_env"]["delly_conda"]
    log:
        "logs/sv/{sample}.geno.log"
    threads: 2
    shell:
        """
        delly call \
        -g {input.ref} \
        -v {input.sitebcf} \
        -o {output.genobcf} \
        {input.bam} \
        2> {log}
        """
        
rule DellySVgenomerge:
    input:
        expand("vcf/sv/germline/{sample}.geno.bcf", sample=config["sample"])
    output:
        genobcf="vcf/sv/merged.geno.bcf"
    log:
        "logs/sv/merged_geno_sv.log"
    shell:
        """
        bcftools merge \
        -m id \
        -O b \
        -o {output.genobcf} \
        {input} \
        2> {log}
        """

rule BCFindex:
    input:
        "vcf/sv/merged.geno.bcf"
    output:
        "vcf/sv/merged.geno.bcf.csi"
    log:
        "logs/sv/merged.geno.bcf.index.log"
    shell:
        """
        bcftools index \
        -c \
        -o {output} \
        {input} \
        2> {log}
        """

rule DellySVfilter:
    input:
        "vcf/sv/merged.geno.bcf",
        "vcf/sv/merged.geno.bcf.csi"
    output:
        "vcf/sv/final.bcf"
    conda:
        config["conda_env"]["delly_conda"]
    log:
        "logs/sv/filter_sv.log"
    shell:
        """
        delly filter \
        -f germline \
        -o {output} \
        {input[0]} \
        2> {log}
        """
        
rule bcf2vcf:
    input:
        "vcf/sv/final.bcf"
    output:
        "vcf/sv/final.vcf.gz"
    log:
        "logs/sv/bcf2vcf.log"
    shell:
        """
        bcftools convert \
        -Oz \
        -o {output} \
        {input} \
        2> {log}
        """
