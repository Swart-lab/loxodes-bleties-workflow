rule all:
    input:
        'bleties/LmagMAC.LmagMAC.milraa_ies_fuzzy.no_overlap_repeats.gff3', # MAC reads on MAC assembly
        'bleties/LmagMIC.LmagMAC.milraa_ies_fuzzy.no_overlap_repeats.gff3',  # MIC reads on MAC assembly (typical use case for BleTIES)
        '/tmp/bleties/bowtie2.comb.LmagMAC.rg_tag.bam',
        'variants/freebayes.LmagMAC.no_overlap_repeats.cmds'


rule freebayes_commands:
    input:
        asm=lambda wildcards: config['asm'][wildcards.asm],
        regions='variants/ref.{asm}.no_overlap_repeats.bed',
        bam='/tmp/bleties/bowtie2.comb.{asm}.rg_tag.bam'
    output:
        'variants/freebayes.{asm}.no_overlap_repeats.cmds'
    run:
        with open(output[0], 'w') as fh_out:
            with open(input.regions, 'r') as fh_in:
                for line in fh_in:
                    [chrom, start, stop] = line.rstrip().split("\t")
                    cmd_out = f'freebayes -f {input.asm} -g 400 --haplotype-length 0 --min-alternate-count 1 --min-alternate-fraction 0 --pooled-continuous {input.bam} --region {chrom}:{start}-{stop} | vcffilter -f \"QUAL > 20\" > variants/split_vcf/{chrom}_{start}_{stop}.naive.q20.vcf'
                    fh_out.write(cmd_out + "\n")



rule no_overlap_repeats:
    input:
        trf_gff=lambda wildcards: config['trf'][wildcards.asm],
        genome='variants/ref.{asm}.contig_lengths'
    output:
        'variants/ref.{asm}.no_overlap_repeats.bed'
    conda: 'envs/bedtools.yml'
    shell:
        r"""
        bedtools complement -i {input.trf_gff} -g {input.genome} > {output}
        """


rule contig_lengths:
    input:
        lambda wildcards: config['ctgbed'][wildcards.asm]
    output:
        'variants/ref.{asm}.contig_lengths'
    shell:
        r"""
        cut -f1,3 {input} > {output}
        """


rule tag_illumina_mappings:
    input:
        mic_reads=lambda wildcards: config['mapping_illumina'][wildcards.asm]['LmagMIC'],
        mac_reads=lambda wildcards: config['mapping_illumina'][wildcards.asm]['LmagMAC'],
    output:
        bam='/tmp/bleties/bowtie2.comb.{asm}.rg_tag.bam',
        bai='/tmp/bleties/bowtie2.comb.{asm}.rg_tag.bam.bai'
    conda: 'envs/variants.yml'
    shell:
        r"""
        bamaddrg -b {input.mac_reads} -s mac -b {input.mic_reads} -s mic > {output};
        samtools index {output.bam}
        """


rule bleties_no_overlap_repeats:
    input:
        trf_gff=lambda wildcards: config['trf'][wildcards.asm],
        fuzzy_gff='bleties/{ccs}.{asm}.milraa_ies_fuzzy.gff3'
    output:
        'bleties/{ccs}.{asm}.milraa_ies_fuzzy.no_overlap_repeats.gff3'
    threads: 2
    conda: 'envs/bedtools.yml'
    shell:
        r"""
        bedtools intersect -a {input.fuzzy_gff} -b {input.trf_gff} -v > {output}
        """


rule run_bleties:
    input:
        'bleties.{ccs}.{asm}.cmds'
    threads: 16
    output:
        gff='bleties/{ccs}.{asm}.milraa_ies.gff3',
        fasta='bleties/{ccs}.{asm}.milraa_ies.fasta',
        fuzzy_gff='bleties/{ccs}.{asm}.milraa_ies_fuzzy.gff3',
        fuzzy_fasta='bleties/{ccs}.{asm}.milraa_ies_fuzzy.fasta'
    log: 'logs/run_bleties.{ccs}.{asm}.log'
    conda: 'envs/bleties.yml'
    shell:
        r"""
        mkdir -p /tmp/bleties/{wildcards.ccs}.{wildcards.asm}
        workflow/ParaFly/bin/ParaFly -c {input} -CPU {threads} -shuffle &> {log}
        cat /tmp/bleties/{wildcards.ccs}.{wildcards.asm}/{wildcards.ccs}.{wildcards.asm}.*.milraa_ies.gff3 | grep -v '^#' > {output.gff}
        cat /tmp/bleties/{wildcards.ccs}.{wildcards.asm}/{wildcards.ccs}.{wildcards.asm}.*.milraa_ies.fasta > {output.fasta}
        cat /tmp/bleties/{wildcards.ccs}.{wildcards.asm}/{wildcards.ccs}.{wildcards.asm}.*.milraa_ies_fuzzy.gff3 | grep -v '^#' > {output.fuzzy_gff}
        cat /tmp/bleties/{wildcards.ccs}.{wildcards.asm}/{wildcards.ccs}.{wildcards.asm}.*.milraa_ies_fuzzy.fasta > {output.fuzzy_fasta}
        """


rule bleties_mappings_commands:
    input:
        asm=lambda wildcards: config['asm'][wildcards.asm],
        bam=lambda wildcards: config['mapping_pb'][wildcards.asm][wildcards.ccs]
    threads: 1
    output:
        'bleties.{ccs}.{asm}.cmds'
    shell:
        r"""
        for CTG in $(grep '>' {input.asm} | cut -f1 -d' ' | sed 's/>//'); do
            echo "bleties milraa --min_break_coverage 3 --min_del_coverage 5 --fuzzy_ies --type ccs --bam {input.bam} --ref {input.asm} --contig $CTG -o /tmp/bleties/{wildcards.ccs}.{wildcards.asm}/{wildcards.ccs}.{wildcards.asm}.$CTG" >> {output}
        done
        """
