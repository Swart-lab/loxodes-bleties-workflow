rule all:
    input:
        'bleties/LmagMAC.LmagMAC.milraa_ies_fuzzy.no_overlap_repeats.gff3', # MAC reads on MAC assembly
        'bleties/LmagMIC.LmagMAC.milraa_ies_fuzzy.no_overlap_repeats.gff3',  # MIC reads on MAC assembly (typical use case for BleTIES)
        expand('variants/freebayes.{mode}.LmagMAC.g{cov_cutoff}.no_overlap_repeats.sort.vcf.gz', mode=['naive','diploid'], cov_cutoff=[200,400]),
        'variants/whatshap.freebayes.diploid.LmagMAC.g400.no_overlap_repeats.phased.vcf.gz'


rule haplotype_phase_whatshap:
    """Phase haplotype blocks with PacBio long reads

    Don't use Illumina for phasing because high coverage is not necessary
    (default max 15x), and whatshap will be confused by duplicate read names
    for paired end Illumina read libraries.
    """
    input:
        vcf='variants/freebayes.{mode}.{asm}.g{cov_cutoff}.no_overlap_repeats.sort.vcf.gz',
        asm=lambda wildcards: config['asm'][wildcards.asm],
        bam='/tmp/bleties/minimap2.comb.{asm}.rg_tag.bam'
    output:
        vcf='variants/whatshap.freebayes.{mode}.{asm}.g{cov_cutoff}.no_overlap_repeats.phased.vcf.gz',
        out='variants/whatshap.freebayes.{mode}.{asm}.g{cov_cutoff}.no_overlap_repeats.phased.out',
    log: 'logs/haplotype_phase_whatshap.{mode}.{asm}.g{cov_cutoff}.log'
    conda: 'envs/variants.yml'
    shell:
        r"""
        whatshap phase -o {output.vcf} --reference={input.asm} {input.vcf} {input.bam} > {output.out} 2> {log};
        bcftools index --threads {output.vcf};
        """


rule sort_index_vcf:
    input:
        'variants/freebayes.{mode}.{asm}.g{cov_cutoff}.no_overlap_repeats.vcf.gz'
    output:
        sort='variants/freebayes.{mode}.{asm}.g{cov_cutoff}.no_overlap_repeats.sort.vcf.gz',
        idx='variants/freebayes.{mode}.{asm}.g{cov_cutoff}.no_overlap_repeats.sort.vcf.gz.csi'
    conda: 'envs/variants.yml'
    shell:
        r"""
        bcftools sort -o {output.sort} -O z {input};
        bcftools index {output.sort}
        """


rule concatenate_vcf:
    input: 'variants/split_vcf.{mode}.{asm}.g{cov_cutoff}/'
    output: 'variants/freebayes.{mode}.{asm}.g{cov_cutoff}.no_overlap_repeats.vcf.gz'
    conda: 'envs/variants.yml'
    threads: 4
    shell:
        r"""
        bcftools concat --threads {threads} -o {output} -O z {input}/*.vcf
        """


rule run_freebayes:
    """Run Freebayes commands with Parafly

    Naive variant calling "hangs" on some regions. This step may have to be
    done partly manually.
    """
    input:
        'variants/freebayes.{mode}.{asm}.g{cov_cutoff}.no_overlap_repeats.cmds'
    output:
        temp(directory('variants/split_vcf.{mode}.{asm}.g{cov_cutoff}/'))
    conda: 'envs/variants.yml'
    threads: 16
    log: 'logs/run_freebayes.{mode}.{asm}.g{cov_cutoff}.log'
    shell:
        # Retry failed commands max 1 time (default is 10); 0 doesn't work for some reason
        r"""
        mkdir -p {output}
        workflow/ParaFly/bin/ParaFly -c {input} -CPU {threads} -max_retry 1 -shuffle &> {log}
        """


rule freebayes_commands_diploid:
    """Call variants assuming diploid genome

    Variant frequencies from naive variant calling show that strain is likely
    to be diploid. Call variants again with diploid mode, which should also
    proceed faster without hanging on some regions. Do not filter Q20 because
    that seems to be overly strict; bowtie2 assigns very conservative MAPQ
    scores for some reason. This could be because apparent PCR duplicates are
    scored 0 or 1, but I'm not sure exactly why.

    Prepare list of commands to parallelize with ParaFly, better job management
    than freebayes-parallel, which leaves some cores unused for long periods if
    individual jobs take longer time than others to run.
    """
    input:
        asm=lambda wildcards: config['asm'][wildcards.asm],
        regions='variants/ref.{asm}.no_overlap_repeats.bed',
        bam='/tmp/bleties/bowtie2.comb.{asm}.rg_tag.bam'
    output:
        'variants/freebayes.diploid.{asm}.g{cov_cutoff}.no_overlap_repeats.cmds'
    run:
        with open(output[0], 'w') as fh_out:
            with open(input.regions, 'r') as fh_in:
                for line in fh_in:
                    [chrom, start, stop] = line.rstrip().split("\t")
                    # Use `timeout` command to terminate jobs that take longer than 2 hours
                    cmd_out = f'timeout -v 2h bash -c \"freebayes -f {input.asm} -g {wildcards.cov_cutoff} {input.bam} --region {chrom}:{start}-{stop} > variants/split_vcf.diploid.{wildcards.asm}.g{wildcards.cov_cutoff}/{chrom}_{start}_{stop}.diploid.q20.vcf\"'
                    fh_out.write(cmd_out + "\n")


rule freebayes_commands_naive:
    """Call variants without assuming ploidy

    We do not know the ploidy a priori, nor can we rule out possibility that
    the population is a mixture of strains. Call variants without assuming
    ploidy then empirically plot variant frequencies to determine whether
    diploid or not.
    """
    input:
        asm=lambda wildcards: config['asm'][wildcards.asm],
        regions='variants/ref.{asm}.no_overlap_repeats.bed',
        bam='/tmp/bleties/bowtie2.comb.{asm}.rg_tag.bam'
    output:
        'variants/freebayes.naive.{asm}.g{cov_cutoff}.no_overlap_repeats.cmds'
    run:
        with open(output[0], 'w') as fh_out:
            with open(input.regions, 'r') as fh_in:
                for line in fh_in:
                    [chrom, start, stop] = line.rstrip().split("\t")
                    # Use `timeout` command to terminate jobs that take longer than 2 hours
                    # freebayes does not terminate after several hours in some genomic regions
                    # and with naive mode, so we set a timeout to prevent those jobs
                    # accumulating in the ParaFly queue.
                    cmd_out = f'timeout -v 2h bash -c \"freebayes -f {input.asm} -g {wildcards.cov_cutoff} --haplotype-length 0 --min-alternate-count 1 --min-alternate-fraction 0 --pooled-continuous {input.bam} --region {chrom}:{start}-{stop} | vcffilter -f \'QUAL > 20\' > variants/split_vcf.naive.{wildcards.asm}.g{wildcards.cov_cutoff}/{chrom}_{start}_{stop}.naive.q20.vcf\"'
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
    """List of contig lengths"""
    input:
        lambda wildcards: config['ctgbed'][wildcards.asm]
    output:
        'variants/ref.{asm}.contig_lengths'
    shell:
        r"""
        cut -f1,3 {input} > {output}
        """


rule tag_pacbio_mappings:
    """Combine MAC and MIC mappings into single BAM file with RG tags"""
    input:
        mic_reads=lambda wildcards: config['mapping_pb'][wildcards.asm]['LmagMIC'],
        mac_reads=lambda wildcards: config['mapping_pb'][wildcards.asm]['LmagMAC'],
    output:
        bam='/tmp/bleties/minimap2.comb.{asm}.rg_tag.bam',
        bai='/tmp/bleties/minimap2.comb.{asm}.rg_tag.bam.bai'
    conda: 'envs/variants.yml'
    shell:
        r"""
        bamaddrg -b {input.mac_reads} -s mac -b {input.mic_reads} -s mic > {output};
        samtools index {output.bam}
        """

rule tag_illumina_mappings:
    """Combine MAC and MIC mappings into single BAM file with RG tags"""
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
    """Filter BleTIES output to remove low-complexity repeat regions"""
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
    """Predict IESs with BleTIES MILRAA in CCS mode"""
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
