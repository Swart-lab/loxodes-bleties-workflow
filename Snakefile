rule all:
    input:
        # 'bleties/LmagMAC.LmagMAC.milraa_ies.gff3', # MAC reads on MAC assembly
        # 'bleties/LmagMIC.LmagMAC.milraa_ies.gff3'  # MIC reads on MAC assembly (typical use case for BleTIES)
        'bleties/LmagMAC.LmagMAC.milraa_ies_fuzzy.no_overlap_repeats.gff3', # MAC reads on MAC assembly
        'bleties/LmagMIC.LmagMAC.milraa_ies_fuzzy.no_overlap_repeats.gff3'  # MIC reads on MAC assembly (typical use case for BleTIES)


rule no_overlap_repeats:
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
        bam=lambda wildcards: config['mapping'][wildcards.asm][wildcards.ccs]
    threads: 1
    output:
        'bleties.{ccs}.{asm}.cmds'
    shell:
        r"""
        for CTG in $(grep '>' {input.asm} | cut -f1 -d' ' | sed 's/>//'); do
            echo "bleties milraa --min_break_coverage 3 --min_del_coverage 5 --fuzzy_ies --type ccs --bam {input.bam} --ref {input.asm} --contig $CTG -o /tmp/bleties/{wildcards.ccs}.{wildcards.asm}/{wildcards.ccs}.{wildcards.asm}.$CTG" >> {output}
        done
        """
