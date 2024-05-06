#!/bin/bash

bcftools --version

# individual filter steps should be piped rather than combined into a single
# command, and the output to be piped should be in uncompressed bcf format; see
# man page for bcftools

bcftools view -O u -s mac variants/freebayes.diploid.LmagMAC.g400.no_overlap_repeats.sort.vcf.gz | bcftools view -O u -g het | bcftools view -O z > variants/freebayes.diploid.LmagMAC.g400.no_overlap_repeats.mac.het.piped.vcf.gz
bcftools stats variants/freebayes.diploid.LmagMAC.g400.no_overlap_repeats.mac.het.piped.vcf.gz > variants/freebayes.diploid.LmagMAC.g400.no_overlap_repeats.mac.het.piped.stats
bcftools view -O u -s mic variants/freebayes.diploid.LmagMAC.g400.no_overlap_repeats.sort.vcf.gz | bcftools view -O u -g het | bcftools view -O z > variants/freebayes.diploid.LmagMAC.g400.no_overlap_repeats.mic.het.piped.vcf.gz
bcftools stats variants/freebayes.diploid.LmagMAC.g400.no_overlap_repeats.mic.het.piped.vcf.gz > variants/freebayes.diploid.LmagMAC.g400.no_overlap_repeats.mic.het.piped.stats

