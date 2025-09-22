# README
*Distortopia: simulation of long-read sequences, high-resolution mapping of meiotic crossovers, and identification of segregation distorters in **Arabidopsis**.* 

# download Athaliana and Alyrata genomes from JGI Phytozome
$ curl --cookie jgi_session=/api/sessions/e68504ec0dc38906b76d27c0e311eeea --output download.20250726.142947.zip -d "{\"ids\":{\"Phytozome-384\":{\"file_ids\":[\"585486967ded5e78cff8c52c\"],\"top_hit\":\"683db57de69ed9082626f62f\"},\"Phytozome-447\":{\"file_ids\":[\"587b0adf7ded5e4229d885ab\"],\"top_hit\":\"683db58de69ed9082626f687\"}},\"api_version\":\"2\"}" -H "Content-Type: application/json" https://files-download.jgi.doe.gov/filedownload/

# generate alt ref haploid genomes
$ python sim_haplotypes.py data/Ath 0.001 123 hap11 > data/ref1_hap1.fa
$ python sim_haplotypes.py data/Ath 0.001 321 hap12 > data/ref1_hap2.fa
$ python sim_haplotypes.py data/Aly 0.001 666 hap21 > data/ref2_hap1.fa
$ python sim_haplotypes.py data/Aly 0.001 999 hap22 > data/ref2_hap2.fa

# create alignments of each genome [minimap2 [options] target query]
$ minimap2 -x asm5 -t 6 data/ref1_hap2.fa data/ref1_hap1.fa > data/ref11xref12.paf

# generate long reads from gametes of an artificial hybrid with uniform recomb
$ sim_gamete_pool_long_reads ref1_hap1.fa ref1_hap2.fa 100_000 1e6 > ref1xref1_gametes_longreads.fq.gz
$ sim_gamete_pool_long_reads ref1_hap1.fa ref1_hap2.fa 100_000 1e6 > ref1xref2_gametes_longreads.fq.gz

# map reads to the reference
$ align_gamete_longreads REF1 ref1xref1_gamete_longreads.fq.gz > ref1xref1_gametes.bam
$ align_gamete_longreads REF1 ref1xref2_gamete_longreads.fq.gz > ref1xref2_gametes.bam

# call variants
$ gamete_variant_calling ref1xref1_gametes.bam > ref1xref1_gametes.vcf
$ gamete_variant_calling ref1xref1_gametes.bam > ref1xref2_gametes.vcf

# run crossover detection
$ colmapper REF1 ref1xref1_gametes.vcf > ref1xref1_gametes_crossovers.tsv
$ colmapper REF1 ref1xref2_gametes.vcf > ref1xref2_gametes_crossovers.tsv
