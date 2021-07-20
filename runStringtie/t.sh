cd '/home/yingh/data/GC-FOXL2-RNA/1_rna_abund'
stringtie -G ~/genome/Gallus_5.0/Gallus_gallus.Gallus_gallus-5.0.91.gff3 -o ./t.gtf -A t.txt ../bam/rmdup/BGC_0_1.rmdup.bam
