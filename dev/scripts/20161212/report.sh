echo "number of reads with mismatch rate <= 10% && length >= 10kb:"
samtools view /share/PI/mrivas/data/nanopore-wgs-consortium/nanopore-wgs.25000.sorted.10k.mapq50.ext.bam|wc -l

echo "number of filtered reads that contain at least one mismatch with Q-value >= 20:"
cat /share/PI/mrivas/data/nanopore-wgs-consortium/nanopore-wgs.25000.sorted.10k.mapq50.ext.snplist|awk '{if(NR>1){print $2}}'|wc -l

echo "total number of mismatch with Q-value >= 20 in filtered reads:"
cat /share/PI/mrivas/data/nanopore-wgs-consortium/nanopore-wgs.25000.sorted.10k.mapq50.ext.snplist|awk '{if(NR>1){print $2}}'|sed -e 's/;/\n/g'|wc -l

echo "total number of mismatches whose positions are present in dbSNPs:"
cat /share/PI/mrivas/data/nanopore-wgs-consortium/nanopore-wgs.25000.sorted.10k.mapq50.ext.w0.snps|grep -v '>' |wc -l
