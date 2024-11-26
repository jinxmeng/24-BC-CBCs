#dwn
cat sample_name | parallel -j 4 run_dwn_sra.sh {} raw_fq

#qc
cat sample_name | parallel -j 4 run_fastp_rmhost.sh raw_fq/{}_1.fastq.gz human clean_fq/{}

#mapping
seqkit replace -p "\s.*" -r "" genomes/GCF_900461125.1_49699_D02_genomic.fna > Blautia_coccoides.fna
minimap2 -d Blautia_coccoides.mmi Blautia_coccoides.fna
cat /share/data1/mjx/meta/20240803_IBD_pub_data_metagenome/SchirmerM_2018.PRJNA389280/clean_fq.filepath | parallel -j 10 --colsep="\t" run_minimap2.sh {2} Blautia_coccoides.mmi SchirmerM_2018.PRJNA389280/{1}
ls SchirmerM_2018.PRJNA389280/*bam | parallel -j 10 --plus samtools coverage -d 0 -o {..}.cvm {}
ls *cvm | parallel -j 20 -q perl -e 'open I, "$ARGV[0]"; $x=0;while(<I>){chomp; next if /numreads/; @s=split/\s+/; $x+=$s[3]}; print "$ARGV[0]\t$x\n"' {} | sed 's/.cvm//g' | csvtk join --left-join -t -f 1 - /share/data1/mjx/meta/20240803_IBD_pub_data_metagenome/SchirmerM_2018.PRJNA389280/clean_fq.lib_size | csvtk add-header -t -n "sample,rc,libs" > ../SchirmerM_2018.PRJNA389280.profile