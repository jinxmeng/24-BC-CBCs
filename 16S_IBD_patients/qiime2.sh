#seqs
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-format PairedEndFastqManifestPhred33V2 --input-path sample_list --output-path 03.seqs/seqs.qza
qiime demux summarize --i-data 03.seqs/seqs.qza --o-visualization 03.seqs/seqs.qzv &

#dada2
qiime dada2 denoise-paired --i-demultiplexed-seqs 03.seqs/seqs.qza --p-trunc-len-f 0 --p-trunc-len-r 0 --o-representative-sequences 05.dada2/rep_seqs.qza --o-table 05.dada2/table.qza --o-denoising-stats 05.dada2/dada2_stats.qza --p-n-threads 100 &
qiime feature-table filter-features --i-table 05.dada2/table.qza --p-min-frequency 1 --p-min-samples 3 --o-filtered-table 05.dada2/table_f.qza
qiime tools export --input-path 05.dada2/table_f.qza --output-path 05.dada2/
biom convert -i 05.dada2/feature-table.biom -o 05.dada2/table_f.tsv --to-tsv
qiime feature-table summarize --i-table 05.dada2/table_f.qza --o-visualization 05.dada2/table_f.qzv --m-sample-metadata-file sample_group
qiime feature-table filter-seqs --i-data 05.dada2/rep_seqs.qza --i-table 05.dada2/table_f.qza --o-filtered-data 05.dada2/rep_seqs_f.qza

#div
qiime diversity core-metrics --i-table 05.dada2/table_f.qza --p-sampling-depth 47875 --m-metadata-file sample_group --p-n-jobs 100 --o-rarefied-table 06.div/table_rarefied.qza --output-dir 06.div/core_diversity
qiime tools export --input-path 06.div/table_rarefied.qza --output-path 06.div/
biom convert -i 06.div/feature-table.biom -o 06.div/table_rarefied.tsv --to-tsv

#tree
qiime phylogeny align-to-tree-mafft-fasttree --i-sequences 05.dada2/rep_seqs_f.qza --o-alignment 08.tree/aligned.qza --o-masked-alignment 08.tree/masked-aligned.qza --o-tree 08.tree/unrooted-tree.qza --o-rooted-tree 08.tree/rooted-tree.qza --p-n-threads 60

#tax
export TMPDIR=/share/data1/mjx/tmp/qiime2/
qiime feature-classifier classify-sklearn --i-classifier /share/data1/database/Silva/341F_806R/silva-138-99-338F_806R_classifier.qza --i-reads 05.dada2/rep_seqs_f.qza --o-classification 07.tax/tax.qza --p-n-jobs -1 &
qiime taxa barplot --i-table 06.div/table_rarefied.qza --i-taxonomy 07.tax/tax.qza --m-metadata-file sample_group --o-visualization 07.tax/tax.qzv