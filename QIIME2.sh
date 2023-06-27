#!/bin/bash

conda activate qiime2-2021.2

qiime tools import --type SampleData[SequencesWithQuality] --input-path forward_only --input-format CasavaOneEightSingleLanePerSampleDirFmt   --output-path forward.all.qza

qiime demux summarize --i-data forward.all.qza --o-visualization forward.all.qzv

qiime dada2 denoise-single --i-demultiplexed-seqs forward.all.qza --p-trim-left 0 --p-trunc-len 295 --p-n-threads 0 --o-table table.forward.all.qza --o-representative-sequences rep-seqs.forward.all.qza --o-denoising-stats denoising-stats.forward.all.qza --verbose

qiime feature-table summarize --i-table table.forward.all.qza --o-visualization table.forward.all.qzv --m-sample-metadata-file sample-metadata.tsv

qiime phylogeny align-to-tree-mafft-fasttree --i-sequences rep-seqs.forward.all.qza --o-alignment aligned-rep-seqs.qza --o-masked-alignment masked-aligned-rep-seqs.qza --o-tree unrooted-tree.qza --o-rooted-tree rooted-tree.qza

qiime diversity core-metrics-phylogenetic --i-phylogeny rooted-tree.qza --i-table table.forward.all.qza --p-sampling-depth 13085 --m-metadata-file sample-metadata.tsv --output-dir core-metrics-results

qiime tools import --type 'FeatureData[Taxonomy]' --input-format HeaderlessTSVTaxonomyFormat --input-path 99_otu_taxonomy.txt --output-path ref-taxonomy.qza

qiime feature-classifier extract-reads --i-sequences 99_otus.qza  --p-f-primer AGAGTTTGATCMTGGCTCAG --p-r-primer GWATTACCGCGGCKGCTG --o-reads ref-seqs.qza

qiime feature-classifier fit-classifier-naive-bayes --i-reference-reads ref-seqs.qza --i-reference-taxonomy ref-taxonomy.qza --o-classifier classifier.qza –verbose

qiime feature-classifier classify-sklearn --i-classifier classifier.qza --i-reads rep-seqs.qza --o-classification taxonomy.qza –verbose

qiime feature-classifier extract-reads  --i-sequences silva-138-99-seqs.qza --p-f-primer AGAGTTTGATCMTGGCTCAG --p-r-primer GWATTACCGCGGCKGCTG --p-n-jobs 12 --o-reads ref-seqs.qza --verbose

qiime feature-classifier fit-classifier-naive-bayes --i-reference-reads ref-seqs.qza --i-reference-taxonomy silva-138-99-tax.qza --o-classifier classifier.qza --verbose

qiime feature-classifier classify-sklearn --i-classifier classifier.qza --i-reads rep-seqs.qza --o-classification taxonomy.qza –verbose

qiime taxa filter-table --i-table table.forward.all.qza --i-taxonomy taxonomy.gg.qza --p-exclude mitochondria,chloroplast --o-filtered-table no-mitochondria-no-chloroplast.table.forward.all.gg.qza

qiime feature-table filter-samples --i-table no-mitochondria-no-chloroplast.table.forward.all.gg.qza --m-metadata-file sample-metadata_env.tsv --o-filtered-table no-mitochondria-no-chloroplast.table.forward.env.gg.qza

qiime feature-table filter-samples --i-table no-mitochondria-no-chloroplast.table.forward.env.gg.qza --m-metadata-file sample-metadata_env_water.tsv --o-filtered-table no-mitochondria-no-chloroplast.table.forward.env.gg_water.qza

qiime feature-classifier classify-consensus-blast --i-query table.forward.all.qza --i-reference-reads silva-138-99-seqs.qza --i-reference-taxonomy silva-138-99-tax.qza --output-dir blast_classification/

qiime metadata tabulate --m-input-file blast_classification/classification.qza --o-visualization blast-classification.qzv

qiime taxa collapse --i-table table.forward.all.qza --i-taxonomy blast_classification/classification.qza \
--p-level 7 --output-dir blast_taxtable/

qiime tools export --input-path blast_taxtable/collapsed_table.qza --output-path blast_taxtable/

biom convert -i blast_taxtable/feature-table.biom -o blast_taxtable/feature_table_collapsed.txt --to-tsv
