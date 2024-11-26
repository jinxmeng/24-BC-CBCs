makeblastdb -in CS.gene.id.ffn -out CS.gene.id.ffn -dbtype nucl -parse_seqids
blastn -query BC/cds_from_genomic.fna -db CS.gene.id.ffn -evalue 1e-5 -max_target_seqs 5 -num_threads 10 -outfmt 6 -out BC.ffn.btn

makeblastdb -in CS.gene.id.faa -out CS.gene.id.faa -dbtype prot
blastp -query BC/protein.faa -db CS.gene.id.faa -evalue 1e-5 -max_target_seqs 5 -num_threads 10 -outfmt 6 -out BC.faa.btp
blastx -query BC/cds_from_genomic.fna -db CS.gene.id.faa -evalue 1e-5 -max_target_seqs 5 -num_threads 10 -outfmt 6 -out BC.ffn.btp
blastx -query BC/GCF_034355335.1_ASM3435533v1_genomic.fna -db CS.gene.id.faa -evalue 1e-5 -max_target_seqs 5 -num_threads 10 -outfmt 6 -out BC.fna.btp

blastp -query EC/protein.faa -db CS.gene.id.faa -evalue 1e-5 -max_target_seqs 5 -num_threads 10 -outfmt 6 -out EC.faa.btp
blastx -query EC/GCF_000005845.2_ASM584v2_genomic.fna -db CS.gene.id.faa -evalue 1e-5 -max_target_seqs 5 -num_threads 10 -outfmt 6 -out EC.fna.btp 
blastx -query EC/cds_from_genomic.fna -db CS.gene.id.faa -evalue 1e-5 -max_target_seqs 5 -num_threads 10 -outfmt 6 -out EC.ffn.btp

seqkit grep -p "CLOSPO_00312|acdA" CS.gene.id.faa > acdA.faa
cat BC.faa.btp | grep acdA | cut -f1 | seqkit grep -f - BC/protein.faa | seqkit replace -p " .*" -r "" | seqkit replace -p "^" -r "B.coccoides|" >> acdA.faa
cat EC.faa.btp | grep acdA | cut -f1 | seqkit grep -f - EC/protein.faa | seqkit replace -p " .*" -r "" | seqkit replace -p "^" -r "E.coli|" >> acdA.faa

cat BC.faa.btp | grep fldH | cut -f1 | seqkit grep -f - BC/protein.faa | seqkit replace -p " .*" -r "" | seqkit replace -p "^" -r "B.coccoides|" >> fldH.faa
cat EC.faa.btp | grep fldH | cut -f1 | seqkit grep -f - EC/protein.faa | seqkit replace -p " .*" -r "" | seqkit replace -p "^" -r "E.coli|" >> fldH.faa

cat BC.faa.btp | grep fldC | cut -f1 | seqkit grep -f - BC/protein.faa | seqkit replace -p " .*" -r "" | seqkit replace -p "^" -r "B.coccoides|" >> fldC.faa
cat EC.faa.btp | grep fldC | cut -f1 | seqkit grep -f - EC/protein.faa | seqkit replace -p " .*" -r "" | seqkit replace -p "^" -r "E.coli|" >> fldC.faa

muscle -align acdA.faa -output acdA.msa
trimal -in acdA.msa -out acdA.msa.trim
