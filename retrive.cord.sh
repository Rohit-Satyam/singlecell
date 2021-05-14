while read a b; 
do 
mod=$(echo $a |tr '-' '_'); 
coord=$(grep $mod PlasmoDB-37_Pfalciparum3D7.gff | grep -w gene | awk {print "pf3d7_"$1":"$4"-"$5}' ); 
echo  -e $a"\t"$coord"\t"$b >> var_genes_with_coord.csv; done < kaust_vargene_list.tsv

## After this edit in excel the PF3D7 ----> Pf3D7 and prepare a bed file.

## samtools -hb -L var_genes.bed possorted.bam > var_gene_16D.bam

## Subsetting the bamswith theuniquely mapping reads obtained from STAR.

samtools view -H var_gene_16D.bam> 16D.header
samtools view var_gene_16D.bam | grep -w -f <(sed 's/.*/\\<&\\>/' var.genes) | grep -w -f 16D_unique_readID.txt > 16Dunique.sam ## because our bam containes read IDs for rifins, surfins etc too.
cat 16D.header 16Dunique.sam | samtools view -bh - > 16D_new.bam

python processing_for_genes.py -b 16D_new.bam -i 16D_new.bam.bai -c 16D_count.csv -p 16D_results.pdf


### References
## https://github.com/alexdobin/STAR/issues/460
## https://www.biostars.org/p/243020/#:~:text=%22too%20short%22%20means%20that%20the,of%20mates)%20should%20be%20mapped
## https://www.ecseq.com/support/ngs-snippets/how-to-extract-a-list-of-specific-read-IDs-from-a-BAM-file



