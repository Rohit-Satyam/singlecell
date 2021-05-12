while read a b; 
do 
mod=$(echo $a |tr '-' '_'); 
coord=$(grep $mod PlasmoDB-37_Pfalciparum3D7.gff | grep -w gene | awk {print "pf3d7_"$1":"$4"-"$5}' ); 
echo  -e $a"\t"$coord"\t"$b >> var_genes_with_coord.csv; done < kaust_vargene_list.tsv

## After this edit in excel the PF3D7 ----> Pf3D7 and prepare a bed file.

## samtools -hb -L var_genes.bed possorted.bam > var_gene_16D.bam

