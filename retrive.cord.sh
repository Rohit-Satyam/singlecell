while read a b; 
do 
mod=$(echo $a |tr '-' '_'); 
coord=$(grep $mod PlasmoDB-37_Pfalciparum3D7.gff | grep -w gene | awk -v var=$mod '{print var":"$4"-"$5}' ); 
echo  -e $a"\t"$coord"\t"$b >> var_genes_with_coord.csv; 
done < kaust_vargene_list.tsv
