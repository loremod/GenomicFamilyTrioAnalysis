function get_alleles(arr,str){gsub(/:[^\t]+/, "", str); split(str, arr, "/");}
function has_allele(who,allele){return (allele == who[1] || allele == who[2])}
function is_variant(allele){return allele!=0}
function is_de_novo(allele,father,mother){return (!has_allele(father,allele) && !has_allele(mother,allele))}
function is_de_novo_variant(allele,father,mother){return (is_variant(allele) && is_de_novo(allele,father,mother))}
function is_heterozygous(alleles){return (alleles[1]!=alleles[2])}
function is_homozygous(alleles){return (alleles[1]==alleles[2])}