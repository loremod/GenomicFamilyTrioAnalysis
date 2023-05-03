@include "../functions.awk"
{
    #father                     mother                      child
    get_alleles(father,$10);    get_alleles(mother,$11);    get_alleles(child,$12);
    if(is_heterozygous(child))
        if(is_de_novo_variant(child[1],father,mother) || is_de_novo_variant(child[2],father,mother))
            print $0
}