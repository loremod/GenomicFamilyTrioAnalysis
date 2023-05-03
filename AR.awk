@include "../functions.awk"
{
    #father                     mother                      child
    get_alleles(father,$10);    get_alleles(mother,$11);    get_alleles(child,$12);
    if(is_variant(child[1]))
        if(is_heterozygous(father) && is_heterozygous(mother) && is_homozygous(child))
            print $0
}