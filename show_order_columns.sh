#CONTROLLA SE LE COLONNE DEI SAMPLE SONO NELL'ORDINE GIUSO
trios=$(cat my_trios)
source colors
for id in $trios
do
    # success "FILE case$id/called_var$id.vcf "
	# head case$id/called_var$id.vcf 
    success "FILE case$id/results$id.vcf "
	cat case$id/results$id.vcf | grep "#CHROM" | cut -f 10,11,12 
done