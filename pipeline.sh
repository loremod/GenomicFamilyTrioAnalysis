#!/bin/bash
source colors

###FUNCTIONS###
AD_prioritization(){
	grep -P '\t0\/0.*?\t0\/0.*?\t0\/1.*' results$1.vcf >> called_var$1.vcf
	# cp called_var$1.vcf out_filtered
}

AR_prioritization(){
	grep -P "\t(0\/[12]).*?\t(0\/[12]).*?\t(1\/1|2\/2).*" results$1.vcf >> called_var$1.vcf
	# cp called_var$1.vcf out_filtered
}

compute_prioritization(){
	local inher_path=$2
	local im=`grep $1 $inher_path | cut -f 2`
    grep "#" results$1.vcf > called_var$1.vcf
	if [ "$im" = "AR" ]
	then
		info "->Sample $1 is AR"
		AR_prioritization $1
		success "->called_var$1.vcf created"
	elif [ "$im" = "AD" ]
	then
		info "->Sample $1 is AD"
		AD_prioritization $1
		success "->called_var$1.vcf created"
	else
		warning "Parameter for genetic model not recognized"
        return 
	fi
    grep "#" called_var$1.vcf > ${1}called_varTG.vcf
    bedtools intersect -a called_var$1.vcf -b /home/BCG2023_genomics_exam/exons16Padded_sorted.bed -u >> ${1}called_varTG.vcf
    success "->called_var$1.vcf created [FINAL VCF FILE]"
    # head ${1}called_varTG.vcf
}

fastqcStep(){
	#Quality controls for each .fq file
	for filename in *.fq
	do
		fastqc $filename
	done
}

alignmentStep(){
	for filename in *.fq
	do
		name=${filename%.*}
		bowtie2 -U $filename -x uni --rg-id $name --rg "SM:$name" | samtools view -Sb | samtools sort -o ${filename%.*}.bam
	done
    for filename in *.bam
	do
		samtools index $filename
	done
}

qualimapStep(){
    for filename in *.bam
	do
		qualimap bamqc -bam $filename -outdir ${filename%.*}
	done
}

pipeline(){
	mkdir case$1
	cd case$1
	cp ../../../BCG2023_genomics_exam/case$1* .
	cp ../../../BCG2023_genomics_exam/universe.fasta . 
    cp ../uni.* . 

	gunzip *gz

    info "Starting FastQC step..."
    fastqcStep
    success "FastQC step completed..."

    info "Starting alignment step..."
	alignmentStep
    success "Alignment step completed..."

    info "Starting qualimap step..."
    qualimapStep
    success "Qualimap step completed..."

    info "Starting multiqc step..."
	multiqc ./
    success "Multiqc step completed..."

    info "Starting freebayes step..."
	freebayes -f universe.fasta -m 20 -C 5 -Q 10 --min-coverage 10 case$1_child.bam case$1_mother.bam case$1_father.bam > results$1.vcf
    success "Freebayes step completed..."

    info "Starting variant prioritization step..."
	compute_prioritization $1 "../inher_mod"
    success "Variant prioritization step completed..."
}

###TESTING FUNCTIONS###
TestForVariantPrioritization(){
    for ID in "$@"; do
        cd case$1
        info "->Starting Variant prioritization of $ID"
        compute_prioritization $ID "./inher_mod"
        success "->$ID variant prioritization successfully completed"
        cd ..
    done
}

###MAIN###
bowtie2-build ../../BCG2023_genomics_exam/universe.fasta uni
for ID in "$@"; do
    info "---STARTING PIPELINE FOR ID = $ID---"
    pipeline $ID
    success "---$ID PIPELINE COMPLETED---"
done

###GIVE EXECUTION PERMISSION###
#chmod +x ./pipeline.sh

###COMMAND TO EXECUTE###
#./pipeline.sh $(cat my_trios)