#!/bin/bash
source colors

###FUNCTIONS###
AD_prioritization(){
	grep -v "#" results$1.vcf | awk -F "\t" -f ../AD.awk >> called_var$1.vcf
}

AR_prioritization(){
	grep -v "#" results$1.vcf | awk -F "\t" -f ../AR.awk >> called_var$1.vcf
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

	cd ..
}

getCoverage(){
	info "---COMPUTING COVERAGE---"
	trios=$(cat my_trios)
	for id in $trios
	do
		info "Starting case ${id}..."
		bedtools genomecov -ibam case${id}/case${id}_father.bam -bg -trackline -trackopts 'name="father${id}"' -max 100 > coverage/fatherCov${id}.bg
		bedtools genomecov -ibam case${id}/case${id}_mother.bam -bg -trackline -trackopts 'name="mother${id}"' -max 100 > coverage/motherCov${id}.bg
		bedtools genomecov -ibam case${id}/case${id}_child.bam -bg -trackline -trackopts 'name="child${id}"' -max 100 > coverage/childCov${id}.bg
		success "Case ${id} done"
	done
	success "---COVERAGE COMPUTED---"
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
getCoverage #get coverage from every bam file


###GIVE EXECUTION PERMISSION###
#chmod +x ./pipeline.sh

###COMMAND TO EXECUTE###
#./pipeline.sh $(cat my_trios)
