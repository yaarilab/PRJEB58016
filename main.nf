$HOSTNAME = ""
params.outdir = 'results'  

evaluate(new File("${params.projectDir}/nextflow_header.config"))
params.metadata.metadata = "${params.projectDir}/tools.json"


if (!params.mate){params.mate = ""} 
if (!params.reads){params.reads = ""} 

Channel.value(params.mate).into{g_1_mate_g_63;g_1_mate_g_16;g_1_mate_g22_14;g_1_mate_g22_12;g_1_mate_g22_10;g_1_mate_g37_9;g_1_mate_g52_8;g_1_mate_g52_1;g_1_mate_g52_0;g_1_mate_g28_15;g_1_mate_g28_19;g_1_mate_g28_12;g_1_mate_g38_11;g_1_mate_g38_9;g_1_mate_g38_12}
if (params.reads){
Channel
	.fromFilePairs( params.reads , size: params.mate == "single" ? 1 : params.mate == "pair" ? 2 : params.mate == "triple" ? 3 : params.mate == "quadruple" ? 4 : -1 )
	.ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
	.set{g_64_reads_g_63}
 } else {  
	g_64_reads_g_63 = Channel.empty()
 }



process unizp {

input:
 set val(name),file(reads) from g_64_reads_g_63
 val mate from g_1_mate_g_63

output:
 set val(name),file("*.fastq")  into g_63_reads0_g38_11

script:

if(mate=="pair"){
	readArray = reads.toString().split(' ')	
	R1 = readArray[0]
	R2 = readArray[1]
	
	"""
	case "$R1" in
	*.gz | *.tgz ) 
	        gunzip -c $R1 > R1.fastq
	        ;;
	*)
	        cp $R1 ./R1.fastq
	        echo "$R1 not gzipped"
	        ;;
	esac
	
	case "$R2" in
	*.gz | *.tgz ) 
	        gunzip -c $R2 > R2.fastq
	        ;;
	*)
	        cp $R2 ./R2.fastq
	        echo "$R2 not gzipped"
	        ;;
	esac
	"""
}else{
	"""
	case "$reads" in
	*.gz | *.tgz ) 
	        gunzip -c $reads > R1.fastq
	        ;;
	*)
	        cp $reads ./R1.fastq
	        echo "$reads not gzipped"
	        ;;
	esac
	"""
}
}


process Mask_Primer_MaskPrimers {

input:
 val mate from g_1_mate_g38_11
 set val(name),file(reads) from g_63_reads0_g38_11

output:
 set val(name), file("*_primers-pass.fast*")  into g38_11_reads0_g_16
 set val(name), file("*_primers-fail.fast*") optional true  into g38_11_reads_failed11
 set val(name), file("MP_*")  into g38_11_logFile2_g38_9
 set val(name),file("out*")  into g38_11_logFile33

script:
method = params.Mask_Primer_MaskPrimers.method
barcode_field = params.Mask_Primer_MaskPrimers.barcode_field
primer_field = params.Mask_Primer_MaskPrimers.primer_field
barcode = params.Mask_Primer_MaskPrimers.barcode
revpr = params.Mask_Primer_MaskPrimers.revpr
mode = params.Mask_Primer_MaskPrimers.mode
failed = params.Mask_Primer_MaskPrimers.failed
fasta = params.Mask_Primer_MaskPrimers.fasta
nproc = params.Mask_Primer_MaskPrimers.nproc
maxerror = params.Mask_Primer_MaskPrimers.maxerror
umi_length = params.Mask_Primer_MaskPrimers.umi_length
start = params.Mask_Primer_MaskPrimers.start
extract_length = params.Mask_Primer_MaskPrimers.extract_length
maxlen = params.Mask_Primer_MaskPrimers.maxlen
skiprc = params.Mask_Primer_MaskPrimers.skiprc
R1_primers = params.Mask_Primer_MaskPrimers.R1_primers
R2_primers = params.Mask_Primer_MaskPrimers.R2_primers
//* @style @condition:{method="score",umi_length,start,maxerror}{method="extract",umi_length,start},{method="align",maxerror,maxlen,skiprc}, {method="extract",start,extract_length} @array:{method,barcode_field,primer_field,barcode,revpr,mode,maxerror,umi_length,start,extract_length,maxlen,skiprc} @multicolumn:{method,barcode_field,primer_field,barcode,revpr,mode,failed,nproc,maxerror,umi_length,start,extract_length,maxlen,skiprc}

method = (method.collect().size==2) ? method : [method[0],method[0]]
barcode_field = (barcode_field.collect().size==2) ? barcode_field : [barcode_field[0],barcode_field[0]]
primer_field = (primer_field.collect().size==2) ? primer_field : [primer_field[0],primer_field[0]]
barcode = (barcode.collect().size==2) ? barcode : [barcode[0],barcode[0]]
revpr = (revpr.collect().size==2) ? revpr : [revpr[0],revpr[0]]
mode = (mode.collect().size==2) ? mode : [mode[0],mode[0]]
maxerror = (maxerror.collect().size==2) ? maxerror : [maxerror[0],maxerror[0]]
umi_length = (umi_length.collect().size==2) ? umi_length : [umi_length[0],umi_length[0]]
start = (start.collect().size==2) ? start : [start[0],start[0]]
extract_length = (extract_length.collect().size==2) ? extract_length : [extract_length[0],extract_length[0]]
maxlen = (maxlen.collect().size==2) ? maxlen : [maxlen[0],maxlen[0]]
skiprc = (skiprc.collect().size==2) ? skiprc : [skiprc[0],skiprc[0]]
failed = (failed=="true") ? "--failed" : ""
fasta = (fasta=="true") ? "--fasta" : ""
def args_values = [];
[method,barcode_field,primer_field,barcode,revpr,mode,maxerror,umi_length,start,extract_length,maxlen,skiprc].transpose().each { m,bf,pf,bc,rp,md,mr,ul,s,el,ml,sk -> {
    
    if(m=="align"){
        s = ""
    }else{
        if(bc=="false"){
            s = "--start ${s}"
        }else{
            s = s + ul
            s = "--start ${s}"
        }
    }
    
    el = (m=="extract") ? "--len ${el}" : ""
    mr = (m=="extract") ? "" : "--maxerror ${mr}" 
    ml = (m=="align") ? "--maxlen ${ml}" : "" 
    sk = (m=="align" && sk=="true") ? "--skiprc" : "" 
    
    PRIMER_FIELD = "${pf}"
    
    // all
    bf = (bf=="") ? "" : "--bf ${bf}"
    pf = (pf=="") ? "" : "--pf ${pf}"
    bc = (bc=="false") ? "" : "--barcode"
    rp = (rp=="false") ? "" : "--revpr"
    args_values.add("${m} --mode ${md} ${bc} ${rp} ${mr} ${s} ${el} ${ml} ${sk} ${pf} ${bf}")
    
    
}}

readArray = reads.toString().split(' ')
if(mate=="pair"){
	args_1 = args_values[0]
	args_2 = args_values[1]
	
  


	R1 = readArray[0]
	R2 = readArray[1]
	
	R1_primers = (method[0]=="extract") ? "" : "-p ${R1_primers}"
	R2_primers = (method[1]=="extract") ? "" : "-p ${R2_primers}"
	
	
	"""
	
	MaskPrimers.py ${args_1} -s ${R1} ${R1_primers} --log MP_R1_${name}.log  --nproc ${nproc} ${failed} ${fasta} 2>&1 | tee -a out_${R1}_MP.log & \
	MaskPrimers.py ${args_2} -s ${R2} ${R2_primers} --log MP_R2_${name}.log  --nproc ${nproc} ${failed} ${fasta} 2>&1 | tee -a out_${R1}_MP.log & \
	wait
	"""
}else{
	args_1 = args_values[0]
	
	R1_primers = (method[0]=="extract") ? "" : "-p ${R1_primers}"
	
	R1 = readArray[0]

	"""
	echo -e "Assuming inputs for R1\n"
	
	MaskPrimers.py ${args_1} -s ${reads} ${R1_primers} --log MP_${name}.log  --nproc ${nproc} ${failed} ${fasta} 2>&1 | tee -a out_${R1}_MP.log
	"""
}

}


process Cluster_UMI {

input:
 set val(name),file(reads) from g38_11_reads0_g_16
 val mate from g_1_mate_g_16

output:
 set val(name),file("*_umi-pass.fastq") optional true  into g_16_reads0_g52_0

script:
umi_field = params.Cluster_UMI.umi_field
umi_cluster_script = params.Cluster_UMI.umi_cluster_script

readArray = reads.toString().split(' ')	

R1 = readArray[0]
R2 = readArray[1]



filename_R1=file(R1).getSimpleName()
filename_R2=file(R2).getSimpleName()



"""
#!/bin/bash

cat ${R1} ${R2} | awk '{if (NR%4==1) {split(\$0,a,"${umi_field}="); print a[2]}}' > pair_pass.umi

python3 ${umi_cluster_script} -i pair_pass.umi -o pair_pass.umi.convert

awk -F'\t' 'NR==FNR && NR>1 {umi[\$1]=\$2;} NR!=FNR {if(/UMI=/){split(\$0,a,"${umi_field}=");split(a[2],b,"|"); print a[1]"UMI="b[1];} else {print}}' pair_pass.umi.convert ${R1} > ${R1}'_umi-pass.fastq'
awk -F'\t' 'NR==FNR && NR>1 {umi[\$1]=\$2;} NR!=FNR {if(/UMI=/){split(\$0,a,"${umi_field}=");split(a[2],b,"|"); print a[1]"UMI="b[1] , b[2];} else {print}}' pair_pass.umi.convert ${R2} > ${R2}'_umi-pass.fastq'
"""


}


process Align_Sets_align_sets {

input:
 set val(name),file(reads) from g_16_reads0_g52_0
 val mate from g_1_mate_g52_0

output:
 set val(name),file("*_align-pass.fastq")  into g52_0_reads0_g22_10
 set val(name), file("AS_*")  into g52_0_logFile1_g52_1
 set val(name),file("*_align-fail.fastq") optional true  into g52_0_reads_failed22
 set val(name), file("out*") optional true  into g52_0_logFile3_g61_0

script:
method = params.Align_Sets_align_sets.method
bf = params.Align_Sets_align_sets.bf
div = params.Align_Sets_align_sets.div
failed = params.Align_Sets_align_sets.failed
nproc = params.Align_Sets_align_sets.nproc

muscle_exec = params.Align_Sets_align_sets.muscle_exec

offset_table = params.Align_Sets_align_sets.offset_table
pf = params.Align_Sets_align_sets.pf
mode = params.Align_Sets_align_sets.mode

primer_file = params.Align_Sets_align_sets.primer_file
reverse = params.Align_Sets_align_sets.reverse

//* @style @condition:{method="muscle",muscle_exec}, {method="offset",offset_table,pf,mode}, {method="table",muscle_exec,primer_file,reverse} @multicolumn:{method,bf,div,nproc},{offset,pf,mode}, {primer_file,reverse}


readArray = reads.toString().split(' ')	

reverse_arg = (reverse=="false") ? "" : "--reverse"
div_arg = (div=="false") ? "" : "--div"
failed_arg = (failed=="true") ? "--failed" : "" 
bf = "--bf ${bf}"

primer_file_argv = ""

if(method=="offset"){
	pf = "--pf ${pf}"
	mode = "--mode ${mode}"
	offset_table_argv = "-d ${offset_table}"
	muscle_exec_argv = ""
}else{
	pf = ""
	mode = ""
	offset_table_argv = ""
	muscle_exec_argv = "--exec ${muscle_exec}"
	
	if(method=="table"){
		primer_file_argv = "-p ${primer_file}"
	}
}

if(mate=="pair"){
	R1 = readArray[0]
	R2 = readArray[1]
	
	
	"""
	AlignSets.py ${method} -s ${R1} ${bf} ${muscle_exec_argv} ${div_arg} ${reverse_arg} ${failed_arg} ${pf} ${offset_table_argv} ${mode} ${primer_file_argv} --log AS_R1_${name}.log --nproc ${nproc} >> out_${R1}_AS.log
	AlignSets.py ${method} -s ${R2} ${bf} ${muscle_exec_argv} ${div_arg} ${reverse_arg} ${failed_arg} ${pf} ${offset_table_argv} ${mode} ${primer_file_argv} --log AS_R2_${name}.log --nproc ${nproc} >> out_${R1}_AS.log
	"""
	
}else{
	R1 = readArray[0]
	"""
	AlignSets.py ${method} -s ${R1} ${bf} ${muscle_exec_argv} ${div_arg} ${reverse_arg} ${failed_arg} ${pf} ${offset_table_argv} ${mode} ${primer_file_argv} --log AS_R1_${name}.log --nproc ${nproc} >> out_${R1}_AS.log
	"""
}

}


process Align_Sets_parse_log_AS {

input:
 set val(name), file(log_file) from g52_0_logFile1_g52_1
 val mate from g_1_mate_g52_1

output:
 set val(name), file("*.tab")  into g52_1_logFile0_g52_8

script:
field_to_parse = params.Align_Sets_parse_log_AS.field_to_parse
readArray = log_file.toString()	

"""
ParseLog.py -l ${readArray}  -f ${field_to_parse}
"""

}


process Align_Sets_report_Align_Sets {

input:
 set val(name), file(log_files) from g52_1_logFile0_g52_8
 val matee from g_1_mate_g52_8

output:
 file "*.rmd"  into g52_8_rMarkdown0_g52_9


shell:

if(matee=="pair"){
	readArray = log_files.toString().split(' ')	
	R1 = readArray[0]
	R2 = readArray[1]

	'''
	#!/usr/bin/env perl
	
	
	my $script = <<'EOF';
	
	
	```{R, message=FALSE, echo=FALSE, results="hide"}
	# Setup
	library(prestor)
	library(knitr)
	library(captioner)
	
	plot_titles <- plot_titles<- c("Read 1", "Read 2")
	if (!exists("tables")) { tables <- captioner(prefix="Table") }
	if (!exists("figures")) { figures <- captioner(prefix="Figure") }
	figures("align_size", 
	        paste("Histogram of UMI read group sizes (reads per UMI) for", plot_titles[1], "(top) and", plot_titles[2],
	        "(bottom). The x-axis indicates the number of reads in a UMI group and the y-axis is the 
	        number of UMI groups with that size."))
	
	```
	
	```{r, echo=FALSE}
	align_log_1 <- loadLogTable(file.path(".", "!{R1}"))
	align_log_2 <- loadLogTable(file.path(".", "!{R2}"))
	```
	
	# Multiple Alignment of UMI Read Groups
	
	Reads sharing the same UMI are multiple aligned using the muscle wrapper in the 
	AlignSets tool.
	
	## Reads per UMI
	
	```{r, echo=FALSE}
	plotAlignSets(align_log_1, align_log_2, style="size", sizing="figure")
	```
	
	`r figures("align_size")`

	
	EOF
	
	open OUT, ">!{name}.rmd";
	print OUT $script;
	close OUT;
	
	'''

}else{
	
	readArray = log_files.toString().split(' ')
	R1 = readArray[0]
	
	'''
	#!/usr/bin/env perl
	
	
	my $script = <<'EOF';
	
	
	```{R, message=FALSE, echo=FALSE, results="hide"}
	# Setup
	library(prestor)
	library(knitr)
	library(captioner)
	
	if (!exists("tables")) { tables <- captioner(prefix="Table") }
	if (!exists("figures")) { figures <- captioner(prefix="Figure") }
	figures("align_size", 
	        "Histogram of UMI read group sizes (reads per UMI). 
	        The x-axis indicates the number of reads in a UMI group and the y-axis is the 
	        number of UMI groups with that size.")
	```
	
	```{r, echo=FALSE}
	align_log <- loadLogTable(file.path(".","!{R1}"))
	```
	
	# Multiple Alignment of UMI Read Groups
	
	Reads sharing the same UMI are multiple aligned using the muscle wrapper in the 
	AlignSets tool.
	
	## Reads per UMI
	
	```{r, echo=FALSE}
	plotAlignSets(align_log, style="size", sizing="figure")
	```
	
	`r figures("align_size")`
	
	EOF
	
	open OUT, ">!{name}.rmd";
	print OUT $script;
	close OUT;
	
	'''
}

}


process Align_Sets_render_rmarkdown {

input:
 file rmk from g52_8_rMarkdown0_g52_9

output:
 file "*.html"  into g52_9_outputFileHTML00
 file "*csv" optional true  into g52_9_csvFile11

"""

#!/usr/bin/env Rscript 

rmarkdown::render("${rmk}", clean=TRUE, output_format="html_document", output_dir=".")

"""
}

boolean isCollectionOrArray_bc(object) {    
    [Collection, Object[]].any { it.isAssignableFrom(object.getClass()) }
}

def args_creator_bc(barcode_field, primer_field, act, copy_field, mincount, minqual, minfreq, maxerror, prcons, maxgap, maxdiv, dep){
	def args_values;
    if(isCollectionOrArray_bc(barcode_field) || isCollectionOrArray_bc(primer_field) || isCollectionOrArray_bc(copy_field) || isCollectionOrArray_bc(mincount) || isCollectionOrArray_bc(minqual) || isCollectionOrArray_bc(minfreq) || isCollectionOrArray_bc(maxerror) || isCollectionOrArray_bc(prcons) || isCollectionOrArray_bc(maxgap) || isCollectionOrArray_bc(maxdiv) || isCollectionOrArray_bc(dep)){
    	primer_field = (isCollectionOrArray_bc(primer_field)) ? primer_field : [primer_field,primer_field]
    	act = (isCollectionOrArray_bc(act)) ? act : [act,act]
    	copy_field = (isCollectionOrArray_bc(copy_field)) ? copy_field : [copy_field,copy_field]
    	mincount = (isCollectionOrArray_bc(mincount)) ? mincount : [mincount,mincount]
    	minqual = (isCollectionOrArray_bc(minqual)) ? minqual : [minqual,minqual]
    	minfreq = (isCollectionOrArray_bc(minfreq)) ? minfreq : [minfreq,minfreq]
    	maxerror = (isCollectionOrArray_bc(maxerror)) ? maxerror : [maxerror,maxerror]
    	prcons = (isCollectionOrArray_bc(prcons)) ? prcons : [prcons,prcons]
    	maxgap = (isCollectionOrArray_bc(maxgap)) ? maxgap : [maxgap,maxgap]
    	maxdiv = (isCollectionOrArray_bc(maxdiv)) ? maxdiv : [maxdiv,maxdiv]
    	dep = (isCollectionOrArray_bc(dep)) ? dep : [dep,dep]
    	args_values = []
        [barcode_field,primer_field,act,copy_field,mincount,minqual,minfreq,maxerror,prcons,maxgap,maxdiv,dep].transpose().each { bf,pf,a,cf,mc,mq,mf,mr,pc,mg,md,d -> {
            bf = (bf=="") ? "" : "--bf ${bf}"
            pf = (pf=="") ? "" : "--pf ${pf}" 
            a = (a=="none") ? "" : "--act ${a}" 
            cf = (cf=="") ? "" : "--cf ${cf}" 
            mr = (mr=="none") ? "" : "--maxerror ${mr}" 
            pc = (pc=="none") ? "" : "--prcons ${pc}" 
            mg = (mg=="none") ? "" : "--maxgap ${mg}" 
            md = (md=="none") ? "" : "--maxdiv ${md}" 
            d = (d=="true") ? "--dep" : "" 
            args_values.add("${bf} ${pf} ${a} ${cf} -n ${mc} -q ${mq} --freq ${mf} ${mr} ${pc} ${mg} ${md} ${d}")
        }}
    }else{
        barcode_field = (barcode_field=="") ? "" : "--bf ${barcode_field}"
        primer_field = (primer_field=="") ? "" : "--pf ${primer_field}" 
        act = (act=="none") ? "" : "--act ${act}" 
        copy_field = (copy_field=="") ? "" : "--cf ${copy_field}" 
        maxerror = (maxerror=="none") ? "" : "--maxerror ${maxerror}" 
        prcons = (prcons=="none") ? "" : "--prcons ${prcons}" 
        maxgap = (maxgap=="none") ? "" : "--maxgap ${maxgap}" 
        maxdiv = (maxdiv=="none") ? "" : "--maxdiv ${maxdiv}" 
        dep = (dep=="true") ? "--dep" : "" 
        args_values = "${barcode_field} ${primer_field} ${act} ${copy_field} -n ${mincount} -q ${minqual} --freq ${minfreq} ${maxerror} ${prcons} ${maxgap} ${maxdiv} ${dep}"
    }
    return args_values
}


process Build_Consensus_build_consensus {

input:
 set val(name),file(reads) from g52_0_reads0_g22_10
 val mate from g_1_mate_g22_10

output:
 set val(name),file("*_consensus-pass.fastq")  into g22_10_reads0_g37_9
 set val(name),file("BC*")  into g22_10_logFile1_g22_12
 set val(name),file("out*")  into g22_10_logFile2_g61_0

script:
failed = params.Build_Consensus_build_consensus.failed
nproc = params.Build_Consensus_build_consensus.nproc
barcode_field = params.Build_Consensus_build_consensus.barcode_field
primer_field = params.Build_Consensus_build_consensus.primer_field
act = params.Build_Consensus_build_consensus.act
copy_field = params.Build_Consensus_build_consensus.copy_field
mincount = params.Build_Consensus_build_consensus.mincount
minqual = params.Build_Consensus_build_consensus.minqual
minfreq = params.Build_Consensus_build_consensus.minfreq
maxerror = params.Build_Consensus_build_consensus.maxerror
prcons = params.Build_Consensus_build_consensus.prcons
maxgap = params.Build_Consensus_build_consensus.maxgap
maxdiv = params.Build_Consensus_build_consensus.maxdiv
dep = params.Build_Consensus_build_consensus.dep
//* @style @condition:{act="none",},{act="min",copy_field},{act="max",copy_field},{act="sum",copy_field},{act="set",copy_field},{act="majority",copy_field} @array:{barcode_field,primer_field,act,copy_field,mincount,minqual,minfreq,maxerror,prcons,maxgap,maxdiv,dep} @multicolumn:{failed,nproc},{barcode_field,primer_field,act,copy_field}, {mincount,minqual,minfreq,maxerror,prcons,maxgap,maxdiv,dep}

args_values_bc = args_creator_bc(barcode_field, primer_field, act, copy_field, mincount, minqual, minfreq, maxerror, prcons, maxgap, maxdiv, dep)

// args 
if(isCollectionOrArray_bc(args_values_bc)){
	args_1 = args_values_bc[0]
	args_2 = args_values_bc[1]
}else{
	args_1 = args_values_bc
	args_2 = args_values_bc
}

failed = (failed=="true") ? "--failed" : "" 


if(mate=="pair"){
	// files
	readArray = reads.toString().split(' ')	
	R1 = readArray[0]
	R2 = readArray[1]
	
	"""
	BuildConsensus.py --version
	BuildConsensus.py -s $R1 ${args_1} --log BC_${name}_R1.log ${failed} --nproc ${nproc} >> out_${R1}_BC.log
	BuildConsensus.py -s $R2 ${args_2} --log BC_${name}_R2.log ${failed} --nproc ${nproc} >> out_${R1}_BC.log
	"""
}else{
	"""
	BuildConsensus.py -s $reads ${args_1} --outname ${name} --log BC_${name}.log ${failed} --nproc ${nproc} >> out_${R1}_BC.log
	"""
}


}


process Pair_Sequence_post_consensus_pair_seq {

input:
 set val(name),file(reads) from g22_10_reads0_g37_9
 val mate from g_1_mate_g37_9

output:
 set val(name),file("*_pair-pass.fastq")  into g37_9_reads0_g28_12
 set val(name),file("out*")  into g37_9_logFile1_g61_0

script:
coord = params.Pair_Sequence_post_consensus_pair_seq.coord
act = params.Pair_Sequence_post_consensus_pair_seq.act
copy_fields_1 = params.Pair_Sequence_post_consensus_pair_seq.copy_fields_1
copy_fields_2 = params.Pair_Sequence_post_consensus_pair_seq.copy_fields_2
failed = params.Pair_Sequence_post_consensus_pair_seq.failed
nproc = params.Pair_Sequence_post_consensus_pair_seq.nproc

if(mate=="pair"){
	
	act = (act=="none") ? "" : "--act ${act}"
	failed = (failed=="true") ? "--failed" : "" 
	copy_fields_1 = (copy_fields_1=="") ? "" : "--1f ${copy_fields_1}" 
	copy_fields_2 = (copy_fields_2=="") ? "" : "--2f ${copy_fields_2}"
	
	readArray = reads.toString().split(' ')	
	R1 = readArray[0]
	R2 = readArray[1]
	"""
	PairSeq.py -1 ${R1} -2 ${R2} ${copy_fields_1} ${copy_fields_2} --coord ${coord} ${act} ${failed} >> out_${R1}_PS.log
	"""
}else{
	
	"""
	echo -e 'PairSeq works only on pair-end reads.'
	"""
}


}


process Assemble_pairs_assemble_pairs {

input:
 set val(name),file(reads) from g37_9_reads0_g28_12
 val mate from g_1_mate_g28_12

output:
 set val(name),file("*_assemble-pass.f*")  into g28_12_reads0_g_65
 set val(name),file("AP_*")  into g28_12_logFile1_g28_15
 set val(name),file("*_assemble-fail.f*") optional true  into g28_12_reads_failed22
 set val(name),file("out*")  into g28_12_logFile33

script:
method = params.Assemble_pairs_assemble_pairs.method
coord = params.Assemble_pairs_assemble_pairs.coord
rc = params.Assemble_pairs_assemble_pairs.rc
head_fields_R1 = params.Assemble_pairs_assemble_pairs.head_fields_R1
head_fields_R2 = params.Assemble_pairs_assemble_pairs.head_fields_R2
failed = params.Assemble_pairs_assemble_pairs.failed
fasta = params.Assemble_pairs_assemble_pairs.fasta
nproc = params.Assemble_pairs_assemble_pairs.nproc
alpha = params.Assemble_pairs_assemble_pairs.alpha
maxerror = params.Assemble_pairs_assemble_pairs.maxerror
minlen = params.Assemble_pairs_assemble_pairs.minlen
maxlen = params.Assemble_pairs_assemble_pairs.maxlen
scanrev = params.Assemble_pairs_assemble_pairs.scanrev
minident = params.Assemble_pairs_assemble_pairs.minident
evalue = params.Assemble_pairs_assemble_pairs.evalue
maxhits = params.Assemble_pairs_assemble_pairs.maxhits
fill = params.Assemble_pairs_assemble_pairs.fill
aligner = params.Assemble_pairs_assemble_pairs.aligner
// align_exec = params.Assemble_pairs_assemble_pairs.// align_exec
// dbexec = params.Assemble_pairs_assemble_pairs.// dbexec
gap = params.Assemble_pairs_assemble_pairs.gap
usearch_version = params.Assemble_pairs_assemble_pairs.usearch_version
assemble_reference = params.Assemble_pairs_assemble_pairs.assemble_reference
head_seqeunce_file = params.Assemble_pairs_assemble_pairs.head_seqeunce_file
//* @style @condition:{method="align",alpha,maxerror,minlen,maxlen,scanrev}, {method="sequential",alpha,maxerror,minlen,maxlen,scanrev,ref_file,minident,evalue,maxhits,fill,aligner,align_exec,dbexec} {method="reference",ref_file,minident,evalue,maxhits,fill,aligner,align_exec,dbexec} {method="join",gap} @multicolumn:{method,coord,rc,head_fields_R1,head_fields_R2,failed,nrpoc,usearch_version},{alpha,maxerror,minlen,maxlen,scanrev}, {ref_file,minident,evalue,maxhits,fill,aligner,align_exec,dbexec}, {gap} 

// args
coord = "--coord ${coord}"
rc = "--rc ${rc}"
head_fields_R1 = (head_fields_R1!="") ? "--1f ${head_fields_R1}" : ""
head_fields_R2 = (head_fields_R2!="") ? "--2f ${head_fields_R2}" : ""
failed = (failed=="false") ? "" : "--failed"
fasta = (fasta=="false") ? "" : "--fasta"
nproc = "--nproc ${nproc}"

scanrev = (scanrev=="false") ? "" : "--scanrev"
fill = (fill=="false") ? "" : "--fill"

// align_exec = (align_exec!="") ? "--exec ${align_exec}" : ""
// dbexec = (dbexec!="") ? "--dbexec ${dbexec}" : ""


ref_file = (assemble_reference!='') ? "-r ${assemble_reference}" : ""



args = ""

if(method=="align"){
	args = "--alpha ${alpha} --maxerror ${maxerror} --minlen ${minlen} --maxlen ${maxlen} ${scanrev}"
}else{
	if(method=="sequential"){
		args = "--alpha ${alpha} --maxerror ${maxerror} --minlen ${minlen} --maxlen ${maxlen} ${scanrev} ${ref_file} --minident ${minident} --evalue ${evalue} --maxhits ${maxhits} ${fill} --aligner ${aligner}"
	}else{
		if(method=="reference"){
			args = "${ref_file} --minident ${minident} --evalue ${evalue} --maxhits ${maxhits} ${fill} --aligner ${aligner}"
		}else{
			args = "--gap ${gap}"
		}
	}
}


readArray = reads.toString().split(' ')	


if(mate=="pair"){
	R1 = readArray[0]
	R2 = readArray[1]
	
	if(R1.contains("."+head_seqeunce_file)){
		R1 = readArray[0]
		R2 = readArray[1]
	}else{
		R2 = readArray[0]
		R1 = readArray[1]
	}
	
	"""
	if [ "${method}" != "align" ]; then
		if  [ "${aligner}" == "usearch" ]; then
			wget -q --show-progress --no-check-certificate https://drive5.com/downloads/usearch${usearch_version}_i86linux32.gz
			gunzip usearch${usearch_version}_i86linux32.gz
			chmod +x usearch${usearch_version}_i86linux32
			mv usearch${usearch_version}_i86linux32 /usr/local/bin/usearch2
			align_exec="--exec /usr/local/bin/usearch2"
			dbexec="--dbexec /usr/local/bin/usearch2"
		else
			align_exec="--exec /usr/local/bin/blastn"
			dbexec="--dbexec /usr/local/bin/makeblastdb"
		fi
	else
		align_exec=""
		dbexec=""
	fi

	AssemblePairs.py ${method} -1 ${R1} -2 ${R2} ${coord} ${rc} ${head_fields_R1} ${head_fields_R2} ${args} \$align_exec \$dbexec ${fasta} ${failed} --log AP_${name}.log ${nproc}  2>&1 | tee out_${R1}_AP.log
	"""

}else{
	
	"""
	echo -e 'AssemblePairs works only on pair-end reads.'
	"""
}

}


process vdjbase_input {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /${chain}$/) "reads/$filename"}
input:
 set val(name),file(reads) from g28_12_reads0_g_65

output:
 file "${chain}"  into g_65_germlineDb00

script:
chain = params.vdjbase_input.chain

"""
mkdir ${chain}
mv ${reads} ${chain}/${name}.fasta
"""

}


process Assemble_pairs_parse_log_AP {

input:
 set val(name),file(log_file) from g28_12_logFile1_g28_15
 val mate from g_1_mate_g28_15

output:
 set val(name),file("*.tab")  into g28_15_logFile0_g28_19, g28_15_logFile0_g28_25

script:
field_to_parse = params.Assemble_pairs_parse_log_AP.field_to_parse
readArray = log_file.toString()	

"""
ParseLog.py -l ${readArray}  -f ${field_to_parse}
"""


}


process Assemble_pairs_report_assemble_pairs {

input:
 set val(name),file(log_files) from g28_15_logFile0_g28_19
 val matee from g_1_mate_g28_19

output:
 file "*.rmd"  into g28_19_rMarkdown0_g28_25



shell:

if(matee=="pair"){
	readArray = log_files.toString().split(' ')
	assemble = readArray[0]
	'''
	#!/usr/bin/env perl
	
	
	my $script = <<'EOF';
	
	```{r, message=FALSE, echo=FALSE, results="hide"}
	# Setup
	library(prestor)
	library(knitr)
	library(captioner)
	
	if (!exists("tables")) { tables <- captioner(prefix="Table") }
	if (!exists("figures")) { figures <- captioner(prefix="Figure") }
	figures("assemble_length", "Histogram showing the distribution assembled sequence lengths in 
	                            nucleotides for the Align step (top) and Reference step (bottom).")
	figures("assemble_overlap", "Histogram showing the distribution of overlapping nucleotides between 
	                             mate-pairs for the Align step (top) and Reference step (bottom).
	                             Negative values for overlap indicate non-overlapping mate-pairs
	                             with the negative value being the number of gap characters between
	                             the ends of the two mate-pairs.")
	figures("assemble_error", "Histograms showing the distribution of paired-end assembly error 
	                           rates for the Align step (top) and identity to the reference germline 
	                           for the Reference step (bottom).")
	figures("assemble_pvalue", "Histograms showing the distribution of significance scores for 
	                            paired-end assemblies. P-values for the Align mode are shown in the top
	                            panel. E-values from the Reference step's alignment against the 
	                            germline sequences are shown in the bottom panel for both input files
	                            separately.")
	```
	
	```{r, echo=FALSE, warning=FALSE}
	assemble_log <- loadLogTable(file.path(".", "!{assemble}"))
	
	# Subset to align and reference logs
	align_fields <- c("ERROR", "PVALUE")
	ref_fields <- c("REFID", "GAP", "EVALUE1", "EVALUE2", "IDENTITY")
	align_log <- assemble_log[!is.na(assemble_log$ERROR), !(names(assemble_log) %in% ref_fields)]
	ref_log <- assemble_log[!is.na(assemble_log$REFID), !(names(assemble_log) %in% align_fields)]
	
	# Build log set
	assemble_list <- list()
	if (nrow(align_log) > 0) { assemble_list[["Align"]] <- align_log }
	if (nrow(ref_log) > 0) { assemble_list[["Reference"]] <- ref_log }
	plot_titles <- names(assemble_list)
	```
	
	# Paired-End Assembly
	
	Assembly of paired-end reads is performed using the AssemblePairs tool which 
	determines the read overlap in two steps. First, de novo assembly is attempted 
	using an exhaustive approach to identify all possible overlaps between the 
	two reads with alignment error rates and p-values below user-defined thresholds. 
	This method is denoted as the `Align` method in the following figures. 
	Second, those reads failing the first stage of de novo assembly are then 
	mapped to the V-region reference sequences to create a full length sequence, 
	padding with Ns, for any amplicons that have insufficient overlap for 
	de novo assembly. This second stage is referred to as the `Reference` step in the
	figures below.
	
	## Assembled sequence lengths
	
	```{r, echo=FALSE, warning=FALSE}
	plot_params <- list(titles=plot_titles, style="length", sizing="figure")
	do.call(plotAssemblePairs, c(assemble_list, plot_params))
	```
	
	`r figures("assemble_length")`
	
	```{r, echo=FALSE, warning=FALSE}
	plot_params <- list(titles=plot_titles, style="overlap", sizing="figure")
	do.call(plotAssemblePairs, c(assemble_list, plot_params))
	```
	
	`r figures("assemble_overlap")`
	
	## Alignment error rates and significance
	
	```{r, echo=FALSE, warning=FALSE}
	plot_params <- list(titles=plot_titles, style="error", sizing="figure")
	do.call(plotAssemblePairs, c(assemble_list, plot_params))
	```
	
	`r figures("assemble_error")`
	
	```{r, echo=FALSE, warning=FALSE}
	plot_params <- list(titles=plot_titles, style="pvalue", sizing="figure")
	do.call(plotAssemblePairs, c(assemble_list, plot_params))
	```

	`r figures("assemble_pvalue")`

	EOF
	
	open OUT, ">AP_!{name}.rmd";
	print OUT $script;
	close OUT;
	
	'''

}else{
	
	"""
	echo -e 'AssemblePairs works only on pair-end reads.'
	"""
}
}


process Assemble_pairs_presto_render_rmarkdown {

input:
 file rmk from g28_19_rMarkdown0_g28_25
 file log_file from g28_15_logFile0_g28_25

output:
 file "*.html"  into g28_25_outputFileHTML00
 file "*csv" optional true  into g28_25_csvFile11

"""

#!/usr/bin/env Rscript 

rmarkdown::render("${rmk}", clean=TRUE, output_format="html_document", output_dir=".")

"""
}


process Build_Consensus_parse_log_BC {

input:
 set val(name),file(log_file) from g22_10_logFile1_g22_12
 val mate from g_1_mate_g22_12

output:
 set val(name),file("*.tab")  into g22_12_logFile0_g22_14

script:
readArray = log_file.toString()

"""
ParseLog.py -l ${readArray} -f BARCODE SEQCOUNT CONSCOUNT PRCONS PRFREQ ERROR
"""

}


process Build_Consensus_report_Build_Consensus {

input:
 set val(name),file(log_files) from g22_12_logFile0_g22_14
 val matee from g_1_mate_g22_14

output:
 file "*.rmd"  into g22_14_rMarkdown0_g22_16




shell:

if(matee=="pair"){
	readArray = log_files.toString().split(' ')	
	R1 = readArray[0]
	R2 = readArray[1]

	'''
	#!/usr/bin/env perl
	
	
	my $script = <<'EOF';
	
	
	
	```{R, message=FALSE, echo=FALSE, results="hide"}
	# Setup
	library(prestor)
	library(knitr)
	library(captioner)
	
	plot_titles <- c("Read 1", "Read 2")
	if (!exists("tables")) { tables <- captioner(prefix="Table") }
	if (!exists("figures")) { figures <- captioner(prefix="Figure") }
	figures("cons_size", 
	        paste("Histogram of UMI read group sizes (reads per UMI) for",  
	              plot_titles[1], "(top) and", plot_titles[2], "(bottom).",
	              "The x-axis indicates the number of reads in a UMI group and the y-axis is the 
	               number of UMI groups with that size. The Consensus and Total bars are overlayed 
	               (not stacked) histograms indicating whether the distribution has been calculated 
	               using the total number of reads (Total) or only those reads used for consensus 
	               generation (Consensus)."))
	figures("cons_prfreq", 
	        paste("Histograms showing the distribution of majority primer frequency for all UMI read groups for",
	              plot_titles[1], "(top) and", plot_titles[2], "(bottom)."))
	figures("cons_prsize", 
	        paste("Violin plots showing the distribution of UMI read group sizes by majority primer for",
	              plot_titles[1], "(top) and", plot_titles[2], "(bottom).",
	              "Only groups with majority primer frequency over the PRFREQ threshold set when running
	               BuildConsensus. Meaning, only retained UMI groups."))
	figures("cons_error", 
	        paste("Histogram showing the distribution of UMI read group error rates for",
	              plot_titles[1], "(top) and", plot_titles[2], "(bottom)."))
	figures("cons_prerror", 
	        paste("Violin plots showing the distribution of UMI read group error rates by majority primer for",
	              plot_titles[1], "(top) and", plot_titles[2], "(bottom).",
	              "Only groups with majority primer frequency over the PRFREQ threshold set when 
	               running BuildConsensus. Meaning, only retained UMI groups."))
	```
	
	```{r, echo=FALSE}
	consensus_log_1 <- loadLogTable(file.path(".", "!{R1}"))
	consensus_log_2 <- loadLogTable(file.path(".", "!{R2}"))
	```
	
	# Generation of UMI Consensus Sequences
	
	Reads sharing the same UMI are collapsed into a single consensus sequence by
	the BuildConsensus tool. BuildConsensus considers several factors in determining
	the final consensus sequence, including the number of reads in a UMI group, 
	Phred quality scores (`Q`), primer annotations, and the number of mismatches 
	within a UMI group. Quality scores are used to resolve conflicting base calls in
	a UMI read group and the final consensus sequence is assigned consensus quality 
	scores derived from the individual base quality scores. The numbers of reads in a UMI
	group, number of matching primer annotations, and error rate (average base mismatches from 
	consensus) are used as strict cut-offs for exclusion of erroneous UMI read groups.
	Additionally, individual reads are excluded whose primer annotation differs from 
	the majority in cases where there are sufficient number of reads exceeding 
	the primer consensus cut-off.
	
	## Reads per UMI
	
	```{r, echo=FALSE, warning=FALSE}
	plotBuildConsensus(consensus_log_1, consensus_log_2, titles=plot_titles, 
	                   style="size", sizing="figure")
	```
	
	`r figures("cons_size")`
	
	## UMI read group primer frequencies
	
	```{r, echo=FALSE, warning=FALSE}
	plotBuildConsensus(consensus_log_1, consensus_log_2, titles=plot_titles, 
	                   style="prfreq", sizing="figure")
	```
	
	`r figures("cons_prfreq")`
	
	```{r, echo=FALSE, warning=FALSE}
	plotBuildConsensus(consensus_log_1, consensus_log_2, titles=plot_titles, 
	                   style="prsize", sizing="figure")
	```
	
	`r figures("cons_prsize")`
	
	## UMI read group error rates
	
	```{r, echo=FALSE, warning=FALSE}
	plotBuildConsensus(consensus_log_1, consensus_log_2, titles=plot_titles, 
	                   style="error", sizing="figure")
	```
	
	`r figures("cons_error")`
	
	```{r, echo=FALSE, warning=FALSE}
	plotBuildConsensus(consensus_log_1, consensus_log_2, titles=plot_titles, 
	                   style="prerror", sizing="figure")
	```
	
	`r figures("cons_prerror")`
	
	EOF
	
	open OUT, ">!{name}.rmd";
	print OUT $script;
	close OUT;
	
	'''

}else{
	
	readArray = log_files.toString().split(' ')
	R1 = readArray[0]
	
	'''
	#!/usr/bin/env perl
	
	
	my $script = <<'EOF';
	
	
	
		
	```{R, message=FALSE, echo=FALSE, results="hide"}
	# Setup
	library(prestor)
	library(knitr)
	library(captioner)
	
	if (!exists("tables")) { tables <- captioner(prefix="Table") }
	if (!exists("figures")) { figures <- captioner(prefix="Figure") }
	figures("cons_size", "Histogram of UMI read group sizes (reads per UMI). 
	                      The x-axis indicates the number of reads 
	                      in a UMI group and the y-axis is the number of UMI groups 
	                      with that size. The Consensus and Total bars are overlayed
	                      (not stacked) histograms indicating whether the distribution
	                      has been calculated using the total number of reads (Total)
	                      or only those reads used for consensus generation (Consensus).")
	figures("cons_error", "Histogram showing the distribution of UMI read group error rates.")
	```
	
	```{r, echo=FALSE}
	consensus_log <- loadLogTable(file.path(".", "!{R1}"))
	```
	
	# Generation of UMI Consensus Sequences
	
	Reads sharing the same UMI are collapsed into a single consensus sequence by
	the BuildConsensus tool. BuildConsensus considers several factors in determining
	the final consensus sequence, including the number of reads in a UMI group, 
	Phred quality scores (`Q`), primer annotations, and the number of mismatches 
	within a UMI group. Quality scores are used to resolve conflicting base calls in
	a UMI read group and the final consensus sequence is assigned consensus quality 
	scores derived from the individual base quality scores. The numbers of reads in a UMI
	group, number of matching primer annotations, and error rate (average base mismatches from 
	consensus) are used as strict cut-offs for exclusion of erroneous UMI read groups.
	Additionally, individual reads are excluded whose primer annotation differs from 
	the majority in cases where there are sufficient number of reads exceeding 
	the primer consensus cut-off.
	
	## Reads per UMI
	
	```{r, echo=FALSE, warning=FALSE}
	plotBuildConsensus(consensus_log, style="size", sizing="figure")
	```
	
	`r figures("cons_size")`
	
	## UMI read group error rates
	
	```{r, echo=FALSE, warning=FALSE}
	plotBuildConsensus(consensus_log, style="error", sizing="figure")
	```
	
	`r figures("cons_error")`
	
	EOF
	
	open OUT, ">!{name}.rmd";
	print OUT $script;
	close OUT;
	
	'''
}

}


process Build_Consensus_render_rmarkdown {

input:
 file rmk from g22_14_rMarkdown0_g22_16

output:
 file "*.html"  into g22_16_outputFileHTML00
 file "*csv" optional true  into g22_16_csvFile11

"""

#!/usr/bin/env Rscript 

rmarkdown::render("${rmk}", clean=TRUE, output_format="html_document", output_dir=".")

"""
}


process make_report_pipeline_cat_all_file {

input:
 set val(name), file(log_file) from g52_0_logFile3_g61_0
 set val(name), file(log_file) from g22_10_logFile2_g61_0
 set val(name), file(log_file) from g37_9_logFile1_g61_0

output:
 set val(name), file("all_out_file.log")  into g61_0_logFile0_g61_2, g61_0_logFile0_g61_10

script:
readArray = log_file.toString()

"""

echo $readArray
cat out* >> all_out_file.log
"""

}


process make_report_pipeline_report_pipeline {

input:
 set val(name), file(log_files) from g61_0_logFile0_g61_2

output:
 file "*.rmd"  into g61_2_rMarkdown0_g61_10


shell:

readArray = log_files.toString().split(' ')
R1 = readArray[0]

'''
#!/usr/bin/env perl


my $script = <<'EOF';


```{r, message=FALSE, echo=FALSE, results="hide"}
# Setup
library(prestor)
library(knitr)
library(captioner)

plot_titles <- c("Read 1", "Read 2")
if (!exists("tables")) { tables <- captioner(prefix="Table") }
if (!exists("figures")) { figures <- captioner(prefix="Figure") }
tables("count", 
       "The count of reads that passed and failed each processing step.")
figures("steps", 
        paste("The number of reads or read sets retained at each processing step. 
               Shown as raw counts (top) and percentages of input from the previous 
               step (bottom). Steps having more than one column display individual values for", 
              plot_titles[1], "(first column) and", plot_titles[2], "(second column)."))
```

```{r, echo=FALSE}
console_log <- loadConsoleLog(file.path(".","!{R1}"))
```

# Summary of Processing Steps

```{r, echo=FALSE}
count_df <- plotConsoleLog(console_log, sizing="figure")

df<-count_df[,c("task", "pass", "fail")]

write.csv(df,"pipeline_statistics.csv") 
```

`r figures("steps")`

```{r, echo=FALSE}
kable(count_df[c("step", "task", "total", "pass", "fail")],
      col.names=c("Step", "Task", "Input", "Passed", "Failed"),
      digits=3)
```

`r tables("count")`


EOF
	
open OUT, ">pipeline_statistic_!{name}.rmd";
print OUT $script;
close OUT;

'''
}


process make_report_pipeline_presto_render_rmarkdown {

input:
 file rmk from g61_2_rMarkdown0_g61_10
 file log_file from g61_0_logFile0_g61_10

output:
 file "*.html"  into g61_10_outputFileHTML00
 file "*csv" optional true  into g61_10_csvFile11

"""

#!/usr/bin/env Rscript 

rmarkdown::render("${rmk}", clean=TRUE, output_format="html_document", output_dir=".")

"""
}


process Mask_Primer_parse_log_MP {

input:
 val mate from g_1_mate_g38_9
 set val(name), file(log_file) from g38_11_logFile2_g38_9

output:
 file "*table.tab"  into g38_9_logFile0_g38_12, g38_9_logFile0_g38_19

script:
readArray = log_file.toString()	

"""
ParseLog.py -l ${readArray}  -f ID PRIMER BARCODE ERROR
"""

}


process Mask_Primer_try_report_maskprimer {

input:
 file primers from g38_9_logFile0_g38_12
 val mate from g_1_mate_g38_12

output:
 file "*.rmd"  into g38_12_rMarkdown0_g38_19


shell:

if(mate=="pair"){
	readArray = primers.toString().split(' ')	
	primers_1 = readArray[0]
	primers_2 = readArray[1]
	name = primers_1 - "_table.tab"
	'''
	#!/usr/bin/env perl
	
	
	my $script = <<'EOF';
	
	
	```{r, message=FALSE, echo=FALSE, results="hide"}
	
	# Setup
	library(prestor)
	library(knitr)
	library(captioner)
	
	
	plot_titles<- c("Read 1", "Read 2")
	print(plot_titles)
	if (!exists("tables")) { tables <- captioner(prefix="Table") }
	if (!exists("figures")) { figures <- captioner(prefix="Figure") }
	figures("primers_count", 
	        paste("Count of assigned primers for",  plot_titles[1], "(top) and", plot_titles[2], "(bottom).",
	              "The bar height indicates the total reads assigned to the given primer,
	               stacked for those under the error rate threshold (Pass) and
	               over the threshold (Fail)."))
	figures("primers_hist", 
	        paste("Distribution of primer match error rates for", plot_titles[1], "(top) and", plot_titles[2], "(bottom).",
	              "The error rate is the percentage of mismatches between the primer sequence and the 
	               read for the best matching primer. The dotted line indicates the error threshold used."))
	figures("primers_error", 
	        paste("Distribution of primer match error rates for", plot_titles[1], "(top) and", plot_titles[2], "(bottom),",
	              "broken down by assigned primer. The error rate is the percentage of mismatches between the 
	               primer sequence and the read for the best matching primer. The dotted line indicates the error
	               threshold used."))
	```
	
	```{r, echo=FALSE}
	primer_log_1 <- loadLogTable(file.path(".", "!{primers_1}"))
	primer_log_2 <- loadLogTable(file.path(".", "!{primers_2}"))
	
	primer_log1_error <- any(is.na(primer_log_1[['ERROR']]))
	primer_log2_error<- any(is.na(primer_log_2[['ERROR']]))
	
	```
	
	# Primer Identification
	
	The MaskPrimers tool supports identification of multiplexed primers and UMIs.
	Identified primer regions may be masked (with Ns) or cut to mitigate downstream
	SHM analysis artifacts due to errors in the primer region. An annotion is added to 
	each sequences that indicates the UMI and best matching primer. In the case of
	the constant region primer, the primer annotation may also be used for isotype 
	assignment.
	
	## Count of primer matches
	
	```{r, echo=FALSE, warning=FALSE}
	if(!primer_log1_error && !primer_log2_error)
		plotMaskPrimers(primer_log_1, primer_log_2, titles=plot_titles,
	                style="count", sizing="figure")
	```
	
	`r figures("primers_count")`
	
	## Primer match error rates
	
	```{r, echo=FALSE, warning=FALSE}
	if(!primer_log1_error && !primer_log2_error)
		plotMaskPrimers(primer_log_1, primer_log_2, titles=plot_titles, 
	                style="hist", sizing="figure")
	```
	
	`r figures("primers_hist")`
	
	```{r, echo=FALSE, warning=FALSE}
	# check the error column exists 
	if(!primer_log1_error && !primer_log2_error)
		plotMaskPrimers(primer_log_1, primer_log_2, titles=plot_titles, 
	                style="error", sizing="figure")
	```
	
	`r figures("primers_error")`
	
	EOF
	
	open OUT, ">!{name}.rmd";
	print OUT $script;
	close OUT;
	
	'''

}else{

	readArray = primers.toString().split(' ')
	primers = readArray[0]
	name = primers - "_table.tab"
	'''
	#!/usr/bin/env perl
	
	
	my $script = <<'EOF';
	
	
	```{r, message=FALSE, echo=FALSE, results="hide"}
	
	# Setup
	library(prestor)
	library(knitr)
	library(captioner)
	
	
	plot_titles<- c("Read")
	print(plot_titles)
	if (!exists("tables")) { tables <- captioner(prefix="Table") }
	if (!exists("figures")) { figures <- captioner(prefix="Figure") }
	figures("primers_count", 
	        paste("Count of assigned primers for",  plot_titles[1],
	              "The bar height indicates the total reads assigned to the given primer,
	               stacked for those under the error rate threshold (Pass) and
	               over the threshold (Fail)."))
	figures("primers_hist", 
	        paste("Distribution of primer match error rates for", plot_titles[1],
	              "The error rate is the percentage of mismatches between the primer sequence and the 
	               read for the best matching primer. The dotted line indicates the error threshold used."))
	figures("primers_error", 
	        paste("Distribution of primer match error rates for", plot_titles[1],
	              "broken down by assigned primer. The error rate is the percentage of mismatches between the 
	               primer sequence and the read for the best matching primer. The dotted line indicates the error
	               threshold used."))
	```
	
	```{r, echo=FALSE}
	primer_log_1 <- loadLogTable(file.path(".", "!{primers}"))
	```
	
	# Primer Identification
	
	The MaskPrimers tool supports identification of multiplexed primers and UMIs.
	Identified primer regions may be masked (with Ns) or cut to mitigate downstream
	SHM analysis artifacts due to errors in the primer region. An annotion is added to 
	each sequences that indicates the UMI and best matching primer. In the case of
	the constant region primer, the primer annotation may also be used for isotype 
	assignment.
	
	## Count of primer matches
	
	```{r, echo=FALSE, warning=FALSE}
	plotMaskPrimers(primer_log_1, titles=plot_titles,
	                style="count", sizing="figure")
	```
	
	`r figures("primers_count")`
	
	## Primer match error rates
	
	```{r, echo=FALSE, warning=FALSE}
	plotMaskPrimers(primer_log_1, titles=plot_titles, 
	                style="hist", sizing="figure")
	```
	
	`r figures("primers_hist")`
	
	```{r, echo=FALSE, warning=FALSE}
	plotMaskPrimers(primer_log_1, titles=plot_titles, 
	                style="error", sizing="figure")
	```
	
	`r figures("primers_error")`
	
	EOF
	
	open OUT, ">!{name}.rmd";
	print OUT $script;
	close OUT;
	
	'''
}
}


process Mask_Primer_presto_render_rmarkdown {

input:
 file rmk from g38_12_rMarkdown0_g38_19
 file log_file from g38_9_logFile0_g38_19

output:
 file "*.html"  into g38_19_outputFileHTML00
 file "*csv" optional true  into g38_19_csvFile11

"""

#!/usr/bin/env Rscript 

rmarkdown::render("${rmk}", clean=TRUE, output_format="html_document", output_dir=".")

"""
}


process metadata {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.json$/) "metadata/$filename"}

output:
 file "*.json"  into g_69_jsonFile00

script:
metadata = params.metadata.metadata
"""
#!/usr/bin/env Rscript

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  install.packages("jsonlite")
}
library(jsonlite)

data <- read_json("${metadata}") 

versions <- lapply(1:length(data), function(i){
	
	docker <- data[i]
	tool <- names(data)[i]
	
	if(grepl("Custom", docker)){
		ver <- "0.0"
	}else{
		ver <- system(paste0(tool," --version"), intern = TRUE)
		ver <- gsub(paste0(tool,": "), "", ver)
	}
	ver
	
})

names(versions) <- names(data)

json_data <- list(
  sample = list(
    data_processing = list(
      preprocessing = list(
        software_versions = versions
	   )
	 )
  )
)

# Convert to JSON string without enclosing scalar values in arrays
json_string <- toJSON(json_data, pretty = TRUE, auto_unbox = TRUE)
print(json_string)
# Write the JSON string to a file
writeLines(json_string, "pre_processed_metadata.json")
"""

}


workflow.onComplete {
println "##Pipeline execution summary##"
println "---------------------------"
println "##Completed at: $workflow.complete"
println "##Duration: ${workflow.duration}"
println "##Success: ${workflow.success ? 'OK' : 'failed' }"
println "##Exit status: ${workflow.exitStatus}"
}
