configfile: "config.yaml"

input_dir=config["input_dir"]
out_dir=config["out_dir"]
refGenome_path=config["refGenome_path"]
fq_id_txt=config["fq_id_txt"]
trim_adaptor_path=config["trim_adaptor_path"]

log_dir=out_dir + "/log"
result_dir=out_dir + "/result"

md5_log_dir=log_dir + "/md5"

trim_res_dir=result_dir + "/trim"
trim_log_dir=log_dir + "/trim"

bwa_res_dir=result_dir + "/bwa"
bwa_log_dir=log_dir + "/bwa"

gatk_md_dir=result_dir + "/gatk_md"
gatk_md_log_dir=log_dir + "/gatk_md"

gatk_gvcf_dir=result_dir + "/gatk_gvcf"
gatk_gvcf_log_dir=log_dir + "/gatk_gvcf"

import os,glob
fq1_list=[os.path.basename(x) for x in glob.glob(input_dir + "/*_1.fq.gz")]
fq1_list_clean=[x.replace("_1.fq.gz","") for x in fq1_list]
fq2_list=[os.path.basename(x) for x in glob.glob(input_dir + "/*_2.fq.gz")]
fq2_list_clean=[x.replace("_2.fq.gz","") for x in fq2_list]
fq_avil=set(fq1_list_clean).intersection(set(fq2_list_clean))
print("%s paired fq samples detected:" %(len(fq_avil)))

fq_quest_raw=[]
with open(fq_id_txt,"r") as f:
	for line in f.readlines():
		fq_quest_raw.append(line.strip("\n"))
fq_quest=set(fq_quest_raw)
f.close()
print("%s paired fq samples quested:" %(len(fq_quest)))

fq_inter=set(fq_avil).intersection(fq_quest)
print("%s paired fq samples quested and available:" %(len(fq_inter)))

print(fq_inter)
fq=fq_inter

# def readTwoColTxtAsDic(TwoColTxt):
# 	with open(TwoColTxt,"r",encoding="utf-8") as f:
# 		my_dict={}
# 		for line in f.readlines():
# 			line=line.strip('\n')
# 			b=line.split('  ',1)
# 			fq=b[1].rsplit('/',1)[-1]
# 			my_dict[fq]=b[0]
# 	return my_dict
		


rule all:
		input:
				expand(input_dir + "/{fq}_1.fq.gz",input_dir=input_dir,fq=fq),
				expand(input_dir + "/{fq}_2.fq.gz",input_dir=input_dir,fq=fq),
				expand(md5_log_dir + "/{fq}.log",fq=fq),
				expand(trim_log_dir + "/{fq}.log",fq=fq),
				expand(bwa_log_dir + "/{fq}.log",fq=fq,allow_missing=True),
				expand(gatk_md_log_dir +"/{fq}.log",fq=fq,allow_missing=True),
				expand(gatk_gvcf_log_dir +"/{fq}.log",fq=fq,allow_missing=True),
				log_dir + "/join_check_md5sum_info.txt",
				log_dir + "/join_trim_info.txt",
				log_dir + "/join_bwa2_info.txt",
				log_dir + "/join_gatk_markDuplicates_info.txt",
				log_dir + "/join_gatk_HaplotypeCaller_info.txt"


rule check_md5sum:
		input:
				fq1= input_dir + "/{fq}_1.fq.gz",
				fq2= input_dir + "/{fq}_2.fq.gz"
		output:
				md5_re= input_dir + "/{fq}.md5_re1.txt",
				md5_log= md5_log_dir + "/{fq}.log"
		params:
				md5_golden = input_dir + "/{fq}.md5.gold.txt"
		priority: 
				1
		threads:
				20
		# resources:
		# 		readDic=readTwoColTxtAsDic
		script:
				"check_md5sum.py"

rule join_check_md5sum_info:
		input:
				expand(md5_log_dir + "/{fq}.log", fq=fq)
		output:
				log_dir + "/join_check_md5sum_info.txt"
		priority: 
				1
		shell:
				"""
				echo $(cat {input} |grep "PASS:" | wc -l | awk -F ' ' '{{print$1}}') " / " $(cat {input} | wc -l | awk -F ' ' '{{print$1}}') " succeed" > {output}
				cat {input} >> {output}
				"""

rule trim:
		input:
				md5_log= md5_log_dir + "/{fq}.log",
				fq1= input_dir + "/{fq}_1.fq.gz",
				fq2= input_dir + "/{fq}_2.fq.gz"
		output:
				trim_paired_fq1= trim_res_dir + "/{fq}_1.paired.fq.gz",
				trim_paired_fq2= trim_res_dir + "/{fq}_2.paired.fq.gz",
				trim_unpaired_fq1= trim_res_dir + "/{fq}_1.unpaired.fq.gz",
				trim_unpaired_fq2= trim_res_dir + "/{fq}_2.unpaired.fq.gz",
				trim_summary= trim_res_dir + "/{fq}_trim.summary",
				trim_log= trim_log_dir + "/{fq}.log"
		params:
				trim_adaptor=trim_adaptor_path
		priority:
				2
		threads:
				10
		shell:
				"""
				md5_status=$(cat {input.md5_log} | grep 'PASS: ' | wc -l)
				( [ $md5_status == 1 ] && \
				trimmomatic PE \
				-threads {threads} \
				-summary {output.trim_summary} \
				{input.fq1} {input.fq2} \
				{output.trim_paired_fq1} {output.trim_unpaired_fq1} \
				{output.trim_paired_fq2} {output.trim_unpaired_fq2} \
				ILLUMINACLIP:{params.trim_adaptor}:2:30:10:8:true \
				LEADING:3 \
				TRAILING:3 \
				SLIDINGWINDOW:4:15 \
				MINLEN:36) && \
				(echo "PASS: " {wildcards.fq} " trim succeed" > {output.trim_log}) || \
				(touch {output.trim_summary}; \
				touch {output.trim_paired_fq1}; \
				touch {output.trim_paired_fq2}; \
				touch {output.trim_unpaired_fq1}; \
				touch {output.trim_unpaired_fq2}; \
				echo "FAIL: " {wildcards.fq} " trim failed" > {output.trim_log})
				"""

rule join_trim_info:
		input:
				expand(trim_log_dir + "/{fq}.log",fq=fq)
		output:
				log_dir + "/join_trim_info.txt"
		priority: 
				1
		shell:
				"""
				echo $(cat {input} |grep "PASS: " | wc -l | awk -F ' ' '{{print$1}}') " / " $(cat {input} | wc -l | awk -F ' ' '{{print$1}}') " succeed" > {output}
				cat {input} >> {output}
				"""

rule bwa2:
		input:
				fq1= trim_res_dir + "/{fq}_1.paired.fq.gz",
				fq2= trim_res_dir + "/{fq}_2.paired.fq.gz",
				trim_log= trim_log_dir + "/{fq}.log",
				ref=refGenome_path
		output:
				bam = bwa_res_dir + "/{fq}.bam",
				bwa_log = bwa_log_dir + "/{fq}.log"
		params:
				ref_index = refGenome_path + ".bwt.2bit.64"
		priority:
				3
		threads:
				10
		shell:
				"""
				if [ ! -f {params.ref_index} ] ; then
					bwa-mem2 index {input.ref};
				fi; 

				trim_status=$(cat {input.trim_log} | grep 'PASS: ' | wc -l);
				( [ $trim_status == 1 ] && \
				bwa-mem2 mem -t {threads} -k 32 -M \
				-R '@RG\tID:{wildcards.fq}\tLB:{wildcards.fq}\tSM:{wildcards.fq}\tPL:Illumina\tPU:{wildcards.fq}' \
				{input.ref} \
				{input.fq1} {input.fq2} | \
				samtools view -@ 10 -h -q 30 -F 4 -F 256 | \
				grep -v XA:Z | grep -v SA:Z | \
				samtools sort -@ 10 \
				-o {output.bam} ) && \
				(echo "PASS: "{wildcards.fq} " bwa_align succeed" > {output.bwa_log}) || \
				(touch {output.bam}; \
				echo "FAIL: " {wildcards.fq} " bwa_align failed" > {output.bwa_log})

				"""

rule join_bwa2_info:
		input:
				expand(bwa_log_dir + "/{fq}.log",fq=fq)
		output:
				log_dir + "/join_bwa2_info.txt"
		priority: 
				1
		shell:
				"""
				echo $(cat {input} |grep "PASS: " | wc -l | awk -F ' ' '{{print$1}}') " / " $(cat {input} | wc -l | awk -F ' ' '{{print$1}}') > {output}
				cat {input} >> {output}
				"""

rule gatk_markDuplicates:
		input:
				bwa_log = bwa_log_dir + "/{fq}.log",
				bam = bwa_res_dir + "/{fq}.bam",
				ref=refGenome_path
		output:
				md_bam = gatk_md_dir + "/{fq}.md.bam",
				md_metrices = gatk_md_dir + "/{fq}.markdup_metrics.txt",
				gatk_md_log = gatk_md_log_dir + "/{fq}.log"
		threads:
				6
		priority: 
				4
		shell:
				"""
				bwa_status=$(cat {input.bwa_log} | grep 'PASS: ' | wc -l);
				gatkIndexFile=$(echo {input.ref} | awk -F'.fa$' '{{print$1".dict"}}');
				if [ ! -f $gatkIndexFile ] ; then 
					echo "GATK indexing reference genome...";
					gatk CreateSequenceDictionary -R {input.ref} -O $gatkIndexFile
				fi;

				( [ $bwa_status == 1 ] ) && \
				(gatk  MarkDuplicates \
				-I {input.bam} \
				-O {output.md_bam} \
				-M {output.md_metrices} ) && \
				(samtools index -@ 10 {output.md_bam}) && \
				(echo "PASS: " {wildcards.fq} " star_pass2_align succeed" > {output.gatk_md_log}) || \
				(touch {output.md_bam} ; echo "FAIL: " {wildcards.fq} " gatk_markDuplicates failed" > {output.gatk_md_log})
				"""

rule join_gatk_markDuplicates_info:
		input:
				expand(gatk_md_log_dir + "/{fq}.log",fq=fq)
		output:
				log_dir + "/join_gatk_markDuplicates_info.txt"
		priority: 
				1
		shell:
				"""
				echo $(cat {input} |grep "PASS: " | wc -l | awk -F ' ' '{{print$1}}') " / " $(cat {input} | wc -l | awk -F ' ' '{{print$1}}') > {output}
				cat {input} >> {output}
				"""



rule gatk_HaplotypeCaller:
		input:
				gatk_md_log = gatk_md_log_dir + "/{fq}.log",
				md_bam = gatk_md_dir + "/{fq}.md.bam",
				ref = refGenome_path
		output:
				gvcf = gatk_gvcf_dir + "/{fq}.g.vcf.gz",
				gvcf_log = gatk_gvcf_log_dir + "/{fq}.log"
		threads:
				6
		priority: 
				5
		shell:
				"""
				gatk_md_status=$(cat {input.gatk_md_log} | grep 'PASS: ' | wc -l);
				( [ $gatk_md_status == 1 ] )&& \
				(gatk HaplotypeCaller  \
				-R {input.ref} \
				-I {input.md_bam} \
				-O {output.gvcf} \
				-ERC GVCF -ploidy 2 -stand-call-conf 30.0 \
				--native-pair-hmm-threads 10) && \
				(echo "PASS: " {wildcards.fq} " rule gatk_HaplotypeCaller succeed" > {output.gvcf_log}) || \
				(touch {output.gvcf} ; echo "FAIL: " {wildcards.fq} " gatk_HaplotypeCaller failed" > {output.gvcf_log})
				"""

rule join_gatk_HaplotypeCaller_info:
		input:
				expand(gatk_gvcf_log_dir + "/{fq}.log",fq=fq)
		output:
				log_dir + "/join_gatk_HaplotypeCaller_info.txt"
		priority: 
				1
		shell:
				"""
				echo $(cat {input} |grep "PASS: " | wc -l | awk -F ' ' '{{print$1}}') " / " $(cat {input} | wc -l | awk -F ' ' '{{print$1}}') > {output}
				cat {input} >> {output}
				"""








