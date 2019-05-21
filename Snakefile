
samples = ['Run29_D8','Run32_Day20','Run32_Day8_TIM3p','Run32_Day8_TIM3n']



# rule DGE_microarray:
# 	input: 
# 	output:
#  	script:
#  		'code/Microarray_DGEA.R'

rule DGE_RNAseq:
	input:
		expand('data/{sample}_Count.tsv',sample=samples),
		expand('data/{sample}_Design.tsv',sample=samples)
	output:
		expand('output/DEG/Filter_DEG_{sample}_WT_KO.xlsx',sample=samples),
		expand('output/DEG/Rlog_{sample}.csv',sample=samples)
	conda: 'envs/DGE_RNA.yaml'
	script:
		'code/DEG.R'


rule generate_heatmaps:
	input:
		expand('output/DEG/Filter_DEG_{sample}_WT_KO.xlsx',sample=samples),
		expand('output/DEG/Rlog_{sample}.csv',sample=samples)
	output:
		expand('output/Heatmaps/Rlog_Scale_{sample}.pdf',sample=samples)
	script:
		'code/Heatmap.R'