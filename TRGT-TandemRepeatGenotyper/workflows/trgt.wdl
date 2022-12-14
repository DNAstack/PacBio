version 1.0

workflow trgt {
	input {
		File ref 
		File repeats 
		File aligned_bam
		File aligned_bai
		String repeat_id
		String container_registry
	}

	String aligned_bam_basename = basename(aligned_bam, ".bam")

	call genotype_repeats {
		input:
			ref = ref,
			repeats = repeats,
			aligned_bam = aligned_bam,
			aligned_bai = aligned_bai,
			aligned_bam_basename = aligned_bam_basename,
			container_registry = container_registry
	}

	call sort_index_vcf {
		input:
			trgt_vcf = genotype_repeats.trgt_vcf,
			aligned_bam_basename = aligned_bam_basename,
			container_registry = container_registry
	}

	call sort_index_spanning_bam {
		input:
			trgt_bam = genotype_repeats.trgt_bam,
			aligned_bam_basename = aligned_bam_basename,
			container_registry = container_registry
	}

	call visualize_repeats {
		input:
			ref = ref,
			repeats = repeats,
			sorted_trgt_vcf = sort_index_vcf.sorted_trgt_vcf,
			sorted_trgt_vcf_index = sort_index_vcf.sorted_trgt_vcf_index,
			sorted_trgt_bam = sort_index_spanning_bam.sorted_trgt_bam,
			sorted_trgt_bam_index = sort_index_spanning_bam.sorted_trgt_bam_index,
			repeat_id = repeat_id,
			container_registry = container_registry
	}

	output {
		File sorted_trgt_vcf = sort_index_vcf.sorted_trgt_vcf
		File sorted_trgt_vcf_index = sort_index_vcf.sorted_trgt_vcf_index

		File sorted_trgt_bam = sort_index_spanning_bam.sorted_trgt_bam
		File sorted_trgt_bam_index = sort_index_spanning_bam.sorted_trgt_bam_index

		File pileup_image = visualize_repeats.pileup_image
	}

	meta {
		author: "PacBio"
		email: "bioinformatics@dnastack.com"
		description: "# PacBio TRGT: Tandem Repeat Genotyper implemented in Workflow Description Language (WDL)\n\n![PacBio logo](https://raw.githubusercontent.com/DNAstack/PacBio/main/pacbio-logo-small.png)\n\nThis repository contains a [WDL workflow](https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md) for processing PacBio HiFi data using the [TRGT tool](https://github.com/pacificBiosciences/trgt/). `trgt` profiles sequence composition, mosaicism, and CpG methylation of analyzed repeats.\n\nDocker images containing the tools used by this workflow can be explored [in DNAstack's image repository](https://github.com/dnastack/bioinformatics-public-docker-images), or [on Dockerhub](https://hub.docker.com/u/dnastack).\n\n\n## Workflow inputs\n\nAn input template file with some defaults predefined can be found [here](./workflows/inputs.json).\nSome example input files can be found [in PacBio's `trgt` repository](https://github.com/PacificBiosciences/trgt/tree/main/example).\n\n| Input | Description |\n| :- | :- |\n| `ref` | The reference genome that was used for read alignment (FASTA) |\n| `aligned_bam`, `aligned_bai` | Aligned HiFi reads (BAM) and index (BAI) |\n| `repeats` | The repeat definition file with reference coordinates and structure of tandem repeats (BED) |\n| `repeat_id` | ID of the repeat to visualize |\n| `container_registry` | Registry that hosts workflow containers. All containers are hosted in [DNAstack's Dockerhub](https://hub.docker.com/u/dnastack) [`dnastack`] |\n\n\n## Workflow outputs\n\n| Output | Description |\n| :- | :- |\n| `sorted_trgt_vcf`, `sorted_trgt_vcf_index` | Sorted VCF file and index that contains repeat genotypes; output by `trgt` |\n| `sorted_trgt_bam`, `sorted_trgt_bam_index` | Sorted BAM file and index that contains pieces of HiFi reads that fully span the repeat sequences; output by `trgt` |\n| `pileup_image` | An SVG file that contains the pileup read image; output by `trvz` |\n\n\n## Running workflows\n\n### Required software\n\n- [Docker](https://docs.docker.com/get-docker/)\n- [Cromwell](https://github.com/broadinstitute/cromwell/releases) & Java (8+) OR [miniwdl](https://github.com/chanzuckerberg/miniwdl/releases) & python3\n\n### Running using Cromwell\n\nFrom the root of the repository, run:\n\n```bash\njava -jar /path/to/cromwell.jar run /path/to/workflow.wdl -i /path/to/inputs.json\n```\n\nOutput and execution files will be located in the `cromwell-executions` directory. When the workflow finishes successfully, it will output JSON (to stdout) specifying the full path to each output file.\n\n\n### Running using miniwdl\n\nThis command assumes you have `miniwdl` available on your command line. If `miniwdl` is not available, try installing using `pip install miniwdl`.\n\n```bash\nminiwdl run /path/to/workflow.wdl -i /path/to/inputs.json\n```\n\nOutput and execution files will be located in a dated directory (e.g. named `20200704_073415_main`). When the workflow finishes successfully, it will output JSON (to stdout) specifying the full path to each output file.\n\n\n## Future work\n\n* Add an optional alignment step so that FASTQ can be an input\n* Improve workflow by changing `repeat_id` to grep the ID in the `repeat.bed` file instead of feeding it a literal single string in order to loop over TRVZ to generate multiple pile-up images\n"
	}
}

task genotype_repeats {
	input {
		File ref
		File repeats
		File aligned_bam
		File aligned_bai

		String aligned_bam_basename
		String container_registry
	}

	Int disk_size = ceil((size(ref, "GB") + size(repeats, "GB") + size(aligned_bam, "GB")) * 2 + 20)

	command <<<
		trgt --genome ~{ref} \
			--repeats ~{repeats} \
			--reads ~{aligned_bam} \
			--output-prefix ~{aligned_bam_basename}
	>>>

	output {
		File trgt_vcf = "~{aligned_bam_basename}.vcf.gz"
		File trgt_bam = "~{aligned_bam_basename}.spanning.bam"
	}

	runtime {
		docker: "~{container_registry}/pacbio_trgt_tools:0.0.1"
		cpu: 1
		memory: "7.5 GB"
		disks: "local-disk " + disk_size + " HDD"
		preemptible: 2
	}
}

task sort_index_vcf {
	input {
		File trgt_vcf

		String aligned_bam_basename
		String container_registry
	}

	Int disk_size = ceil(size(trgt_vcf, "GB") * 2 + 20)


	command <<<
		bcftools sort \
			-Ob -o "~{aligned_bam_basename}.sorted.vcf.gz" \
			~{trgt_vcf}

		bcftools index \
			"~{aligned_bam_basename}.sorted.vcf.gz" \
			-o "~{aligned_bam_basename}.sorted.vcf.gz.tbi"
	>>>

	output {
		File sorted_trgt_vcf = "~{aligned_bam_basename}.sorted.vcf.gz"
		File sorted_trgt_vcf_index = "~{aligned_bam_basename}.sorted.vcf.gz.tbi"
	}

	runtime {
		docker: "~{container_registry}/pacbio_trgt_tools:0.0.1"
		cpu: 1
		memory: "7.5 GB"
		disks: "local-disk " + disk_size + " HDD"
		preemptible: 2
	}
}


task sort_index_spanning_bam {
	input {
		File trgt_bam

		String aligned_bam_basename
		String container_registry
	}

	Int disk_size = ceil(size(trgt_bam, "GB") * 2 + 20)

	command <<<
		samtools sort \
			-o "~{aligned_bam_basename}.spanning.sorted.bam" \
			~{trgt_bam}

		samtools index \
			"~{aligned_bam_basename}.spanning.sorted.bam" \
			-o "~{aligned_bam_basename}.spanning.sorted.bam.bai"
	>>>

	output {
		File sorted_trgt_bam = "~{aligned_bam_basename}.spanning.sorted.bam"
		File sorted_trgt_bam_index = "~{aligned_bam_basename}.spanning.sorted.bam.bai"
	}

	runtime {
		docker: "~{container_registry}/pacbio_trgt_tools:0.0.1"
		cpu: 1
		memory: "7.5 GB"
		disks: "local-disk " + disk_size + " HDD"
		preemptible: 2
	}
}

task visualize_repeats {
	input {
		File ref
		File repeats
		File sorted_trgt_vcf
		File sorted_trgt_vcf_index
		File sorted_trgt_bam
		File sorted_trgt_bam_index

		String repeat_id
		String container_registry
	}

	Int disk_size = ceil((size(ref, "GB") + size(repeats, "GB") + size(sorted_trgt_vcf, "GB") + size(sorted_trgt_bam, "GB")) * 2 + 20)

	command <<<
		trvz --genome ~{ref} \
			--repeats ~{repeats} \
			--vcf ~{sorted_trgt_vcf} \
			--spanning-reads ~{sorted_trgt_bam} \
			--repeat-id ~{repeat_id} \
			--image "~{repeat_id}.svg"
	>>>

	output {
		File pileup_image = "~{repeat_id}.svg"
	}

	runtime {
		docker: "~{container_registry}/pacbio_trgt_tools:0.0.1"
		cpu: 1
		memory: "7.5 GB"
		disks: "local-disk " + disk_size + " HDD"
		preemptible: 2
	}
}
