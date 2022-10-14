version 1.0

workflow run_cosa {
	input {
		String accession
		File primer_trimmed_fastq

		File reference
		File reference_index

		String container_registry
	}

	call align {
		input:
			samplename = accession,
			primer_trimmed_fastq = primer_trimmed_fastq,
			reference = reference,
			container_registry = container_registry
	}

	call deepvariant {
		input:
			samplename = accession,
			aligned_bam = align.aligned_bam,
			aligned_bam_index = align.aligned_bam_index,
			reference = reference,
			reference_index = reference_index,
			container_registry = container_registry
	}

	call VCF_consensus_deepvariant {
		input:
			samplename = accession,
			aligned_bam = align.aligned_bam,
			aligned_bam_index = align.aligned_bam_index,
			vcf = deepvariant.vcf,
			reference = reference,
			container_registry = container_registry
	}

	output {
		File aligned_bam = align.aligned_bam
		File aligned_bam_index = align.aligned_bam_index
		File high_quality_variants = VCF_consensus_deepvariant.high_quality_variants
		File high_quality_variants_index = VCF_consensus_deepvariant.high_quality_variants_index
		File consensus_fa = VCF_consensus_deepvariant.consensus_fa
	}

	meta {
		author: "Heather Ward"
		email: "heather@dnastack.com"
		description: <<<
		# Analysis of PacBio SARS-CoV-2 data using the CoSA pipeline

		This repository provides a [WDL wrapper](https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md) for running ["PacBio's CoSA pipeline"](https://github.com/PacificBiosciences/CoSA) to process Pacific Biosciences SARS-CoV-2 long read HiFi data.

		The CoSA pipeline produces an assembly and a variants file. This pipeline uses DeepVariant for variant calling, but [other variant callers are available](https://github.com/PacificBiosciences/CoSA/wiki/Variant-calling-using-PacBio-HiFi-CCS-data#4-variant-calling).

		This \workflow starts from primer-trimmed fastq files; for instructions on how to convert a PacBio `subreads.bam` file to a primer-trimmed fastq or bam, see steps 1-3 [here](https://github.com/PacificBiosciences/CoSA/wiki/Variant-calling-using-PacBio-HiFi-CCS-data#1-generate-ccs-data).


		## Workflow inputs

		An \input template file with some defaults pre-defined can be found [here](https://github.com/DNAstack/PacBio/blob/main/CoSA-CoronavirusSequencingAnalysis/workflows/inputs.json).

		| Input | Description |
		|:-|:-|
		| `accession` | Sample ID |
		| `primer_trimmed_fastq` | A demultiplexed, single-sample fastq file containing primer-trimmed reads |
		| `reference`, `reference index` | [The SARS-CoV-2 reference genome](https://www.ncbi.nlm.nih.gov/nuccore/MN908947.3) and its index |
		| `container_registry` | Registry that hosts \workflow containers. All containers are hosted in ["DNAstack's Dockerhub"](https://hub.docker.com/u/dnastack) [`dnastack`] |


		## Workflow outputs

		| Output | Description |
		|:-|:-|
		| `aligned_bam`, `aligned_bam_index` | Reads aligned to the SARS-CoV-2 reference genome |
		| `high_quality_variants`, `high_quality_variants_index` | Variant calls and index |
		| `consensus_fa` | Genome assembly |


		## Test data

		PacBio SARS-CoV-2 runs are available from [the SRA](https://www.ncbi.nlm.nih.gov/sra?term=(((txid2697049%5BOrganism%3Anoexp%5D%20NOT%200%5BMbases))%20AND%20((txid2697049%5BOrganism%3Anoexp%5D%20NOT%200%5BMbases)%20AND%20%22platform%20pacbio%20smrt%22%5BProperties%5D))%20AND%20%22pacbio%20smrt%22%5BPlatform%5D). Most of this data is available as primer-trimmed fastqs, which can be used as \input directly.


		## Containers

		Docker image definitions can be found in our [bioinformatics-public-docker-images](https://github.com/DNAstack/bioinformatics-public-docker-images) repo.

		All containers are publicly hosted in ["DNAstack's container registry"](https://hub.docker.com/u/dnastack).
		>>>
	}
}

task align {
	input {
		String samplename
		File primer_trimmed_fastq
		File reference

		String container_registry
	}

	command {
		minimap2 \
			-a ~{reference} \
			~{primer_trimmed_fastq} \
			> ~{samplename}.aligned.sam

		samtools view \
			-bS ~{samplename}.aligned.sam \
		| samtools sort \
			> ~{samplename}.aligned.bam

		samtools index ~{samplename}.aligned.bam
	}

	output {
		File aligned_bam = "~{samplename}.aligned.bam"
		File aligned_bam_index = "~{samplename}.aligned.bam.bai"
	}

	runtime {
		docker: "~{container_registry}/cosa:0.0.1"
		cpu: 2
		memory: "7.5 GB"
		disks: "local-disk 100 HDD"
	}
}

task deepvariant {
	input {
		String samplename
		File aligned_bam
		File aligned_bam_index
		File reference
		File reference_index

		String container_registry
	}

	Int threads = 16

	command {
		run_deepvariant \
			--model_type PACBIO \
			--ref ~{reference} \
			--reads ~{aligned_bam} \
			--output_vcf ~{samplename}.vcf \
			--num_shards ~{threads}
	}

	output {
		File vcf = "~{samplename}.vcf"
	}

	runtime {
		docker: "google/deepvariant:1.1.0"
		cpu: threads
		memory: "32 GB"
		disks: "local-disk 100 HDD"
	}
}

task VCF_consensus_deepvariant {
	input {
		String samplename
		File aligned_bam
		File aligned_bam_index
		File vcf
		File reference

		String container_registry
	}

	command {
		samtools mpileup \
			--min-BQ 1 \
			-f ~{reference} \
			-s ~{aligned_bam} \
			> ~{samplename}.aligned.bam.mpileup

		samtools depth \
			-q 0 \
			-Q 0 \
			~{aligned_bam} \
			> ~{samplename}.aligned.bam.depth

		VCFCons.py \
			~{reference} \
			~{samplename} \
			-c 4 \
			-f 0.5 \
			--vcf_type deepvariant \
			-q 0 \
			--input_depth ~{samplename}.aligned.bam.depth \
			--input_vcf ~{vcf}

		# change fasta header name to just be the samplename
		sed "1s/.*/>~{samplename}/" ~{samplename}.vcfcons.fasta > ~{samplename}.fasta

		minimap2 \
			-a ~{reference} \
			~{samplename}.vcfcons.frag.fasta \
			> ~{samplename}.vcfcons.frag.fasta.sam

		samtools view \
			-bS ~{samplename}.vcfcons.frag.fasta.sam \
			> ~{samplename}.vcfcons.frag.fasta.bam

		samtools sort \
			~{samplename}.vcfcons.frag.fasta.bam \
			> ~{samplename}.vcfcons.frag.fasta.sorted.bam

		samtools index \
			~{samplename}.vcfcons.frag.fasta.sorted.bam

		# fix sample name and GQ field type (to be able to import via variant transforms)
		bcftools reheader --samples <(echo "~{samplename}") ~{samplename}.vcfcons.vcf -o ~{samplename}.sample_fixed.vcf
		bcftools view -h --no-version ~{samplename}.sample_fixed.vcf > header.txt
		sed -i 's/##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Conditional genotype quality">/##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Conditional genotype quality">/' header.txt
		bcftools reheader -h header.txt ~{samplename}.sample_fixed.vcf -o ~{samplename}.vcf

		bgzip ~{samplename}.vcf
		tabix ~{samplename}.vcf.gz
	}

	output {
		File high_quality_variants = "~{samplename}.vcf.gz"
		File high_quality_variants_index = "~{samplename}.vcf.gz.tbi"
		File consensus_fa = "~{samplename}.fasta"
	}

	runtime {
		docker: "~{container_registry}/cosa:0.0.1"
		cpu: 2
		memory: "7.5 GB"
		disks: "local-disk 100 HDD"
	}
}
