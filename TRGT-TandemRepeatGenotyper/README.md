# PacBio TRGT: Tandem Repeat Genotyper implemented in Workflow Description Language (WDL)

This repository contains a [WDL workflow](https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md) for processing PacBio HiFi data using the [TRGT tool](https://github.com/pacificBiosciences/trgt/). `trgt` profiles sequence composition, mosaicism, and CpG methylation of analyzed repeats.

Docker images containing the tools used by this workflow can be explored [in DNAstack's image repository](https://github.com/dnastack/bioinformatics-public-docker-images), or [on Dockerhub](https://hub.docker.com/u/dnastack).


## Workflow inputs

An input template file with some defaults predefined can be found [here](./workflows/inputs.json).
Some example input files can be found [in PacBio's `trgt` repository](https://github.com/PacificBiosciences/trgt/tree/main/example).

| Input | Description |
| :- | :- |
| `ref` | The reference genome that was used for read alignment (FASTA) |
| `aligned_bam`, `aligned_bai` | Aligned HiFi reads (BAM) and index (BAI) |
| `repeats` | The repeat definition file with reference coordinates and structure of tandem repeats (BED) |
| `repeat_id` | ID of the repeat to visualize |
| `container_registry` | Registry that hosts workflow containers. All containers are hosted in [DNAstack's Dockerhub](https://hub.docker.com/u/dnastack) [`dnastack`] |


## Workflow outputs

| Output | Description |
| :- | :- |
| `sorted_trgt_vcf`, `sorted_trgt_vcf_index` | Sorted VCF file and index that contains repeat genotypes; output by `trgt` |
| `sorted_trgt_bam`, `sorted_trgt_bam_index` | Sorted BAM file and index that contains pieces of HiFi reads that fully span the repeat sequences; output by `trgt` |
| `pileup_image` | An SVG file that contains the pileup read image; output by `trvz` |


## Running workflows

### Required software

- [Docker](https://docs.docker.com/get-docker/)
- [Cromwell](https://github.com/broadinstitute/cromwell/releases) & Java (8+) OR [miniwdl](https://github.com/chanzuckerberg/miniwdl/releases) & python3

### Running using Cromwell

From the root of the repository, run:

```bash
java -jar /path/to/cromwell.jar run /path/to/workflow.wdl -i /path/to/inputs.json
```

Output and execution files will be located in the `cromwell-executions` directory. When the workflow finishes successfully, it will output JSON (to stdout) specifying the full path to each output file.


### Running using miniwdl

This command assumes you have `miniwdl` available on your command line. If `miniwdl` is not available, try installing using `pip install miniwdl`.

```bash
miniwdl run /path/to/workflow.wdl -i /path/to/inputs.json
```

Output and execution files will be located in a dated directory (e.g. named `20200704_073415_main`). When the workflow finishes successfully, it will output JSON (to stdout) specifying the full path to each output file.


## Future work

* Add an optional alignment step so that FASTQ can be an input
* Improve workflow by changing `repeat_id` to grep the ID in the `repeat.bed` file instead of feeding it a literal single string in order to loop over TRVZ to generate multiple pile-up images
