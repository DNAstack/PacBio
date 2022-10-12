# PacBio tools and workflows

This repo contains DNAstack-managed tools and workflows used for processing PacBio data. Workflows are written in the [Workflow Description Language (WDL)](https://github.com/openwdl/wdl) and are containerized using [Docker](https://docs.docker.com/). Docker images are maintained in [DNAstack's container registry](https://hub.docker.com/u/dnastack).

Workflows found here are published in [DNAstack's Dockstore](https://dockstore.org/organizations/DNAstack).


## Workflows

See individual workflow directories for more information on workflow inputs and outputs.


### [Coronavirus (SARS-CoV-2) Sequencing Analysis (CoSA)](CoSA-CoronavirusSequencingAnalysis)

Process PacBio SARS-CoV-2 HiFi data to produce assembled viral genomes and VCFs.

### [Tandem Repeat Genotyper (TRGT)](TRGT-TandemRepeatGenotyper)

Perform targeted genotyping on tandom repeats from PacBio HiFi data.
