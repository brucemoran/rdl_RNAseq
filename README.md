# Analysis of Endometrial Carcinoma in BD Han II Rats

### Dependencies
The repo contains [NextFlow](https://www.nextflow.io/index.html#GetStarted) pipelines run inside [Singularity](https://sylabs.io/guides/3.5/user-guide/quick_start.html) containers.

This is relatively easily set up following respective installation instructions through the links.

### Steps to Reproduce Analysis
#### 1: Clone repo into desired location
```
github clone https://github.com/brucemoran/rdl_RNAseq && cd rdl_RNAseq
```

#### 2: Run RNAseq_STAR-kallisto_references.simg.nf to generate reference data
```
nextflow run refs.nf \
          --outDir "refs" \
          --set "rat" \
          -profile standard,refs
```

#### 3: Download sample data and generate sampleCsv input

#### 4: Run RNAseq_STAR-kallisto.simg.nf to generate output
```
nextflow run main.nf \
          --sampleCsv "data/sample.csv" \
          --projectBase "./" \
          --refsDir "refs" \
          --noStar \
          -profile standard
```
