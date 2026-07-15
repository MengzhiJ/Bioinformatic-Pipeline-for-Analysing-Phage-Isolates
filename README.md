# Bioinformatic Pipeline for Analysing Phage Isolates (Fan Lab)
**Last updated:** 15 July 2026, Mengzhi Ji, Xiangyu Fan
## 1. Genome Assembly
### 1.1 Quality Control and Adapter Trimming
**Software:** fastp (https://github.com/OpenGene/fastp)
```bash
$ fastp \
    -i ${sample}_R1.fastq.gz \
    -I ${sample}_R2.fastq.gz \
    -o ${sample}_R1.clean.fastq.gz \
    -O ${sample}_R2.clean.fastq.gz \
    -q 20 \
    -l 36 \
    -w 10
 ```
### 1.2 Assembly
**Software:** SPAdes (https://github.com/ablab/spades)
```bash
$ spades.py \
    -1 ${sample}_R1.clean.fastq.gz \
    -2 ${sample}_R2.clean.fastq.gz \
    --isolate \
    --careful \
    -t 10 \
    -o ${sample}_assembly
```
>If SPAdes fails to generate a complete genome assembly, alternative assemblers can be used.   

**Alternative Assembler:** shovill (https://github.com/tseemann/shovill)
```bash
$ shovill \
    --R1 ${sample}_R1.clean.fastq.gz \
    --R2 ${sample}_R2.clean.fastq.gz \
    --outdir ${sample}_assembly \
    --force \
    --cpus 10 \
    --assembler skesa
```
## 2. Genome Quality Evaluation
**Software:** CheckV (https://bitbucket.org/berkeleylab/checkv/src/master/)
```bash
$ checkv end_to_end \
    ${sample}_assembly_virus.fasta \
    ${sample}_assembly_virus.checkv \
    -t 10
```
## 3. Coverage Calculation (if needed, mainly for prophages induced from host strains)
**Software** 
- Bowtie2 (https://github.com/benlangmead/bowtie2) 
- Samtools (https://github.com/samtools/samtools) 
- CoverM (https://github.com/wwood/CoverM)
  
```bash
$ bowtie2-build \
    ${sample}_assembly_contigs.fasta \
    contigs \
    --large-index \
    --threads 10 
$ bowtie2 \
    -x contigs \
    -1 ${sample}_R1.clean.fastq.gz \
    -2 ${sample}_R2.clean.fastq.gz \
    -S ${sample}.sam \
    --no-unal \
    -p 10
$ samtools view -bS ${sample}.sam -@ 10 | \
    samtools sort -o ${sample}.sorted.bam -@ 10 -
$ TMPDIR=./ coverm contig \
    --bam-files ${sample}.sorted.bam \
    --min-read-percent-identity 95 \
    --min-read-aligned-percent 90 \
    -m trimmed_mean \
    --trim-min 10 \
    --trim-max 90 \
    -t 10 \
    -o ${sample}.coverage.tsv 
```
## 4. ANI Calculation and Comparison
**Software** 
- dRep (https://github.com/MrOlm/drep)
- FastANI (https://github.com/ParBLiSS/FastANI)
  
```bash
$ dRep dereplicate \
    ANI95_within_species \
    -g All_phages_reference_and_query/*.fasta \
    --ignoreGenomeQuality \
    -pa 0.9 \
    -sa 0.95 \
    -p 10
$ fastani \
    --ql query.txt \
    --rl reference.txt \
    -o fastani.txt \
    --matrix \
    --minFraction 0.1 \
    -t 10
```
## 5. Taxonomic Assignment using VICTOR based on the GBDP Algorithm
**Software:** VICTOR online (https://ggdc.dsmz.de/victor.php#)

>**Note:** The output includes a proteome-based phylogenetic tree of the query and reference genomes, along with the recommended taxonomic assignments at the genus and subfamily levels.

## 6. Genus/Family-Level AAI Comparison and Protein Sharing
**Software:** MGV script (https://github.com/snayfach/MGV/tree/master/aai_cluster)
```bash
$ diamond makedb \
    --in All_viruses.faa \
    --db viral_proteins \
    --threads 10
$ diamond blastp \
    --query All_viruses.faa \
    --db viral_proteins \
    --out blastp.tsv \
    --outfmt 6 \
    --evalue 1e-5 \
    --max-target-seqs 10000 \
    --query-cover 50 \
    --subject-cover 50
$ python amino_acid_identity.py \
    --in_faa All_viruses.faa \
    --in_blast blastp.tsv \
    --out_tsv aai.tsv
```
>**Note:** The `aai.tsv` file contains the AAI values and shared gene content for each genome pair. The results can be further processed using R or Python to generate a matrix for visualisation.

## 7. Functional Annotation and Genome Visualisation
>**Note:** There is no standard method for functional annotation, and no single approach is clearly superior to the others. However, considering the overall workflow and ease of use, I recommend the following two approaches.

### 7.1 Functional Annotation using Protein Language Models (PLMs) and Traditional Functional Databases
#### 7.1.1 VPF-PLM annotation
**Software:** VPF-PLM (https://github.com/kellylab/viral-protein-function-plm)
```bash
$ python scripts/embed_faa.py \
    -faa virus_protein.faa \
    -out PLM_annotation
$ python scripts/predict_function.py \
    -faa virus_protein.faa \
    -out PLM_annotation \
    --output_predictions \
    --prediction_heatmap \
    --output_embeddings \
    --efam_calibration_threshold
```
>**Note:** VPF-PLM assigns functional categories to viral proteins based on the PHROGs database. To obtain detailed functional annotations, the predicted proteins should be further searched against functional databases using sequence alignment tools.

#### 7.1.2 Functional Annotation by Sequence Alignment
**Software**
- Alignment tools: DIAMOND, BLASTP, MMseqs2, etc.
- Functional databases: NCBI nr, eggNOG, KEGG, COGDB, CYCDB, NCycDB, SCycDB, DRAM-v, etc.

**Recommended parameters**
- E-value <1e-05
- Bit score >50
- Sequence identity >40%
- Query coverage >60%

### 7.2 Functional Annotation using Pharokka, Phold and Phynteny
**Software**
- Pharokka (https://github.com/gbouras13/pharokka) 
- Phold (https://github.com/gbouras13/phold)
- Phynteny (https://github.com/susiegriggo/Phynteny_transformer)

**Recommended parameters**
- Official protocol provided by the developers (https://currentprotocols.onlinelibrary.wiley.com/doi/10.1002/cpz1.70405)

## 8. Genomic Comparison between Isolates and Closely Related Reference Phages
**Software:** Pharokka (v1.10.1, generating gbk file), Clinker (v0.0.32, for genomic comparison)
```bash
$ pharokka run \
    -i virus_id.fasta \
    -o pharokka_annotation \
    -d /Users/mengzhiji/database/pharokka_v1.8.0_database \
    -t 10 \
    -e 1E-05 \
    -g prodigal \
    --force
$ clinker gbk_files/*.gbk \
    -cm gene_function_color.csv \
    -i 0.4 \
    -j 10 \
    --force \
    -p plot_comparison.html 
```
## Citation
- Chen S. Ultrafast one‐pass FASTQ data preprocessing, quality control, and deduplication using fastp[J]. Imeta, 2023, 2(2): e107.
- Bankevich A, Nurk S, Antipov D, et al. SPAdes: a new genome assembly algorithm and its applications to single-cell sequencing[J]. Journal of computational biology, 2012, 19(5): 455-477.
- Souvorov A, Agarwala R, Lipman D J. SKESA: strategic k-mer extension for scrupulous assemblies[J]. Genome biology, 2018, 19(1): 153.
- Langmead B, Salzberg S L. Fast gapped-read alignment with Bowtie 2[J]. Nature methods, 2012, 9(4): 357-359.
- Li H, Handsaker B, Wysoker A, et al. The sequence alignment/map format and SAMtools[J]. bioinformatics, 2009, 25(16): 2078-2079.
- Aroney S T N, Newell R J P, Nissen J N, et al. CoverM: read alignment statistics for metagenomics[J]. Bioinformatics, 2025, 41(4): btaf147.
- Nayfach S, Camargo A P, Schulz F, et al. CheckV assesses the quality and completeness of metagenome-assembled viral genomes[J]. Nature biotechnology, 2021, 39(5): 578-585.
- Olm M R, Brown C T, Brooks B, et al. dRep: a tool for fast and accurate genomic comparisons that enables improved genome recovery from metagenomes through de-replication[J]. The ISME journal, 2017, 11(12): 2864-2868.
- Jain C, Rodriguez-R L M, Phillippy A M, et al. High throughput ANI analysis of 90K prokaryotic genomes reveals clear species boundaries[J]. Nature communications, 2018, 9(1): 5114.
- Meier-Kolthoff J P, Göker M. VICTOR: genome-based phylogeny and classification of prokaryotic viruses[J]. Bioinformatics, 2017, 33(21): 3396-3404.
- Buchfink B, Xie C, Huson D H. Fast and sensitive protein alignment using DIAMOND[J]. Nature methods, 2015, 12(1): 59-60.
- Nayfach S, Páez-Espino D, Call L, et al. Metagenomic compendium of 189,680 DNA viruses from the human gut microbiome[J]. Nature microbiology, 2021, 6(7): 960-970.
- Flamholz Z N, Biller S J, Kelly L. Large language models improve annotation of prokaryotic viral proteins[J]. Nature Microbiology, 2024, 9(2): 537-549.
- Bouras G, Nepal R, Houtak G, et al. Pharokka: a fast scalable bacteriophage annotation tool[J]. Bioinformatics, 2023, 39(1): btac776. 
- Bouras G, Grigson S R, Mirdita M, et al. Protein structure-informed bacteriophage genome annotation with Phold[J]. Nucleic Acids Research, 2026, 54(1): gkaf1448.
- Grigson S R, Bouras G, Papudeshi B, et al. Synteny-aware functional annotation of bacteriophage genomes with Phynteny[J]. bioRxiv, 2025: 2025.07. 28.667340.
- Bouras G, Grigson S R, Durr L, et al. Decoding viral dark matter: Metagenomic prokaryotic virus characterization with pharokka, phold, and phynteny[J]. Current Protocols, 2026, 6(7): e70405.
