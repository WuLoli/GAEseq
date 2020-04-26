# GAEseq: A Graph Auto-Encoder for Haplotype Assembly and Viral Quasispecies Reconstruction

Reconstructing components of a genomic mixture from data obtained by means of DNA sequencing is a challenging problem encountered in a variety of applications including single individual haplotyping and studies of viral communities. High-throughput DNA sequencing platforms oversample mixture components to provide massive amounts of reads whose relative positions can be determined by mapping the reads to a known reference genome; assembly of the components, however, requires discovery of the reads' origin -- an NP-hard problem that the existing methods struggle to solve with the required level of accuracy. In this paper, we present a learning framework based on a graph auto-encoder designed to exploit structural properties of sequencing data. The algorithm is a neural network which essentially trains to ignore sequencing errors and infers the posteriori probabilities of the origin of sequencing reads. Mixture components are then reconstructed by finding consensus of the reads determined to originate from the same genomic component. Results on realistic synthetic as well as experimental data demonstrate that the proposed framework reliably assembles haplotypes and reconstructs viral communities, often significantly outperforming state-of-the-art techniques.

Installation:
-----------------
1. Create a new folder and download the source code from repository
2. Enter Haplotype_assembly(or Viral_Quasispecies_Reconstruction) directory and run make (In terminal, go to Haplotype_assembly(or Viral_Quasispecies_Reconstruction) directory and make)

Software requirement:
-----------------
1. Python 3.6 or higher (https://www.python.org/)
2. Tensorflow package (Initial experiment uses version 1.14) (https://www.tensorflow.org/)
3. Biopython package (https://biopython.org/)
4. Numpy package (http://www.numpy.org/)
5. C++ (http://www.cplusplus.com/)

Data Preparation:
-----------------
Check config file format (configure to your setting)

* config file included in the package is configured to sample set using only CPU
* the aligned paired-end reads file (SAM format) and corresponding reference file (FASTA format) are required for haplotype assembly
* the order of rows in config file matters and the space + colon + space (" : ") to seperate variable name and value is required

Config file example:<br/> 
filename of reference sequence (FASTA) : K2_Cov5_Sample1.fasta  
filname of the aligned reads (sam format) : K2_Cov5_Sample1.sam<br/>
SNP_thres (0.5 / ploidy) : 0.25<br/>
reconstruction_start : 1<br/>
reconstruction_stop : 5000<br/>
min_mapping_qual : 60<br/>
min_read_length : 70<br/>
max_insert_length : 3000<br/>
characteristic zone name : K2_Cov5_Sample1<br/>
K (ploidy) : 2<br/>
Graph Convolutional layer number (must be even) : 2<br/>
Number of experiments : 200<br/>
dropout probability of read nodes to SNP nodes layer : 0<br/>
dropout probability of SNP nodes to read nodes layer : 0<br/>
dropout probability of dense layer : 0<br/>
Which GPU to use (set to -1 if only use CPU; otherwise set to 0, 1, 2 or 3 ... if multiple GPUs are available; check instructions for details) : -1<br/>

1. filename of reference sequence (FASTA) represents the name of the reference genome in .fasta or .fa format.
2. filname of the aligned reads (sam format) represents the name of the aligned paired-end reads file in .sam format.
3. SNP_thres represents the threshold for SNP calling. The lower it is, the more nucleotide positions would be considered as SNPs. (it's set to 0.5 / ploidy in initial experiments)
4. reconstruction_start represents the index in reference genome that is corresponding to the starting position of the genome region of interest (1 corresponds to the first nucleotide of the reference genome).
5. reconstruction_stop represents the index in reference genome that is corresponding to the ending position of the genome region of interest.
6. min_mapping_qual represents minimum read quality score (reads with score lower than it would be discarded).
7. min_read_length represents the minimum read length (reads shorter than it would be discarded).
8. max_insert_length represents the maximum insert length of reads (reads with insert longer than it would be discarded)
9. characteristic zone name can be set to any name of your choice.
10. K (ploidy) represents the ploidy.
11. Graph Convolutional layer number represents the number of graph convolutional layer you want to use and must to an even number.
12. Number of experiments represents the number of experiments for which GAEhaP would be implemented and the reconstructed haplotypes corresponding to the lowest MEC score would be the stored in haplotypes.txt file.  
13. dropout probability of read nodes to SNP nodes layer
14. dropout probability of SNP nodes to read nodes layer
15. dropout probability of dense layer
16. Which GPU to use represents the GPU you want to use. Use (from tensorflow.python.client import device_lib device_lib.list_local_devices()) to check the index of GPUs that are available
17. MEC improvement threshold(This is used for viral quasispecies reconstruction) : 0.09

Running GAEseq:
-----------------
1. First set up config file.
2. Command : <br/>
    (1) ./ExtractMatrix config<br/>
    (2) python GAEseq_haplo.py(or GAEseq_viral.py)<br/>
    Output : haplotypes.txt

**Alignment can be done via BWA MEM (http://bio-bwa.sourceforge.net/bwa.shtml) to generate the SAM file from FASTQ (reads) and FASTA (reference genome) files. SAM file can be sorted using "samtools sort".**

Reference:
-----------------
Ziqi Ke and Haris Vikalo. "A graph auto-encoder for haplotype assembly and viral quasispecies reconstruction," The 34th AAAI Conference on Artificial Intelligence (AAAI-20), New York, NY, February 7-12, 2020.<br/>
<br/>
Simulated dataset used in the paper can be downloaded at (https://drive.google.com/open?id=1v84CcCwCcKWc8zqsCKKT64-fMwLLFWBV)
<br/>
Realistic HIV-1 dataset used in the paper can be downloaded at (https://drive.google.com/file/d/1qMKQcFLsFegibTuZohKz8XpuTvPN8GSw/view?usp=sharing)
