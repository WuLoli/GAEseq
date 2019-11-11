GAEseq

Installation:
-----------------
1. Create a GAEseq directory and download the source code from repository
2. Enter GAEseq directory and run make (In terminal, go to GAEseq directory and make)

Software requirement:
-----------------
1. Python 3.6 or higher (https://www.python.org/)
2. Tensorflow package (Initial experiment uses version 1.10.0) (https://www.tensorflow.org/)
3. Biopython package (https://biopython.org/)
4. Numpy package (http://www.numpy.org/)
5. C++ (http://www.cplusplus.com/)

Data Preparation:
-----------------
Check config file format (configure to your setting)

* config file included in the package is configured to sample set using only CPU
* the aligned paired-end reads file (SAM format) and corresponding reference file (FASTA format) are required for haplotype assembly and viral quasispecies reconstruction
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
K (ploidy) or Initial population size : 2<br/>
Graph Convolutional layer number (must be even) : 2<br/>
Number of experiments : 200<br/>
dropout probability of read nodes to SNP nodes layer : 0<br/>
dropout probability of SNP nodes to read nodes layer : 0<br/>
dropout probability of dense layer : 0<br/>
Which GPU to use (set to -1 if only use CPU; otherwise set to 0, 1, 2 or 3 ... if multiple GPUs are available; check instructions for details) : -1<br/>
MEC improvement threshold : 0.09

1. filename of reference sequence (FASTA) represents the name of the reference genome in .fasta or .fa format.
2. filname of the aligned reads (sam format) represents the name of the aligned paired-end reads file in .sam format.
3. SNP_thres represents the threshold for SNP calling. The lower it is, the more nucleotide positions would be considered as SNPs. (it's set to 0.5 / ploidy in haplotype assembly and 0 for viral quasispecies reconstruction)
4. reconstruction_start represents the index in reference genome that is corresponding to the starting position of the genome region of interest (1 corresponds to the first nucleotide of the reference genome).
5. reconstruction_stop represents the index in reference genome that is corresponding to the ending position of the genome region of interest.
6. min_mapping_qual represents minimum read quality score (reads with score lower than it would be discarded).
7. min_read_length represents the minimum read length (reads shorter than it would be discarded).
8. max_insert_length represents the maximum insert length of reads (reads with insert longer than it would be discarded)
9. characteristic zone name can be set to any name of your choice.
10. K (ploidy) represents the ploidy or initial population size.
11. Graph Convolutional layer number represents the number of graph convolutional layer you want to use and must to an even number.
12. Number of experiments represents the number of experiments for which GAEseq would be implemented and the reconstructed haplotypes corresponding to the lowest MEC score would be the stored in haplotypes.txt file.  
13. dropout probability of read nodes to SNP nodes layer
14. dropout probability of SNP nodes to read nodes layer
15. dropout probability of dense layer
16. Which GPU to use represents the GPU you want to use. Use (from tensorflow.python.client import device_lib device_lib.list_local_devices()) to check the index of GPUs that are available.
17. MEC improvement threshold. The threshold used to estimate the number of clusters of a viral quasispecies.

Running GAEseq:
-----------------
1. First set up config file.
2. Command : ./ExtractMatrix config<br/>
             python GAEseq_haplo.py

Output : reconstructed haplotypes or viral quasispecies and inferred frequencies
