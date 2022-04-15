# Downsample fastq

A set of scripts for randomly downsampling fastq files.
-

Randomly selects and outputs reads from single or paired read files based on either estimated genome size and desired fold-coverage or a given number of reads   

### **REQUIRED:** 

* perl  
* gzip  

### **USAGE:**

**downsample_paired**  
`perl downsample_paired.pl [options] \<read file 1 [fastq or fastq.gz]> \<read file 2 [fastq or fastq.gz]>`  
Proper pairing of reads is maintained in the output

**downsample_single**  
`perl downsample_single.pl [options] \<read file [fastq or fastq.gz]>`  

**Options:**  
  `-p`    Output prefix (default "downsampled")  
  `-g`    Estimated genome size, in Mbp (i.e. enter 3.5 for 3,500,000 bp)  
  `-f`    Desired average coverage, in fold (i.e. 100 for 100-fold coverage)  
  `-l`    Average length of reads in file1 and file2, separated by a comma (i.e. 101,99). If no values given, the average read length will be calculated (takes longer). If read lengths are given, be aware that the script does not check whether they are accurate.  
  `-r`    Number of reads to output per file (i.e. 1000 will output 1000 R1 reads and 1000 R2 reads). If -r is given, -g, -f, and -l will be ignored. If -r is not given, values must be entered for -g and -f  
  `-o`    Force output of all reads into new files if number of reads in original files is less than number of reads desired. (default is to not output any reads in this case)  
  
Output files will be gzipped.
  
