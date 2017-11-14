# RNAditing

RNAditing is a software used for the RNA editing sites identifiaction on paired end or single end RNA-Seq.

It combine two open access softwares for the prediction and the research of RNA editing sites already known on databases:

- **RNAEditor**: This software is used for the prediction of editing sites on RNASeq experiments.  
Furthermore a list of "editing island" is returned. (Cluster of editing sites)

*John, D., Weirick, T., Dimmeler, S., & Uchida, S. (2016). RNAEditor: easy detection of RNA editing events and the introduction of editing islands. Briefings in bioinformatics, bbw087.*

https://github.com/djhn75/RNAEditor   


-  **REDItools**: This software is composed of different scripts implemented on Python. In our software only REDItoolsKnown.py is used.

This specific tool is used for the detection of BAM files editing sites. It make a comparison with specific editing sites databases for each site present on the BAM.
It explores known editing events in RNA-Seq experiments.

*Picardi, E., & Pesole, G. (2013). REDItools: high-throughput RNA editing detection made easy. Bioinformatics, btt287.*

*Picardi, E., D'Erchia, A. M., Montalvo, A., & Pesole, G. (2015). Using REDItools to detect RNA editing events in NGS datasets. Current protocols in bioinformatics, 12-12.*

This pipeline allowed to make analysis by population. In fact, comparisons beetween data samples of a same population is realised for keep only representative sites of this population.

In this README is explained:
- Requirements for our pipeline
- how to download the software and install it
- how to get files necessaries for our software
- how to use this software  

## Requirements

For this pipeline is important to create two specificals conda environnement:

### Environnement 1 (RNAEditor)   

**With conda installation:**   

- Python Packages

  - matplotlib
  - pyqt =4.11.4 (important: v4)
  - pysam =0.8.4 pre0
  - matplotlib - venn
  - pybedtools
  - mygene
  
- Softwares
  - blat =35
  - bwa =0.7.15
  - samtools =1.4
  - bedtools =2.26.0
  - trimmomatic
  - fastqc =0.11.4
  - java
  - GCC
  - Make   
  
**Without conda installation:**   
  
  - Picard tools   

    https://sourceforge.net/projects/picard/files/picard-tools/1.119/picard-tools-1.119.zip/
    ```
    $ unzip picard-tools-1.119.zip
    $ cd picard-tools-1.119
    $ mkdir /PATH/TO/.conda/envs/RNAEditor/bin/picard-tools/
    $ mv *.jar /PATH/TO/.conda/envs/RNAEditor/bin/picard-tools/
    ```

  - GenomeAnalyseToolKit3.5   

    https://software.broadinstitute.org/gatk/download/archive   
    
    Choose 3.5 version

    ```
    $ mkdir /PATH/TO/.conda/envs/RNAEditor/bin/GATK/
    $ mv GenomeAnalysisTK-3.5.tar.bz2/PATH/TO /.conda/envs/RNAEditor/bin/GATK/
    $ cd /PATH/TO/.conda/envs/RNAEditor/bin/GATK/
    $ bunzip2 GenomeAnalysisTK -3.5.tar.bz2
    $ tar -xvf GenomeAnalysisTK -3.5.tar 
    ```

  -  Samstat 1.5.1   

    https://sourceforge.net/projects/samstat/
    ```
    $ gunzip samstat -XXX.tar.gz
    $ tar -xvf samstat -XXX.tar
    $ cd samstat -XXX
    $ ./ configure
    $ make
    $ make check
    $ cp src/samstat/PATH/TO/.conda/envs/RNAEditor/bin/
    ```

  -  Homer v4.9   

    ```
    $ wget http://homer.ucsd.edu/homer/configureHomer.pl
    $ mkdir PATH/TO/Homer
    $ cd PATH/TO/Homer
    $ mv PATH/TO/configureHomer.pl
    $ perl configureHomer.pl -install
    $ perl configureHomer.pl -install hg19
    ```

    PATH to .bashrc   

    ```
    $ gedit PATH/TO/.bashrc
    $ PATH =" $PATH:PATH/TO/Homer/bin"
    ```

### Environnement 2 (REDItools)   

**With conda installation:**   

- Python Packages
  - pybedtools =0.6.9
  - numpy =1.12.1
  - mygene
  - pysam =0.6 (important : v0.6)

- Sofwares
  - samtools
  - matplotlib
  - matplotlib-venn   

**With conda installation:**   

  - psl_splicesites   

    ```
    $ wget http://research-pub.gene.com/gmap/src/gmap-gsnap-2017-04-13.tar.gz
    $ gunzip gmap-gsnap-2017-04-13.tar.gz
    $ tar -xvf gmap-gsnap-2017-04-13.tar
    $ cd gmap-2017-04-13
    $ ./ configure
    $ make
    $ make check
    $ make install
    $ cp util /psl_splicesites/PATH/TO/miniconda/envs/REDItools/bin/
    ```


## Download Pipeline

```$ git clone https://github.com/WYSNI/RNAditing.py```

## Pipeline Installation   

On .bashrc   

   ```
   $ export RNAditing_PATH_RNAEditor=/PATH/TO/RNAditing/RNAEditor/RNAEditor/
   $ export RNAditing_PATH_REDItools=/PATH/TO/RNAditing/REDItools/
   $ export RNAditing_PATH_RNAditing=/PATH/TO/RNAditing/RNAditing/
   $ alias RNAditing =" python /PATH/TO/RNAditing/RNAditing/RNAditing.py"
   $ alias TAB_modification="python /PATH/TO/RNAditing/RNAditing/TAB_modification.py"
   ```
## Download files

### For RNAEditor

  ```
  $ wget http://141.2.194.197/rnaeditor_annotations/GRCH37.tar.gz
  ```

### For REDItoolsKnown  

- Spicing Sites recovery
  ```
  $ wget http :// hgdownload .cse. ucsc .edu/ goldenPath / hg19 / database / refGene . txt.gz
  $ gunzip -c refGene .txt .gz | psl_splicesites -s 1 > mysplicesites
  $ gawk -F" " ’{ split ($2 ,a,":"); split (a[2] ,b,".");
  if (b[1] >b [3])
  print a[1] ,b[3] ,b[1] , toupper ( substr ($3 ,1 ,1)) ,"-";
  else
  print a[1] ,b[1] ,b[3] , toupper ( substr ($3 ,1 ,1)) ,"+"}’ mysplicesites
  > mysplicesites .ss
  ```

- Editing sites knows on database are used by REDItoolsKnown.

  Exist two major databases for the RNA Editing:

  + **Radar**:

  http://rnaedit.com/download/

  **Download editing sites knowns on Radar Database**

  All sites:

  ```
  wget http://lilab.stanford.edu/GokulR/database/Human_AG_all_hg19_v2.txt
  ```

  Download editing sites by region:

  ```
  wget http://lilab.stanford.edu/GokulR/database/Human_AG_Alu_hg19_v2.txt
  wget http://lilab.stanford.edu/GokulR/database/Human_AG_RepetitiveNonAlu_hg19_v2.txt
  wget http://lilab.stanford.edu/GokulR/database/Human_AG_NonRepetitive_hg19_v2.txt
  ```

  + **REDIportal:**

  http://srv00.recas.ba.infn.it/atlas/

  **Download editing sites known on REDIportal Database**

  ```
  wget http://srv00.recas.ba.infn.it/webshare/rediportalDownload/table1_full.txt.gz
  ```

  **Create files with editing sites by region:**

  ```
  gunzip -c REDIportal_full.txt.gz| awk '{if ($7 == "ALU"){print $0}}'>REDIportal_ALU.txt
  gunzip -c REDIportal_full.txt.gz| awk '{if ($7 == "NONREP"){print $0}}'>REDIportal_NONREP.txt
  gunzip -c REDIportal_full.txt.gz| awk '{if ($7 == "REP"){print $0}}'>REDIportal_REP.txt
  ```

  **File creation with 3 columns: Chromosome Position Strand**

  The REDItoolsKnown tool use files in a specific form.   
  So it's necessary to make modification of this file. The format is Chromosome Position Strand   
  Use our TAB_modification script for the modifications
  ```
  #EXAMPLE TAB_modification PATH/TO/file_in PATH/OUT/file_out
  $ TAB_modification PATH/TO/REDIportal_full.txt.gz PATH/OUT/REDIportal.tab
  ```
- dbSNP
  
    http://genome.ucsc.edu/cgi-bin/hgTables

## How to use RNAditing 

Exist two ways to use the pipeline

### All analysis

For populations analysis it's important to store samples files by population. Each population must be a folder containing its associated files.

For clarify analysis results, it's require to make a tabulate samples file used by the pipeline.
This file must contain 3 columns: 

Population, File PATH, sample name.   

Example:    

POP1 /PATH/TO/DATASAMPLE1 Sample1_POP1   

POP1 /PATH/TO/DATASAMPLE2 Sample2_POP1

```
$ source activate RNAEditor
$ 
# Example
$ RNAditing -A -a -c PATH/TO/configuration.txt -f PATH/TO/FASTQ_FOLDER/
-s PATH/TO/samples_files.txt -o /PATH/TO/BAM_FOLDER/ -p 0.75 -e paired
```
In this step:   
-  Samples editing sites prediction and representative editing sites for each population are identified.   
-  Editing sites annotation, edited genes recovery and graph creation are realized.
```
$ source deactivate
$ source activate REDItools
$
# Example
$ RNAditing -E -a -i PATH/TO/BAM_folder -d PATH/TO/database_file/
-r PATH/TO/Genome_ref.fasta -o PATH/TO/RESULTS_FOLDER/ -p 0.75 -c /PATH/TO/RNAEditor_Result_folder/
```

In this step:   
-  Samples editing sites discovery and representative editing sites for each population are identified. 
-  The common sites identified by the two programs are retained.   
-  De novo editing sites are retained.   
-  Editing sites annotation, edited genes recovery and graph creation are realized.

### Step by Step

For an utilisation of the pipeline step by step please check help of each RNAEditing Software.

```
$ RNAditing 
$ RNAditing -A
$ RNAditing -R
$ RNAditing -E
```


