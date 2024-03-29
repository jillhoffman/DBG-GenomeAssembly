# De Bruijn Graphs for Genome Assembly

The repository implements a De Bruijn Graph (DBG) approach for genome assembly. The 
required inputs are two FATSA files, one of a target query and another of all reads.
The additional required input is the desired k-mer size to be used in DBG assembly.
This tool is useful for de-novo genome assembly and can be particularly relevant to
microbiome analysis.

## Setup
To install from git: 

    git clone https://github.com/jillhoffman/DBG-GenomeAssembly.git
    cd DBG-GenomeAssembly

Set up virtual environment and install requirements:    
*Please note that BioPython is only used for parsing of fasta files*

    virtualenv venv
    source venv/bin/activate
    pip install -r install/requirements.txt 

## Description of Scripts

**main.py:** Main script to run     
**makekmers.py:** Cuts strings into k-mers of desired size
**matchkmers.py:** Finds reads with >= 50% k-mer similarity to the reference read
**order.py:** Functions to put unique k-mers into the correct order for DBG assembly *(Note: still in development)*

## To Run

**Required Input Files:**   
The following inputs are required and are passed as arguments in the following order through command line.

* **query.fasta:** fasta file containing all target reads
* **query.fasta:** fasta file containing all reference reads
* **kmer_size** desired size of k-mers to be used 

Example data is provided in the *example_data* folder. This folder contains
a large *READS.fasta* file of all reference reads and smaller 
*READS.##.fasta* files, each containing ~5,000 sequences. The smaller files are recommended to 
use for testing as run time is only 17s, compared to 4 minutes for all sequences. Below shows statistics 
of identified reads based on different k-mer sizes. The recommended k-mer size for the example data is 4.

![](kmerstats.png)

**The example can be run by:**

```
python scripts/main.py ./example_data/QUERY.fasta ./example_data/READS.01.fasta  4
```
In this example the first set of 5,000 reads are being used and a 4-mers are being used. The "01" in *READS.01.fasta* 
can be switched with any number 01-23, or for the whole file *READS.fasta*

**Output:**
```
creating 4-mers...
matching 5184 reads to initial query...
matched 578 reads, the highest at (76.74418604651163%) kmer overlap
checking for more relevent reads...
found 106 additional reads overlapped
time elapsed: 17 seconds
```
Further into development, an output files describing the longest overlapped sequence will also be provided.

## Testing

Unittesting for the complete modules has been implemented. Manually curated ground truth examples were created for input to test
basic functionality of the makekmers.py and matchkmers.py functions. To run unit test script:
    
    cd ./scripts 
    python -m unittest unittesting.py

Output:

    Ran 2 tests in 0.021s
    OK

