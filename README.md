# K-mer Counter
K-mer counter counts the most frequent `n` `k-mers` in a given FASTQ file. This project has been developed and tested with Python 3.5.2.

## Getting Started

### Algorithms
K-mer counter uses two different algorithms:
1. [https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-333](BFCounter)
   BFCounter uses Bloom Filter for eliminating the unique k-mers in the FASTQ file. This operation reduces the number of items to be put in the hash table drastically. Thus, minimizes the memory usage.
2. [http://minia.genouest.org/dsk/](DSK)
   DSK uses hash values of the k-mers for grouping and partitioning the k-mers according to the specified target disk and memory spaces. K-mers are written in separate files, according to their hash value, for making the memory usage within the target value. For making the disk usage within the target value, this operation is done in multiple iterations.

### Prerequisites
K-mer counter uses Python 3.5.2, and usage of virtual environment is encouraged.
If you don't have virtual environment, you can use
```
pip install virtualenv
```
for installing the virtual environment module. Then create and activate the virtual environment.
```
virtualenv -p ${PYTHON_PATH} venv
source venv\bin\activate
```
For installing required modules, use following commands while virtual environment is activated.
```
pip install -r requirements.txt
```

#### Used Packages
1. **pybloomfiltermmap3:** Bloom filter implementation for Python 3
2. **mmh3:** MurmurHash implementation for Python
3. **progressbar2:** Progress bar used in verbose mode
4. **Cython:** Required by *pybloomfiltermmap3* module

There are also other packages used for code styling an linter purposes, such as
*flake8*, *pep8* etc.

## Running
Project can be run as follows
```
python kmer.py --file-name ${FASTQ_FILE} --kmer-size ${KMER_SIZE} --most-frequent ${MOST_FREQUENT} --error-rate ${ERROR_RATE} --target-disk ${TARGET_DISK} --target-memory ${TARGET_MEMORY} --algorithm ${ALGORITHM} --verbose
```

### Parameters
* **FASTQ_FILE:** Name of the .fastq file in which k-mers will be counted. (required)
* **KMER_SIZE:** K (required)
* **MOST_FREQUENT:** N (required)
* **ERROR_RATE:** Bloom Filter error rate, used only for *BFCounter* algorithm. (DEFAULT=0.001)
* **TARGET_DISK:** Target disk space that will be used in Gigabytes, used only for *DSK* algorithm. (DEFAULT=25)
* **TARGET_MEMORY:** Target memory space that will be used in Gigabytes, used only for *DSK* algorithm. (DEFAULT=4)
* **ALGORITHM:** For choosing the algorithm to be used. If the algorithm is not set (strongly recommended), the program will decide the algorithm to be used based on the criterion given in the [http://minia.genouest.org/dsk/](*DSK paper*).
  - `BF` or `bf` for *BFCounter*
  - `DSK` or `dsk` for *DSK*
* **VERBOSE:** For printing elapsed time and the hash table memory usage.

## Results
This program is tested with the `ERR047698.filt.fastq` and `ERR055763_1.filt.fastq` files, which can be found  [http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG01595/sequence_read/](here) for finding the most frequent 30 `25-mers`. Program has chosen BFCounter for *ERR047698.filt.fastq* and DSK for *ERR055763_1.filt.fastq*.

| File Name | Algorithm | Memory | Iteration | Partition | Duration |
| :---: | :---: | :---: | :---: | :---: | :---: |
| ERR047698.filt.fastq | BFCounter | 24.0MB | N\A | N\A | 25sec |
| ERR055763_1.filt.fastq | DSK | 1536.0MB | 4 | 9 | 2h18min |
