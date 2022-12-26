# FWASS
FWASS repository computes the symmetric frequency-weighted allele-sharing similarity (FWASS) measure from Greenbaum, Gili, Alan R. Templeton, and Shirli Bar-David. "Inference and analysis of population structure using genetic data and network theory." Genetics 202, no. 4 (2016): 1299-1312.


## What this repo do?

This repo compute similarity between individuals based on genetic variation data.
It can compute either a frequency-weighted (Greenbaum et al., 2016),
or non frequency-weighted metric (Li and Horvitz, 1953). 
This computation method uses an assumption that the maximum number of alleles per site is low, and benefit from matrix multiplication fast algorithm to speed up the computation.


## How to use it?


### input files
You can currently use the FWASS with 2 file formats. VCF (compressed by gzip or not), or GenePop.

You don't need to specify your input file type.
If it ends with .csv or .xlsx we treat it as Genepop, and if it ends with .vcf or .vcf.gz we treat it as not compressed or compressed VCF file respectively.

The program is using multi-threading computations in order to compute faster on large files. This is able only with VCF files.

#### GenePop
Files should be in a cvs or xlsx formats. 
GenePop files can be only deal with a single thread, and should be able to be read as a whole in memory.
This is good for a small dataset, and not GWAS whole sequence datasets. 
For large datasets, use VCF format, or you might face with memory problems.

GenePop format is a table, which every line is an individual, and every column is a site.
First column should be named "ID", and under it should be the individual names. "ID" is an illegal name for a regular locus.

After 'ID' column, every column header is a locus name. The information in row x column y, is the alleles individual x holds at site y divided by '/'.
Do not use '/' in anywhere else, and allele names should not have '/' in it. 

#### VCF

Files should end with .vcf if not compressed and .vcf.gz if compressed with gzip.
Files can be arbitrary large, and are not load simultaneously in memory.
The only limitation we added to the general VCF format, is asking for every single line (locus)
format, to start with GT (genotype). We find it is common to all dataset we were working with,
and handling with the more general case will slow down the computation time.
Don't worry, you don't have to check it yourself! If there is a line which the format is not in the correct order,
the program will fail with a proper massage.

If you do have such a dataset where GT is not the first information in the line format, and the program failed,
you are welcome to talk with us (Greenbaum lab, The hebrew University of Jerusalem) and we will elaborate the program functionality.

### Requirements

This is a very simple repository, not using any special packages. The only non-standard packages it uses are:

Numpy, Pandas, tqdm

There is a requirements.txt file in the repository for specific versions. 
The repository supposed to work just fine with other versions as well.

#### Make sure to work with Python3. This was not tested with Python2!

### Running command line

If the repository is your working directory, you can execute a computation with the following command line:

`python3 compute_similarity -i <input file path> -o <output directory path> `

output  directory has to be a name that is not yet exist. The program will open this directory automatically.
If an exist directory name is given, the program will fail.

Optional flags:

  -w, --weighted: &emsp; If used, compute weighted metric, from Greenbaum et al., 2016. If not use, compute unweighted metric from Li and Horvitz, 1953
  
  --max_memo: &emsp; Max number of cells (individuals multiply by sites) to use in a single matrix (in millions).
  If you don't know, don't touch. If there are memory failures, reduce it. Default is 10
  
  --max_threads: &emsp; Maximum number of threads to compute small matrices simultaneously. Default is 8.
  
  --max_sites: &emsp; If assigned, compute similarity only based on the first n sites.

### output

There are 2 outputs of the program in the output directory given with -o:
1. similarity.csv (similarity_weighted.csv if used -w/--weighted). The result of the similarity computation. Every row is an individual, and every column is individual.
2. counts.csv. As well a matrix where every row and every column are individuals. The data in the cells are the number of sites the computation was based on. 
This differs between pairs of individuals because of missing data.


Good luck! If you encounter with problems, you can reach us in shaharmazia@gmail.com or gili.greenbaum@gmail.com