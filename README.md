# STRONG Taxa Finder
This code was developed to iterate through the bins created as output from [STRONG](https://github.com/chrisquince/STRONG "STRONG GitHub repo") to query the cog genes identified for each STRONG bin. The cog genes are queried against NCBI's ref_prok_rep_genomes via remote BLAST. Results for each cog queried as well as a bin summary file are produced. This code was written for and tested in Python 3.8.3.

# Execute Code:
To execute the code, you must specify two parameters:
1. Input Directory. This is the "results" directory created by STRONG. Either provide the full path or call the Python script from the "results" directory's parent directory.
2. Output Directory. This is the user specified directory name (or full path including directory name) for the results to be written. If the directory already exists, the script will append a "_1" to the directory name and will continue running.

```python parse_STRONG.py [-h] [-in INPUT_DIR] [-out OUTPUT_DIR]```

# Output:
Within the output directory, two files per bin will be created. The first is *.xml. This is the output from the remote BLAST. The second *_read_me.txt is a text file of the parsed BLAST results. The top two hits (meeting the threshold specified in the script [E-value < 1e-50]) for each cog is listed. Once all bins are processed, a summary file called "bin-compositions.txt" will be written to the directory. This lists the taxonomies identified for each bin with percentages reported. These percentages represent the % of cogs with a top hit to this taxa.
