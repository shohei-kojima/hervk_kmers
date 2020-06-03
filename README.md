# hervk_kmers

# About
Calculate k-mer occurrence from BAM file. Currently only BAM file is supported (does not support CRAM).

# git clone
git clone https://github.com/shohei-kojima/hervk_kmers


# Usage for impatient

### Calculate k-mer freq from all reads (e.g. all reads in BAM, including HERV-K)
```
python main.py -b tmp.bam -overwrite
```

### Calculate k-mer freq only from non-LTR region of HERV-K
```
python main_hervk.py -b tmp.bam -overwrite
```


# Parameters
Parameters can be changed by directly modifying a file `scripts/load_parameters.py`.

### self.max_mut [default=10]
If a read contained more than 10 (11 or more) mutations, this read will be discarded and NOT be used for analysis.

### self.max_clip_len [default=10]
If a read contained more than 10 (11 or more) soft cliped region(s), this read will be discarded and NOT be used for analysis.

### self.k [default=50]
Length of k of k-mer.

# Options
See `python main.py -h`.
