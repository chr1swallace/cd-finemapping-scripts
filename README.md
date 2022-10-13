Scripts to run susie analysis on Crohn's disease CD_DeLange_28067908_1-hg38.tsv.gz downloaded from the genome catalog.

blocks.txt contains the lddetect block definitions in hg38

for each block b in blocks.txt, run

``` bash
./impute_block.R --args block=$b
./run_susie.R --args block=$b
```

Then collate results together

``` bash
Rscript ./collate.R
```

