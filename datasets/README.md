# Datasets

## Catalog

Datasets are large so they are on gitignore list. Decompress the .gz files before use.

NOTE: soc-pokec-relationships.txt.gz is too large for GitHub, even compressed. Follow the link below to download.

| Dataset | Vertices (V) | Edges (E) | Density (E/V²) |
|---------|-------------|-----------|----------------|
| [amazon0601](https://snap.stanford.edu/data/amazon0601.html) | 403,394 | 3,387,388 | 0.00002082 |
| [soc-pokec-relationships](https://snap.stanford.edu/data/soc-Pokec.html) | 1,632,803 | 30,622,564 | 0.0000115 |
| [twitter_combined](https://snap.stanford.edu/data/ego-Twitter.html) | 81,306 | 1,768,149 | 0.000267 |
| [web-BerkStan](https://snap.stanford.edu/data/web-BerkStan.html) | 685,230 | 7,600,595 | 0.0000162 |
| [web-Google](https://snap.stanford.edu/data/web-Google.html) | 875,713 | 5,105,039 | 0.00000667 |
| [web-NotreDame](https://snap.stanford.edu/data/web-NotreDame.html) | 325,729 | 1,497,134 | 0.0000141 |
| [web-Stanford](https://snap.stanford.edu/data/web-Stanford.html) | 281,903 | 2,312,497 | 0.0000291 |
| [wiki-Talk](https://snap.stanford.edu/data/wiki-Talk.html) | 2,394,385 | 5,021,410 | 0.000000877 |
| [soc-Epinions1](https://snap.stanford.edu/data/soc-Epinions1.html) | 75,879 | 508,837 | 0.0000884 |

Find more datasets at: https://snap.stanford.edu/data/index.html#web

## Downloading

To download soc-pokec-relationships.txt.gz: `wget https://snap.stanford.edu/data/soc-pokec-relationships.txt.gz`

## Decompressing

To unzip all `.txt.gz` archives: `gunzip -k *.txt.gz`
