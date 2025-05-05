# Datasets

## Catalog

Datasets are large so they are on gitignore list. Decompress the .gz files before use.

NOTE: gplus_combined.txt.gz and soc-pokec-relationships.txt.gz are too large for GitHub, even compressed. Follow the links below to download.

| Dataset | Vertices (V) | Edges (E) | Density (E/V²) |
|---------|-------------|-----------|----------------|
| [amazon0601](https://snap.stanford.edu/data/amazon0601.html) | 403,394 | 3,387,388 | 0.00002082 |
| [gplus_combined](https://snap.stanford.edu/data/ego-Gplus.html) | 107,614 | 13,673,453 | 0.0011807 |
| [soc-pokec-relationships](https://snap.stanford.edu/data/soc-Pokec.html) | 1,632,803 | 30,622,564 | 0.0000115 |
| [twitter_combined](https://snap.stanford.edu/data/ego-Twitter.html) | 81,306 | 1,768,149 | 0.000267 |
| [web-BerkStan](https://snap.stanford.edu/data/web-BerkStan.html) | 685,230 | 7,600,595 | 0.0000162 |
| [web-Google](https://snap.stanford.edu/data/web-Google.html) | 875,713 | 5,105,039 | 0.00000667 |
| [web-NotreDame](https://snap.stanford.edu/data/web-NotreDame.html) | 325,729 | 1,497,134 | 0.0000141 |
| [web-Stanford](https://snap.stanford.edu/data/web-Stanford.html) | 281,903 | 2,312,497 | 0.0000291 |
| [wiki-Talk](https://snap.stanford.edu/data/wiki-Talk.html) | 2,394,385 | 5,021,410 | 0.000000877 |

Find more datasets at: https://snap.stanford.edu/data/index.html#web

## Downloading

To download gplus_combined.txt.gz and soc-pokec-relationships.txt.gz: `wget https://snap.stanford.edu/data/gplus_combined.txt.gz https://snap.stanford.edu/data/soc-pokec-relationships.txt.gz`

## Decompressing

To unzip all `.txt.gz` archives: `gunzip -k *.txt.gz`
