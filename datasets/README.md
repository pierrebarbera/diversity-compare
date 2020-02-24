
# Drop-in directory for datasets.

## Expected Structure

Folder contains the preprocessing scripts, which have these assumptions about the structure of the dataset folders ($DF):

Path | Description
--- | ---
$DF/data/samples | per-sample (unaligned) fasta-files
$DF/data/reference.phylip | Reference MSA, phylip format
$DF/tree | raxml-ng treesearch result
$DF/tree/reference.raxml.bestTree | result tree of raxml-ng tree search
$DF/tree/reference.raxml.bestModel | result model of raxml-ng tree search

## Using the scripts

```<script>.sh <dataset folder name>```

Example:
```./otu_clustering.sh bv```

### Order of execution
```
otu_clustering.sh
chunkify.sh
align.sh
place.sh
unchunkify.sh
```