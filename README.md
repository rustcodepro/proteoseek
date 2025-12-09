# proteoseek

- a sequence and a annotated data approach to classify the PTM using the kmer and the annotated peak data.
-  association criteria: algorithmic view: Number of the kmers occurences from the previous
   sequences using the PTM generate and then the peak data clustering. If the value is there and 
   more than the threshold bumping up the cluster and if the value is less then lowering the cluster, 
   so that the binary classification can be configured.


```
cargo build
```

```
proteoseek for PTM protein machine learning classification based on the kmer abundance and the annotated peak
       ************************************************
       Gaurav Sablok
       codeprog@icloud.com
      ************************************************

Usage: proteoseek <COMMAND>

Commands:
  generate-ptm  Generate PTM classification
  ptm-classify  PTM machine learning
  help          Print this message or the help of the given subcommand(s)

Options:
  -h, --help     Print help
  -V, --version  Print version
```

Gaurav Sablok \
codeprog@icloud.com
