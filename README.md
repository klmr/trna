tRNA gene regulation downstream analysis
========================================

This repository contains the downstream analysis pipeline of the Schmitt,
Rudolph & al. paper. You can find information on how to run the pipeline or use
individual components in the following sections. Enjoy!

Dependencies
------------

This pipeline requires Python 2.7 and R 3.0. Dependencies on R packages are
listed in the `DEPENDS` file. Further dependencies include:

* Meme (4.9.0)
* [klmr/rcane][] (…)
* BioPython (>= 1.59)
* [klmr/pygoo][] (…)

[klmr/rcane]: https://github.com/klmr/rcane/tree/…
[klmr/pygoo]: https://github.com/klmr/pygoo/tree/…

How to run
----------

Before running the pipeline, please ensure that all the necessary upstream data
is inside the `common/data` directory. This data can be obtained from the
supplementary material of the publication.

To re-generate all results, run the command

```shell
Rscript --vanilla common/scripts/generate-all.R
```

From within the base directory. This will take a while. Individual results can
be generated by running the respective script inside either the `rna` or `chip`
folder. For instance, to generate just the PCA plot for the ChIP data, execute
the following:

```shell
cd chip
Rscript --vanilla scripts/pca.R
```

Result files
------------

This section gives an overview over the result files and folders, with
instructions how to generate the file and how to interpret it.

### Outline

This is what the  directory tree of the project folder looks like:

    <tRNA project>
    ├── chip
    │   ├── plots
    │   │   ├── colocalization
    │   │   ├── compensation
    │   │   ├── correlation
    │   │   ├── de
    │   │   ├── de-hist
    │   │   ├── distribution
    │   │   ├── replicates
    │   │   └── usage
    │   ├── results
    │   │   ├── de-acc
    │   │   ├── de-genes
    │   │   ├── de-type
    │   │   └── meme
    │   └── scripts
    ├── common
    │   ├── cache
    │   ├── data
    │   └── scripts
    └── rna
        ├── plots
        │   ├── correlation
        │   ├── de
        │   ├── distribution
        │   └── usage
        ├── results
        │   └── de
        └── scripts

The **`*/scripts`** directories contain the source code of the pipeline contained in
this project.

The **`common/data`** directory contains the source data files from the
supplementary materials.

The **`common/cache`** directory contains cache files which speed up re-generation
of the results by caching some intermediate results. These files are *not
checked for staleness* – changing parameters or input data may require manually
deleting these files, otherwise results may not correctly update.

<!--
    Curious bug: replacing <code>…</code> with `…` in the following paragraph
    causes the intial ``**`plots`**`` to be rendered with verbatim `**` rather
    than in bold. GFM FTW.
-->

The **`{chip,rna}/plots`** directories contain result plots. Plots which
contrast tRNA ChIP data with RNA-seq derived data are found in the
<code>chip</code> directory. The <code>rna</code> directory contains only result
which are *solely* based on the RNA-seq data and not on ChIP-seq data.

The **`{chip,rna}/results`** directories contain any additional results which
are not in the form of plots. These are tables and raw text files, as well as
HTML files for the Meme analysis.
