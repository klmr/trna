tRNA gene regulation downstream analysis
========================================

This repository contains the downstream analysis pipeline of the Schmitt,
Rudolph & al. paper. You can find information on how to run the pipeline or use
individual components in the following sections. Enjoy!

Please consult the online methods for more details.

Update
------

A previous version of the code plotted the rotation eigenvectors of the PCAs,
rather than the rotated data. The current version amends this. The changes are
highlighted in [commit 4860f48](commit/4860f48).

Dependencies
------------

This pipeline requires Python 2.7 and R 3.0. Dependencies on R packages are
listed [in the `DEPENDS` file](blob/master/DEPENDS).
Further dependencies are:

* Meme (4.9.0)
* [klmr/rcane][] (included as submodule)
* BioPython (>= 1.59)
* [klmr/pygoo][] (included as submodule)

[klmr/rcane]: https://github.com/klmr/rcane/tree/trna-project
[klmr/pygoo]: https://bitbucket.org/klmr/pygoo/src/?at=trna-project

After pulling the repository, please do

```shell
git submodule init
git submodule update
```

to initialise the submodules.

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

<!--
tree -d -L 3 | tail -n +2 | sed '$d' | tr -d '│─' | tr '├└' '*' | sed 's/[[:space:]]/ /g' | sed 's/^\( *\* \)\(.*\)/\1`\2`/'
-->

* `chip`
   * `plots`
      * [`colocalization`](#colocalization)
      * [`compensation`](#compensation)
      * [`correlation`](#correlation)
      * [`de`](#chip-de)
      * [`de-hist`](#chip-de-histograms)
      * [`distribution`](#chip-distribution)
      * [`replicates`](#chip-replicates)
      * [`usage`](#chip-usage)
      * [`usage-sampling`](#chip-usage-sampling)
   * `results`
      * [`active-genes`](#active-genes)
      * [`compensation`](#compensation)
      * `de-acc`
      * `de-genes`
      * `de-type`
      * [`meme`](#meme)
      * [`usage-sampling`](#chip-usage-sampling)
   * `scripts`
* `common`
   * `cache`
   * `data`
      * `meme`
   * `scripts`
       * `rcane`
* `rna`
    * `plots`
       * `correlation`
       * `de`
       * `distribution`
       * `usage`
       * `usage-sampling`
    * `results`
       * `de`
       * `usage-sampling`
    * `scripts`
        * `goo`

<!-- -->

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

The following contains a breakdown of the individual directories.

### Colocalization

> Figures are generated by `chip/scripts/colocalization.R`.

Each PDF file corresponds to one colocalisation test. The parameters for the
test are given in the title of the plot. The <i>p</i> value of the significance
test is labelled in the plot.

### Compensation

> Figures and tables are generated by `chip/scripts/compensation.R`.

`results/compensation/tested-isoacceptors.tsv` is a table of the anticodon
isoacceptor families for which compensation was tested, and their respective
(raw and adjusted) <i>p</i>-values.

### Correlation

### ChIP DE

### ChIP DE Histograms
