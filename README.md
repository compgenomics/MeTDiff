### MeTDiff

Version: 1.1.0

Date: 2014-12-18

Author: Xiaodong Cui <xiaodong.cui@outlook.com>, Jia Meng <jia.meng@hotmail.com>
  
  Maintainer: Xiaodong Cui <xiaodong.cui@outlook.com>
  
    The package is developed for the differential analysis for MeRIP-seq data of two experimental conditions to unveil the dynamics in post-transcriptional regulation of the RNA methylome. The MeTDiff R-package explicitly models the reads variation in data and also devices a more power likelihood ratio test for differential methylation site prediction (Xiaodong Cui, et al."MeTDiff: a Novel Differential RNA Methylation Analysis for MeRIP-Seq Data." Computational Biology and Bioinformatics, IEEE/ACM Transactions on  (Volume:PP ,Issue: 99)　). It accepts and statistically supports multiple biological replicates, internally removes PCR artifacts and multi-mapping reads, outputs exome-based binding sites (RNA methylation sites) and detects differential post-transcriptional RNA modification sites between two experimental conditions in term of percentage rather the absolute amount. The package is still under active development, and we welcome all biology and computation scientist for all kinds of collaborations and communications. Please feel free to contact Dr. Xiaodong Cui <xiaodong.cui@outlook.com> if you have any questions. Many thanks to Jia Meng for the exomePeak R package (Meng, Jia, et al. "Exome-based analysis for RNA epigenome sequencing data." Bioinformatics 29.12 (2013): 1565-1567.), which provides the base for MeTDiff. 

License: GPL-2

Depends: Rsamtools, GenomicFeatures (>= 1.0.0), rtracklayer


### Installation

`MeTDiff` depends on `Rsamtools`, `GenomicFeatures`, `rtracklayer` R packages and please make sure install them before installing `MeTDiff`
These three packages can be installed as:
  ```
source("http://bioconductor.org/biocLite.R")
biocLite("Rsamtools")
biocLite("GenomicFeatures")
biocLite("rtracklayer")
```
In addtion, R package `devtools` is required for MeTDiff to be installed from GitHub.
```
install.packages("devtools")
```
At last, `MeTDiff` can be installed as:
  
  ```
library("devtools")
install_github("compgenomics/MeTDiff")
```

### Toy Example
```
# using the data included in the package
library(MeTDiff)

# in the real case, change the gtf to what you need
gtf <- system.file('extdata','example.gtf',package='MeTDiff')

ip1 <- system.file('extdata','IP1.bam',package='MeTDiff')
ip2 <- system.file('extdata','IP2.bam',package='MeTDiff')
ip3 <- system.file('extdata','IP3.bam',package='MeTDiff')
input1 <- system.file('extdata','Input1.bam',package='MeTDiff')
input2 <- system.file('extdata','Input2.bam',package='MeTDiff')
input3 <- system.file('extdata','Input3.bam',package='MeTDiff')
treated_ip <- system.file('extdata','treated_IP1.bam',package='MeTDiff')
treated_input <- system.file('extdata','treated_Input1.bam',package='MeTDiff')

IP_BAM <- c(ip1,ip2,ip3)
INPUT_BAM <- c(input1,input2,input3)
TREATED_IP_BAM <- c(treated_ip)
TREATED_INPUT_BAM <- c(treated_input)

metdiff(GENE_ANNO_GTF=gtf,IP_BAM = IP_BAM,INPUT_BAM = INPUT_BAM,
        TREATED_IP_BAM = TREATED_IP_BAM,TREATED_INPUT_BAM=TREATED_INPUT_BAM,
        EXPERIMENT_NAME="example")
