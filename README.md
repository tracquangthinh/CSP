# DCSP

Discovery of Druggable Cancer-Specific Pathway (DCSP)

For the interactive resource, please go to [DCSP website](https://www.meb.ki.se/shiny/truvu/DCSP/).

- `pas`: Source code for pathway activation score (PAS) generation.

- `shiny`: Source code for [DCSP website](https://www.meb.ki.se/shiny/truvu/DCSP/).

PAS is calculated for Genomics of Drug Sensitivity in Cancer (GDSC), The Cancer Genome Atlas (TCGA) and BeatAML. 
You can download PAS of GDSC from [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6787033.svg)](https://doi.org/10.5281/zenodo.6787033). For TCGA and BeatAML, we do not publish our results due to the limitation of storage space.
For further requests, please send an email to `quang.thinh.trac {at} ki.se`

To be able to access PAS of GDSC, you need to use `qs` package:

```
library(qs)
pas <- qread("/path/to/gdsc_pas.qs")
```



