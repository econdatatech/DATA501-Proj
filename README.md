# baselineARPLss

This is a repository for project of the DATA 501 class 2024 from Victoria University Wellington
https://www.wgtn.ac.nz/courses/data/501/2024/offering?crn=33170

The assignment is to create R package project distributed on Github.
The aim of the package will be to implement a semi-automatic baseline correction using asymmetrically reweighted penalized least squares smoothing [Baek et al., 2015].

The package can be installed from Github using:

```{r, eval=FALSE}
install.packages("devtools")
devtools::install_github("econdatatech/baselineARPLss") 
```

A Vignette is available from here: https://github.com/econdatatech/baselineARPLss/blob/main/vignettes/arPLS-vignette.pdf

A Testplan from here:
https://github.com/econdatatech/baselineARPLss/blob/main/TestPlan.pdf

[Baek et al., 2015]: Baek, S.-J., Park, A., Ahn, Y.-J., and Choo, J. (2015). Baseline correction
using asymmetrically reweighted penalized least squares smoothing. Analyst, 140:250–257.

The sample data in this package is (with the kind permission from Bob Downs) taken from https://rruff.info/repository/sample_child_record_raman_full/by_minerals/Abelsonite__R070007__Broad_Scan__532__0__unoriented__Raman_Data_RAW__13756.txt

Reference: Lafuente B, Downs R T, Yang H, Stone N (2015) The power of databases: the RRUFF project. 
In: Highlights in Mineralogical Crystallography, 
T Armbruster and R M Danisi, eds. Berlin, Germany, W. De Gruyter, pp 1-30


