This is an R package to refactor the association testing code in GENESIS into something easier to maintain and optimize. The code in this package will eventually be incorporated back into GENESIS.

This package will contain the null model code and basic association testing functions, and will take only simple objects (data.frames for phenotypes, matrices for covariance and genotypes). A separate package will deal with preparing data (reading from GDS, etc) to pass to the functions in this package.
