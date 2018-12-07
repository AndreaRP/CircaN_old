# CircaN

Nonlinear least squares model for accurate detection of circadian expression patterns.

## Setup
If you don't already have, download and install the devtools package.

```
install_github("AndreaRP/CircaN")
```

## Running CircaN

To help you get an idea of the type of data CircaN uses as input, we have included a toy dataset, 
along with it's metadata file in the package.

```
library("CircaN")
# Load expression data and metadata
expression_example <- CircaN::expression_example
metadata_example <-CircaN::metadata_example
# Now simply run the circan function
p <- circan(data=expression_example, s2c=metadata_example)
```
