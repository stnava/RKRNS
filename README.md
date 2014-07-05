Within R:

```
library(devtools)
install_github('RKRNS','stnava')
```

Then check the help / vignette

```
help(package = "RKRNS", help_type = "html")
browseVignettes("RKRNS") # doesnt work yet
% render( paste(RKRNSsrcdir,"/vignettes/RKRNS.Rmd",sep='') ,'all')
```
