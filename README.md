Assumes you have cmake, git and the ability to compile code ( xcode, gcc, whatever is right for your system )

Within R:

```
library(devtools)
pckg = try(require(ANTsR))
if(!pckg) 
  {
  install_github("ANTsR",user="stnava")
  }
install_github("RKRNS",user="stnava")
```

Not the fastest or most updatable but few lines.

Then check the help / vignette:

```
help(package = "RKRNS", help_type = "html")
browseVignettes("RKRNS") #  may not yet work reliably
# render( paste(RKRNSsrcdir,"/vignettes/RKRNS.Rmd",sep='') ,'all')
```

If the vignette works, then great, b/c it tests basic functionality.
