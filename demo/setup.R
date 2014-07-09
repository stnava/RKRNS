## ---- ext1
print("#########setup#########TODO: need to incorporate motion parameters!")
istest<-FALSE
subject<-"111157"
datadir<-paste("/Users/stnava/data/KRNS/",subject,"/",sep='')
tr<-as.numeric(0.5)
responselength<-8/tr # e.g. 15 seconds div by 0.5 tr => 30 volumes
labs<-as.numeric(1:90) # label numbers to use ... need to know which label set is at hand
labs<-as.numeric( c(13,79,81,89) ) # lang network
throwaway<-8
ncompcor<-10
compcorvarval<-0.95
filterlowfrequency<-0.1 # 0.05 if you multiply by TR
filterhighfrequency<-0.5 # 0.4 # because of expected bold response < 25secs, > 5 seconds
trendfrequency<-3
winsorval<-0.01
eigsentbasislength<-100
aalfn<-paste(datadir,"aal/",subject,"_aal2.nii.gz",sep='')
if ( file.exists(aalfn) ) aalimg<-antsImageRead( aalfn , 3 )
bmaskfn<-paste(datadir,"ref/",subject,"_mask.nii.gz",sep='')
if ( file.exists(bmaskfn) ) bmask<-antsImageRead( bmaskfn , 3 )
reffn<-paste(datadir,"ref/",subject,"_mocoref.nii.gz",sep='')
if ( file.exists(reffn) ) ref<-antsImageRead( reffn , 3 )
imagedir<-paste(datadir,"moco/",sep='')
imagepostfix<-"_moco.nii.gz"
data("aal",package="ANTsR")
blocksCSVlist<-Sys.glob(paste(datadir,"design/*csv",sep=''))   # INPUT csv list
dfn<-paste(datadir,"assembly/assembled_design_",labs[1],"_",labs[length(labs)],"test.csv",sep='') # INPUT out csv name
afn<-paste(datadir,"assembly/assembled_aal_",labs[1],"_",labs[length(labs)],"test.mha",sep='') # INPUT out img na
