analyseCAGEtable=function(raw.dir="./", cage.dir="./", data)
{
 print('example data="0609_6_R1_t50"')
 print("loading files.")
 chr1=read.delim(paste(raw.dir,data,".chr1",sep=''),header=F)
 chr2=read.delim(paste(raw.dir,data,".chr2",sep=''),header=F)
 chr3=read.delim(paste(raw.dir,data,".chr3",sep=''),header=F)
 chr1=chr1[,1]
 chr2=chr2[,1]
 chr3=chr3[,1]
 chr1p=chr1[1:(length(chr1)/2)]
 chr1m=chr1[((length(chr1)/2)+1):(length(chr1))]
 chr2p=chr2[1:(length(chr2)/2)]
 chr2m=chr2[((length(chr2)/2)+1):(length(chr2))]
 chr3p=chr3[1:(length(chr3)/2)]
 chr3m=chr3[((length(chr3)/2)+1):(length(chr3))]
 rm(chr1,chr2,chr3)

 test=read.delim(paste(cage.dir,data,".5mis.ALLchr.col.CAGEtable",sep=''),stringsAsFactors=F,header=F)
#test=read.delim("0609_4_R1_t50.5mis.ALLchr.col.fused.filtered.paracluOUT.filtered.map.table",stringsAsFactors=F,header=F)
 colnames(test)=c("name","start","end","strand","chr","annot","c_name","c_chr","c_strand","c_start","c_end","sites","signal","min_dens","max_dens","hit2","hit1")
 test=test[,1:15]
#clust_length=test[,11]-test[,10]
 toDo=c(1:nrow(test))

 clust_length=test[toDo,11]-test[toDo,10]
 tss=vector()
 utr=vector()
 print("extracting tss and utr length.")
##
##extract TSS as the first point in the cluster equal to the highest hit per basepair value of the cluster.
##
 for (i in toDo)
 { 
  if((is.na(test[i,"c_name"]))|(test[i,"chr"] > 3))
  {
   tss[i]=NA
   utr[i]=NA
  }

  ##plus strand
  else if(test[i,"strand"] == "+")
  {
   if(clust_length[i] > 1)
   {
    chr=get(paste("chr",test[i,"chr"],"p",sep=''))
    a=chr[test[i,"c_start"]:test[i,"c_end"]]
    h=min(which(a==max(a)))
    tss[i]=test[i,"c_start"]+h
    utr[i]=test[i,"start"]-tss[i]
#   print(utr[i])
    rm(a,h)
   }
   else
   {
    tss[i]=test[i,"c_start"]
    utr[i]=test[i,"start"]-tss[i]
   }
  }

  ##minus strand
  else if (test[i,"strand"] == "-")
  {
   if(clust_length[i] > 0)
   {
    chr=get(paste("chr",test[i,"chr"],"m",sep=''))
    a=chr[test[i,"c_start"]:test[i,"c_end"]]
    h=max(which(a==max(a)))
    tss[i]=test[i,"c_start"]+h
    utr[i]=tss[i]-test[i,"end"]
#   print(utr[i])
    rm(a,h)
   }
   else
   {
    tss[i]=test[i,"c_start"]
    utr[i]=tss[i]-test[i,"end"]
   }
  }
 }
##
##
##

 test=data.frame(test[toDo,],clust_length,tss,utr,stringsAsFactors=F)
 colnames(test)=c(colnames(test)[1:15],"c_length","TSS","utr_length")

 map.out=read.delim(paste(cage.dir,data,".5mis.ALLchr.col.CAGEmap",sep=''),stringsAsFactors=F,header=F)

 tss.best=vector()
 utr.best=vector()
 clust_best=data.frame(stringsAsFactors=F)
 print("looking for closest clusters.")
###
#Assignes closest cluster to gene (may be usefull if read2 is not informative)
###
 for (i in toDo)
 {
#  print(i)
  chr=test[i,"chr"]
  st=test[i,"strand"]
  clust_best[i,"c_name"]=NA
  clust_best[i,"c_chr"]=NA
  clust_best[i,"c_strand"]=NA
  clust_best[i,"c_start"]=NA
  clust_best[i,"c_end"]=NA
  clust_best[i,"sites"]=NA
  clust_best[i,"signal"]=NA
  clust_best[i,"min_dens"]=NA
  clust_best[i,"max_dens"]=NA
  clust_best[i,"c_length"]=NA
  clust_best[i,"TSS"]=NA
  clust_best[i,"utr_length"]=NA

##
##
 
 if(test[i,"chr"] < 4)
  {
  
   totest=map.out[which((map.out[,2]==chr)&(map.out[,3]==st)),]
  
   if(test[i,"strand"]=="+")
   {
    c_best=test[i,3]-totest[,5]
##
    if(length(c_best)>0)
    {  
     if(length(c_best[which(c_best >= 0)]) > 0)
     {
      c_best=c_best[which(c_best >= 0)]
      c_best=min(c_best)
     }
     else
     {
      c_best=max((c_best))
     }
     tolink=which((test[i,3]-totest[,5])==c_best)
    }
   }
   else if (test[i,"strand"]=="-")
   {
    c_best=totest[,5]-test[i,3]
    if(length(c_best) > 0)
    {
     if(length(c_best[which(c_best >= 0)]) > 0)
     {
      c_best=c_best[which(c_best >= 0)]
      c_best=min(c_best)
     }
     else
     {
      c_best=max((c_best))
     }
     tolink=which((totest[,5]-test[i,3])==c_best)
    }
   }
#print(c_best) 
#print(tolink) 
   if(length(tolink) > 1)
   {
    print("####################")
    tolink=tolink[1]
   }
   if(length(tolink)!=1)
   {
    next
   }
   clust_best[i,"c_name"]=totest[tolink,1]
   clust_best[i,"c_chr"]=totest[tolink,2]
   clust_best[i,"c_strand"]=totest[tolink,3]
   clust_best[i,"c_start"]=totest[tolink,4]
   clust_best[i,"c_end"]=totest[tolink,5]
   clust_best[i,"sites"]=totest[tolink,6]
   clust_best[i,"signal"]=totest[tolink,7]
   clust_best[i,"min_dens"]=totest[tolink,8]
   clust_best[i,"max_dens"]=totest[tolink,9]
   clust_best[i,"c_length"]=totest[tolink,5]-totest[tolink,4]
   tolink=integer(0)
  
##plus strand
   if(test[i,"strand"] == "+")
   {
    if(clust_best[i,"c_length"] > 1)
    {
     chr=get(paste("chr",test[i,"chr"],"p",sep=''))
     a=chr[clust_best[i,"c_start"]:clust_best[i,"c_end"]]
     h=min(which(a==max(a)))
     clust_best[i,"TSS"]=clust_best[i,"c_start"]+h
     clust_best[i,"utr_length"]=test[i,"start"]-clust_best[i,"TSS"]
#   print(utr[i])
     rm(a,h)
    }
    else
    {
     clust_best[i,"TSS"]=clust_best[i,"c_start"]
     clust_best[i,"utr_length"]=test[i,"start"]-clust_best[i,"TSS"]
    }
   }
##minus strand
   else if (test[i,"strand"] == "-")
   {
    if(clust_best[i,"c_length"] > 0)
    {
     chr=get(paste("chr",test[i,"chr"],"m",sep=''))
     a=chr[clust_best[i,"c_start"]:clust_best[i,"c_end"]]
     h=max(which(a==max(a)))
     clust_best[i,"TSS"]=clust_best[i,"c_start"]+h
     clust_best[i,"utr_length"]=clust_best[i,"TSS"]-test[i,"end"]
     rm(a,h)
    }
    else
    {
     clust_best[i,"TSS"]=clust_best[i,"c_start"]
     clust_best[i,"utr_length"]=clust_best[i,"TSS"]-test[i,"end"]
    }
   }
  }
 }
##
##
##

 colnames(clust_best)=c("cb_name","cb_chr","cb_strand","cb_start","cb_end","cb_sites","cb_signal","cb_min_dens","cb_max_dens","cb_length","cb_TSS","cb_utr_length")

###
#Annotate table
###

 test=cbind(test,clust_best)
 print("computing some features.")

##Checks if clusters assigned based on read2 is also the cluster assigned based on distance alone (if true "c_same" = "YES").
 c_same=rep("NO",nrow(test))
 c_same[which(test[,"c_name"]==test[,"cb_name"])]="YES"
##Checks if cluster assigned on distance alone has been assiged to any gene based on read2 (if true "c_assigned" = "YES").
 c_assigned=rep("NO",nrow(test))
 c_assigned[which(test[,"cb_name"] %in% test[,"c_name"])]="YES"
 test=cbind(test,c_same,c_assigned)

##Checks whether a gene has been assigned more than one cluster by scan (dupl="YES"). In that case the gene has multiple lines.
 dupl=rep("NO",nrow(test))
 for (i in 1:nrow(test))
 {
 if (length(which(test[,1]==test[i,1])) > 1)
  {
   dupl[i]="YES"
  }
 }
 test=cbind(test,dupl)

##Flags bad ("dupl_bad"="YES") duplicated entries which are inside ORF if another entry for the same gene outside ORF exists.
 print("formatting summary table and removing in ORF duplicates.")
 d=unique(test[which(test[,"dupl"]=="YES"),1])
 dupl_bad=rep("NO",nrow(test))
 for (i in 1:nrow(test))
 {
  j=i-1
  if(j==0)
  {
   j=1
  }
  k=i+1
  if(k > nrow(test))
  {
   k=nrow(test)
  }
  if (test[i,1] %in% d)
  {
   if (test[j,1]==test[i,1])
   {
    if (is.na(test[i,"utr_length"])==T)
    {
     dupl_bad[i]="YES"
    }
    else if ((is.na(test[j,"utr_length"])==F)&((test[i,"utr_length"] < 0)&(test[j,"utr_length"] >= 0)))
    {
     dupl_bad[i]="YES"
    }
   }
   if (test[k,1]==test[i,1])
   {
    if (is.na(test[i,"utr_length"])==T)
    {
     dupl_bad[i]="YES"
    }
    else if ((is.na(test[k,"utr_length"])==F)&((test[i,"utr_length"] < 0)&(test[k,"utr_length"] >= 0)))
    {
     dupl_bad[i]="YES"
    }
   }
  }
 }

##Assigns a general flag based on all the above flags. Genes flaged "YES" are likelly to be worng (no scan TSS and clostest TSS inside ORF or faraway from ATG).
 bad=rep("NO",nrow(test))


 for (i in 1:nrow(test))
 {
  if (dupl_bad[i]=="YES")
  {
   bad[i]="YES"
  }
  if ((is.na(test[i,"cb_name"])==FALSE)|(is.na(test[i,"c_name"])==FALSE))
  {
   if((is.na(test[i,"c_name"])==T) & ((test[i,"cb_utr_length"] < 0)|(abs(test[i,"cb_utr_length"]) > 3000)))
   {
    bad[i]="YES"
   }
  }
  else
  {
   bad[i]="YES"
  }   
 }

##retrieves best TSS and utr in this order: good call assigned based on read2, or if not available good closest call. The scan column flags good calls assigned based on read2.
 best=rep(NA,nrow(test))
 best_utr=rep(NA,nrow(test))
 scan=rep("NO",nrow(test))
 for (i in 1:nrow(test))
 {
  if(bad[i]=="NO")
  {
   if(is.na(test[i,"TSS"])==FALSE)
   {
    best[i]=test[i,"TSS"]
    best_utr[i]=test[i,"utr_length"]
    scan[i]="YES"
   }
   else
   {
    best[i]=test[i,"cb_TSS"]
    best_utr[i]=test[i,"cb_utr_length"]
   }
  }
 }

 test=cbind(test,dupl_bad,bad,best,best_utr,scan)
 test1=test[,c(1:7,17:20,28:38)]

#print("3")
 assign(paste("all_",data,sep=''),test)
#print("4")
 assign(paste("short_",data,sep=''),test1)
#print("5")
 assign(paste("map_",data,sep=''),map.out)
#print("6")

save(list=c(paste("all_",data,sep=''),paste("short_",data,sep=''),paste("map_",data,sep='')),file=paste(cage.dir,data,".5mis.ALLchr.col.CAGEtable.rda",sep=''))
return(test1)
}




























