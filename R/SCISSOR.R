## -----------------------------------------------------------------------------
## R functions from SCISSOR
#
# This re-implementation does not import or depend on the SCISSOR R package.
# The method details follow below:
# - Choi, H.Y., Jo, H., Zhao, X. et al. SCISSOR: a framework for identifying structural changes in RNA transcripts. Nat Commun 12, 286 (2021).
# - https://github.com/hyochoi/SCISSOR
## -----------------------------------------------------------------------------


#' @references https://github.com/hyochoi/SCISSOR/blob/master/R/get_Ranges.R
#' @noRd

get_Ranges <- function(Gene=NULL,regions=NULL,GTF.file=NULL,hg.ref=c("hg19","hg38"),outputType="part_intron") {

  if (!is.null(Gene)) {
    if (length(Gene)>1) {
      warning("More than one gene was provided. The first gene will be used.")
      Gene=Gene[1]
    }
  }
  if (is.null(regions)) {
    stop("regions must be specified")
    # if (is.null(Gene)) {
    # stop("Either Gene or regions must be specified")
    # } else {
    #   hg.ref = match.arg(hg.ref, choices=c("hg19","hg38"))
    #   regions = build_gaf(Gene=Gene,GTF.file=GTF.file,hg.ref=hg.ref)
    # }
  }
  chr = strsplit(regions,":")[[1]][1]
  strtend = do.call(rbind,strsplit(strsplit(strsplit(regions,":")[[1]][2],",")[[1]],"-"))
  strnd = strsplit(regions,":")[[1]][3]
  strtend.num=matrix(as.numeric(strtend),ncol=2)

  if (outputType=="whole_intron") {
    ep.new = find.exon.hy(regions,is.intron=TRUE,num.intron=NULL) ;
  } else if (outputType=="part_intron") {
    intron.len = ceiling(len.intron.hy(exon=regions)*0.5);
    ep.new = find.exon.hy(regions,is.intron=TRUE,num.intron=intron.len) ;
  } else if (outputType=="only_exon") {
    intron.len = ceiling(len.intron.hy(exon=regions)*0.5);
    ep.new = find.exon.hy(regions,is.intron=TRUE,num.intron=0) ;
  }
  lRanges0=ep.new$ep
  gRanges0=cbind(strtend.num[,1]-(lRanges0[,2]-lRanges0[,1]),strtend.num,strtend.num[,2]+(lRanges0[,4]-lRanges0[,3]))
  cRanges0=matrix(find.exon.hy(regions,is.intron=TRUE,num.intron=0)$ep[,c(2,3)],ncol=2)
  new.regions=paste0(chr,":",paste(apply(matrix(gRanges0[,c(1,4)],ncol=2),1,function(x) paste(x,collapse="-")),collapse=","),":",strnd)

  if (strnd=="+") {
    gRanges=gRanges0
    lRanges=lRanges0
    cRanges=cRanges0
  } else {
    gRanges=matrix(gRanges0[rev(1:nrow(gRanges0)),rev(1:4)],ncol=4)
    lRanges=matrix(max(lRanges0)-lRanges0[rev(1:nrow(lRanges0)),rev(1:4)]+1,ncol=4)
    cRanges=matrix(max(cRanges0)-cRanges0[rev(1:nrow(cRanges0)),2:1]+1,ncol=2)
  }
  rm(gRanges0,lRanges0,cRanges0)
  colnames(gRanges)=c("ip.start","e.start","e.end","ip.end") # e:exon; ip: intronic part;
  rownames(gRanges)=paste("exon",1:nrow(gRanges),sep="")
  colnames(lRanges)=c("ip.start","e.start","e.end","ip.end")
  rownames(lRanges)=paste("exon",1:nrow(gRanges),sep="")
  colnames(cRanges)=c("e.start","e.end")
  rownames(cRanges)=paste("exon",1:nrow(gRanges),sep="")

  output=list(Gene=Gene,gRanges=gRanges,lRanges=lRanges,cRanges=cRanges,chr=chr,strand=strnd,regions=regions,new.regions=new.regions)
  output
}


#' @references https://github.com/hyochoi/SCISSOR/blob/master/R/find.exon.hy.R
#' @noRd

find.exon.hy <- function(exon,is.intron=FALSE,num.intron=NULL){
  ## % First version: find.exon.hy.2; Second version: find.exon.hy.3;
  ## % Last updated: 1/14/2016
  ## % Find the locatons of the exon boundaries in plot
  ## % If there are only exons, then is.intron=FALSE
  ## % If the coverage file contains some intron, is.intron=TRUE and the number of intron is indicated as intron.num
  ## % If intron.num=NULL, whole introns are contained.
  a <- exon ;
  b <- strsplit(a,":")[[1]][2] ;
  c <- gsub("-",",",b) ;
  d <- as.numeric(strsplit(c,",")[[1]]) ;
  epl <- d[seq(1,length(d),by=2)] ;
  epr <- d[seq(2,length(d),by=2)] ;
  ep <- cbind(epl,epr) ;
  if (!is.intron){
    ep2 <- ep - ep[1,1] + 1;
    ep3 <- ep - ep[,1] + 1 ;
    if (nrow(ep3)>1){ new.ep <- ep3 + c(0,cumsum(ep3[1:(nrow(ep3)-1),2])) } else { new.ep <- ep3 }
    coverage.col <- c();
    for (iep in 1:nrow(ep2)){
      coverage.col <- c(coverage.col, ep2[iep,1]:ep2[iep,2])
    }
    return(list(coverage.col=coverage.col, epl=new.ep[,1], epr=new.ep[,2]))
  } else {
    ep2 <- ep - ep[1,1] + 1;
    if (is.null(num.intron)){
      nep2 = nrow(ep2)
      if (nep2>1) {
        ep2 <- cbind(ep2[,1],ep2,c(ep2[2:nep2,1]-1,ep2[nep2,2]))
      } else {
        ep2 <- matrix(c(1,ep2,ep2[2]),ncol=4)
      }
      new.ep <- ep2;
      # return(list(epl=new.ep[,1],epr=new.ep[,2]))
    } else {
      ep2 <- cbind(ep2[,1]-num.intron,ep2,ep2[,2]+num.intron); ep2[1,1] <- 1;
      if (nrow(ep2) > 1){
        overlap.id <- which(ep2[1:(nrow(ep2)-1),4]>ep2[2:nrow(ep2),1]) ;
        if (length(overlap.id)>0){
          ep2[overlap.id,4] <- ep2[overlap.id,3];
          ep2[(overlap.id+1),1] <- ep2[overlap.id,3] + 1;
        }
        ep2[nrow(ep2),4] <- ep2[nrow(ep2),3];
        ep3 <- ep2 - ep2[,1] + 1;
        new.ep <- ep3 + c(0,cumsum(ep3[1:(nrow(ep3)-1),4]));
        overlap.id.2 <- which(new.ep[1:(nrow(new.ep)-1),4]>new.ep[2:nrow(new.ep),2]) ;
        if (length(overlap.id.2)>0){
          print(paste("!!!! Warning: Some exons and introns are overlapped!!!!  Overlapped at:  ",overlap.id.2))
        }
      } else {
        ep2[nrow(ep2),4] <- ep2[nrow(ep2),3];
        new.ep <- ep2
      }
    }
    #	Positions will be used for introns and exons
    coverage.col <- c();
    for (iep in 1:nrow(ep2)){
      coverage.col <- c(coverage.col, ep2[iep,1]:ep2[iep,4])
    }
    return(list(coverage.col=coverage.col,epl=new.ep[,2],epr=new.ep[,3],ep=new.ep));
    #	coverage.col = bp positions that will be used among the whole transcript (exon+intron)
    #	epl, epr = left and right exon positions in coverage.col

  }
}


#' @references https://github.com/hyochoi/SCISSOR/blob/master/R/len.intron.hy.R
#' @noRd

len.intron.hy <- function(exon){
  # 2nd version. (undated on 4/2/2017)
  # Input: exon info.
  # Produces the length of introns whose total length is equal to the total length of exons
  a <- exon ;
  b <- strsplit(a,":")[[1]][2] ;
  c <- gsub("-",",",b) ;
  d <- as.numeric(strsplit(c,",")[[1]]) ;
  epl <- d[seq(1,length(d),by=2)] ;
  epr <- d[seq(2,length(d),by=2)] ;
  ep <- cbind(epl,epr) ;

  n.ex <- nrow(ep);
  if (n.ex==1){
    return(0)
  } else {
    n.int <- nrow(ep)-1;
    len.int <- epl[2:n.ex]-epr[1:(n.ex-1)] - 1 ; # Intron lengths between exons
    srt.len.int <- sort(len.int) ; # sorted intron lengths between exons
    ttl.ex <- sum(epr - epl + 1) ; # total length of exons
    mean.int <- floor(ttl.ex/n.int) ; # Expected mean of intron lengths

    if (srt.len.int[1]>mean.int){
      output <- mean.int ;
    } else if (srt.len.int[n.int] <= mean.int){
      output <- srt.len.int[n.int] ;
    } else if (sum(len.int)<=ttl.ex){
      output <- srt.len.int[n.int];
    } else {
      jp <- 0; jc <- 1;
      while (jc > jp){
        jp <- jc ;
        crr.ttl.int <- ttl.ex - sum(srt.len.int[1:jc]) ; # total lengths of introns at current step
        crr.mean.int <- floor(crr.ttl.int/(n.int-jc)) ;

        if (srt.len.int[jc+1]>=crr.mean.int){
          output <- crr.mean.int ;
          jc <- jp ;
        } else {
          jc <- jp + 1 ;
        }
      }
    }
    return(output) ;
  }
}


#' @references https://github.com/hyochoi/SCISSOR/blob/master/R/build_pileup.R
#' @noRd

build_pileup <- function(Pileup,caseIDs=NULL,regions,
                        inputType="part_intron",
                        outputType="part_intron") {
  # inputType  = "whole_intron", "part_intron", "only_exon"
  # output.ytpe = "whole_intron", "part_intron", "only_exon"

  if (missing(Pileup)) {
    stop("Pileup is missing")
  }
  if (missing(regions)) {
    stop("Regions should be specified")
  }

  strnd = strsplit(regions,":")[[1]][3]
  if (strnd=="-") {
    rawPileup = Pileup[rev(1:nrow(Pileup)),]
  } else {
    rawPileup = Pileup
  }
  if (is.null(caseIDs)) {
    caseIDs = paste0("case-",1:ncol(rawPileup),sep="")
  }
  if (inputType=="whole_intron") {

    if (outputType=="whole_intron") {
      intron.len = NULL
      ep.new = find.exon.hy(regions,is.intron=TRUE,num.intron=intron.len) ;
      covPileup = rawPileup;
    } else if (outputType=="part_intron") {
      intron.len = ceiling(len.intron.hy(regions)*0.5);
      ep.new = find.exon.hy(regions,is.intron=TRUE,num.intron=intron.len) ;
      covPileup = rawPileup[ep.new$coverage.col,] ;    #   Area to be included.
    } else if (outputType=="only_exon") {
      intron.len = 0;
      ep.new = find.exon.hy(regions,is.intron=TRUE,num.intron=intron.len) ;
      covPileup = rawPileup[ep.new$coverage.col,] ;    #   Area to be included.
    } else {
      stop(outputType," is not an option for outputType.")
    }

  } else if (inputType=="part_intron") {

    if (outputType=="whole_intron") {
      stop(outputType," is not an option when inputType=part_intron.")
    } else if (outputType=="part_intron") {
      intron.len = ceiling(len.intron.hy(regions)*0.5);
      ep.new = find.exon.hy(regions,is.intron=TRUE,num.intron=intron.len) ;
      covPileup = rawPileup ;    #   Area to be included.
      ##### If the dimension of the given pileup data is not equal to the expected one
      ##### from the intron.len calculated, should display an error. Fix this!
    } else if (outputType=="only_exon") {
      intron.len.temp = ceiling(len.intron.hy(regions)*0.5);
      ep.new.temp = find.exon.hy(regions,is.intron=TRUE,num.intron=intron.len.temp) ;

      exonic.region = c()
      for (i in 1:nrow(ep.new.temp$ep)) {
        exonic.region = c(exonic.region,(ep.new.temp$epl[i]:ep.new.temp$epr[i]));
      }
      covPileup = rawPileup[exonic.region,];
      rm(intron.len.temp,ep.new.temp,exonic.region);

      intron.len = 0;
      ep.new = find.exon.hy(regions,is.intron=TRUE,num.intron=intron.len) ;
    } else {
      stop(outputType," is not an option for outputType.")
    }

  } else if (inputType=="only_exon") {

    if (outputType=="whole_intron") {
      stop(outputType," is not an option when inputType==only_exon.")
    } else if (outputType=="part_intron") {
      stop(outputType," is not an option when inputType==only_exon.")
    } else if (outputType=="only_exon") {
      intron.len = 0;
      ep.new = find.exon.hy(regions,is.intron=TRUE,num.intron=intron.len) ;
      covPileup = rawPileup ;
    } else {
      stop(outputType," is not an option for outputType.")
    }

  } else {
    stop(inputType," is not an option for inputType.")
  }

  Ranges=get_Ranges(regions=regions,outputType=outputType)
  new.regions=Ranges$new.regions
  chr = strsplit(new.regions,":")[[1]][1]
  strtend = do.call(rbind,strsplit(strsplit(strsplit(new.regions,":")[[1]][2],",")[[1]],"-"))
  strnd = strsplit(new.regions,":")[[1]][3]
  strtend.num=matrix(as.numeric(strtend),ncol=2)
  allPos = unlist(sapply(1:nrow(strtend.num), function(x) strtend.num[x,1]:strtend.num[x,2]))

  rownames(covPileup) = allPos
  colnames(covPileup) = caseIDs
  if (strnd=="+") {
    return(covPileup)
  } else {
    return(covPileup[rev(1:nrow(covPileup)),])
  }
}


#' @references https://github.com/hyochoi/SCISSOR/blob/master/R/process_pileup.R
#' @noRd

process_pileup <- function(pileupData,Ranges,
                          logshiftVal=NULL,
                          plotNormalization=FALSE, ...) {

  # Log transform data
  data.log = logtransform_data(inputData=pileupData,logshiftVal=logshiftVal)

  # Center and normalize data
  data.normalized = normalize_data(inputData=data.log$outputData,
                                   pileup=pileupData,Ranges=Ranges,
                                   makePlot=plotNormalization, ...)

  readconstr=20
  exonset=Ranges$lRanges
  exon.rm=c(); intron.rm=c();
  if (nrow(exonset)>1) {
    for (ie in 1:nrow(exonset)) {
      b.exon=c(exonset[ie,2]:exonset[ie,3])
      if (max(apply(pileupData[b.exon,],1,max))<readconstr) {
        exon.rm=c(exon.rm,ie)
      }
    }

    for (ie in 1:(nrow(exonset)-1)) {
      b.exon1=c(exonset[ie,2]:exonset[ie,3])
      b.exon2=c(exonset[ie+1,2]:exonset[ie+1,3])
      q.exon=apply(cbind(apply(pileupData[b.exon1,],2,FUN=function(x){quantile(x,probs=0.5)}),
                         apply(pileupData[b.exon2,],2,FUN=function(x){quantile(x,probs=0.5)})),1,max)

      b.intron=c((exonset[ie,3]+6):(exonset[ie+1,2]-6)) # remove each boundary of length 4
      lcarea=length(b.intron)
      if (lcarea<22) { # remove introns less than 30 (22+4+4)
        intron.rm=c(intron.rm,ie)
      } else {
        q.intron=apply(cbind(apply(pileupData[b.intron[6:ceiling(lcarea/3)],],2,FUN=function(x){quantile(x,probs=0.75)}),
                             apply(pileupData[b.intron[ceiling(lcarea/3):(2*ceiling(lcarea/3))],],2,FUN=function(x){quantile(x,probs=0.75)}),
                             apply(pileupData[b.intron[(2*ceiling(lcarea/3)):(lcarea-6)],],2,FUN=function(x){quantile(x,probs=0.75)})),1,max)

        if ((max(apply(pileupData[b.intron,],1,max))<readconstr) &
            (length(which((q.intron>0.2*q.exon) & (q.exon>10)))==0)) {
          intron.rm=c(intron.rm,ie)
        }
      }
      # cat(ie,"\n")
    }

    new.exonset=exonset
    for (ie in 1:(nrow(exonset)-1)) {
      if (exonset[ie+1,2]==exonset[ie+1,1]) {
        baselen=floor(0.5*(exonset[ie+1,1]-exonset[ie,3]))
        new.exonset[ie,4]=exonset[ie,3]+baselen
        new.exonset[ie+1,1]=exonset[ie,3]+baselen+1
      }
    }
    # new.exonset

    if (length(exon.rm)>0) {
      for (ie in exon.rm) {
        data.normalized$outputData[c(new.exonset[ie,1]:new.exonset[ie,4]),]=0
      }
    }
    if (length(intron.rm)>0) {
      for (ie in intron.rm) {
        data.normalized$outputData[c((exonset[ie,3]+1):(exonset[ie+1,2]-1)),]=0
      }
    }
  }

  return(list(logData=data.log$outputData,
              normalizedData=data.normalized$outputData,
              logshiftVal=data.log$logshiftVal,
              msf=data.normalized$msf,
              g1.offset=data.normalized$g1.offset,g2.offset=data.normalized$g2.offset))
}


#' @references https://github.com/hyochoi/SCISSOR/blob/master/R/logtransform_data.R
#' @noRd

logtransform_data <- function(inputData, logshiftVal=NULL, param.grid=NULL, draw.plot=FALSE) {
  # Find exon bases out of intron-contained coverage
  #     & Take log-transformation.
  if (is.null(logshiftVal)) {
    logshiftVal = getShiftParam(X=inputData, param.grid=param.grid, draw.plot=draw.plot)$optim.param
  }
  outputData = log10(inputData + logshiftVal) - log10(logshiftVal) ;
  return(list(outputData=outputData,logshiftVal=logshiftVal))
}


#' @references https://github.com/hyochoi/SCISSOR/blob/master/R/getShiftParam.R
#' @noRd

getShiftParam <- function(X, param.grid=NULL, draw.plot=FALSE) {

  if (is.null(param.grid)){
    param.grid = seq(1,10,by=1)
  }
  res = apply(matrix(param.grid,ncol=1),1,
              FUN=function(x){SlogADstat.hy(X=X,shift.val=x)$ADval})
  optim.idx = which.min(res)
  optim.param = param.grid[optim.idx]
  optim.stat = res[optim.idx]

  if (is.infinite(optim.stat)) {
    warning("Minimum = Inf")
  }

  if (draw.plot){
    plot(param.grid,res,type="o",ylab="A-D statistic",xlab="Shift value")
    abline(v=optim.param,col="red");
    text(x=optim.param+(0.2*(max(param.grid)-min(param.grid))),
         y=(min(res)+0.8*(max(res)-min(res))),labels=round(optim.param,digits=3),col="red")
  }
  return(list(optim.param=optim.param,optim.stat=optim.stat))
}


#' @references https://github.com/hyochoi/SCISSOR/blob/master/R/SlogADstat.hy.R
#' @noRd

SlogADstat.hy <- function(X,shift.val=1){
  # log transform
  datalog = log10(X + shift.val) - log10(shift.val) ;
  # calculate the median scale factor
  msf = scale.factor(X=datalog,average="mean",trim=0.1,adjval=NULL)$msf

  ADval = ADstatWins.hy(msf);
  return(list(msf=msf,ADval=ADval));
}


#' @references https://github.com/hyochoi/SCISSOR/blob/master/R/ADstatWins.hy.R
#' @noRd

ADstatWins.hy <- function(x,trim=NULL){
  # Winsorize and calculate Anderson Darling test statistic
  if (!is.null(trim)) {
    n = length(x)
    m = floor(0.5*trim*n)
    if (m>0) {
      y = x[order(x,decreasing=F)[-c(1:m)]]
      x = y[order(y,decreasing=T)[-c(1:m)]]
    }
  }
  n = length(x);
  mabd = mean(abs(x-median(x)));

  if (mabd==0){
    xfinal = rep(0,n);
  } else {
    xk = (x-median(x))/mabd ;
    # Winsorise data
    p95 = 1.0972;   # Extreme value inverse with p=95%, mu=0,sigma=1
    an = (2*log(n))^(-0.5);
    bn = sqrt(2*log(n)-log(log(n))-log(4*pi));
    L = p95*an+bn;
    xk[which(xk>L)] = rep(L,length(which(xk>L)));
    xk[which(xk<(-L))] = rep((-L),length(which(xk<(-L))));

    # Re-standardize
    xfinal = (xk-mean(xk))/sd(xk);
  }

  # Calculate evaluation stat
  ADval = ADstat.hy(xfinal);
  return(ADval);
}


#' @references https://github.com/hyochoi/SCISSOR/blob/master/R/ADstat.hy.R
#' @noRd

ADstat.hy <- function(x){
  n = length(x);
  if (n < 7){
    return(print("Sample size must be greater than 7"));
  } else {
    xs = sort(x);
    f = pnorm(xs,mean(xs),sd(xs));
    i = 1:n ;
    S = sum(((2*i-1)/n)*(log(f)+log(1-rev(f))));
    ADval = -n-S;
    return(ADval);
  }
}


#' @references https://github.com/hyochoi/SCISSOR/blob/master/R/scale.factor.R
#' @noRd

scale.factor <- function(X,average="mean",trim=0.1,adjval=NULL,robustLM=FALSE) {
  mean.vec = apply(X,1,FUN=function(x){adj.center(x,average=average,trim=trim,adjval=adjval)})
  if (sum(mean.vec^2)>0) {
    if (robustLM) {
      msf = apply(X,2,FUN=function(x){MASS::rlm(x ~ mean.vec - 1,psi=psi.bisquare)$coefficients})
      if (length(which(is.na(msf)==1))>0) {
        temp_id=which(is.na(msf))
        msf[temp_id]=as.vector(t(mean.vec)%*%X[,temp_id]/(sum(mean.vec^2)))
      }
    } else {
      msf = as.vector(t(mean.vec)%*%X/(sum(mean.vec^2)));
    }
  } else {
    msf = rep(0,ncol(X))
  }
  NewX = X - mean.vec%*%t(msf)
  return(list(msf=msf,mean.vec=mean.vec,NewX=NewX))
}


#' @references https://github.com/hyochoi/SCISSOR/blob/master/R/adj.center.R
#' @noRd

adj.center <- function(x,average="mean",trim=0.1,adjval=NULL) {
  ##  % Obtain adjusted center (mean or median)
  if (is.null(adjval)) {
    if (average=="median") {
      center.val = median(x);
    } else {
      center.val = mean(x,trim=trim);
    }
  } else {
    nonzero.pos = which(x>adjval)
    if (length(nonzero.pos)>0) {
      if (average=="median") {
        center.val = median(x[nonzero.pos])
      } else {
        center.val = mean(x[nonzero.pos])
      }
    } else {
      center.val = 0
    }
  }
  return(center.val);
}


#' @references https://github.com/hyochoi/SCISSOR/blob/master/R/normalize_data.R
#' @noRd

normalize_data <- function(inputData,
                          pileupData,Ranges,
                          smoothness=0.7,
                          makePlot=FALSE, ...) {

  centerDataResult = center_data(inputData=inputData,Ranges=Ranges)
  exonset = Ranges$lRanges
  Gene = Ranges$Gene
  g1.offset = estimate_offset(centerDataResult=centerDataResult,
                              msf=NULL,centeredData=NULL,
                              pileupData=pileupData,
                              exonset=exonset,
                              smoothness=smoothness,makePlot=FALSE)
  msf = centerDataResult$msf
  normalizedData = sweep(x=centerDataResult$outputData,2,g1.offset$g,FUN="/")
  goodcase = g1.offset$goodcase
  if (makePlot) {
    if (is.null(Gene)) {
      plot_offset(offset.obj=g1.offset,draw.legend=T,
                  main="Before normalization",...)
    } else {
      plot_offset(offset.obj=g1.offset,draw.legend=T,
                  main=paste0(Gene," | Before normalization"),...)
    }
  }

  ## loop until no different variations
  if (sum((g1.offset$g-1)^2)<1e-10) {
    cat("No adjustment applied. Data do not have enough expression.","\n")
    cat("........Plots are omitted.","\n")
    g2.offset=g1.offset
  } else {
    g2.offset = g1.offset
    k = 1
    while ((k<5) & (max(g2.offset$g)>1.05)) {
      g2.offset = estimate_offset(centerDataResult=NULL,
                                  msf=msf,centeredData=normalizedData,
                                  pileupData=pileupData,
                                  exonset=exonset,
                                  smoothness=smoothness,makePlot=FALSE)
      normalizedData = sweep(x=normalizedData,2,g2.offset$g,FUN="/")
      k = k+1
    }
    if (makePlot) {
      if (is.null(Gene)) {
        plot_offset(offset.obj=g2.offset,draw.legend=T,
                    main="After normalization",...)
      } else {
        plot_offset(offset.obj=g2.offset,draw.legend=T,
                    main=paste0(Gene," | After normalization"),...)
      }
    }

  }
  return(list(outputData=normalizedData,msf=msf,
              g1.offset=g1.offset,g2.offset=g2.offset,
              goodcase=goodcase))
}


#' @references https://github.com/hyochoi/SCISSOR/blob/master/R/center_data.R
#' @noRd

center_data <- function(inputData, Ranges, average="mean",trim=0.1,adjval=NULL, ...) {

  d = dim(inputData)[1]
  nexons = nrow(Ranges$lRanges)
  mean.vec = apply(inputData,1,FUN=function(x){adj.center(x,average=average,trim=trim,adjval=adjval)})

  exonbaseMat = matrix(0,ncol=nexons,nrow=d)
  colnames(exonbaseMat) = rownames(Ranges$lRanges)
  zeromean = c()
  for (ie in 1:nexons) {
    x = mean.vec[Ranges$lRanges[ie,2]:Ranges$lRanges[ie,3]]
    if (sum(x^2)>1e-05) {
      exonbaseMat[Ranges$lRanges[ie,2]:Ranges$lRanges[ie,3],ie] = x/sum(x^2)
    } else {
      zeromean = c(zeromean,ie)
    }
  }
  msf.exons = t(exonbaseMat)%*%inputData
  keepexons = which(! c(1:nrow(msf.exons)) %in% zeromean)
  if (length(keepexons)>1) {
    msf = apply(msf.exons[which(! c(1:nrow(msf.exons)) %in% zeromean),],2,median)
  } else if (length(keepexons)==1) {
    msf = msf.exons[which(! c(1:nrow(msf.exons)) %in% zeromean),]
  } else {
    msf = rep(0,dim(inputData)[2])
  }
  NewX = inputData - mean.vec%*%t(msf)
  # msf2=rep(1,length(msf))
  # msf2[which(msf>1)]=1/msf[which(msf>1)]
  # NewX = sweep(NewX,2,msf2,"*")
  return(list(outputData=NewX,msf=msf,data.center=mean.vec))
}


#' @references https://github.com/hyochoi/SCISSOR/blob/master/R/estimate_offset.R
#' @noRd

estimate_offset <- function(centerDataResult=NULL,msf=NULL,centeredData=NULL,
                           pileupData,exonset,
                           smoothness=0.7,
                           makePlot=FALSE, ...) {

  if (is.null(centerDataResult)) {
    if ((is.null(msf)) | (is.null(centeredData))) {
      stop("either of msf or centered matrix should be specified.")
    }
  } else {
    msf = centerDataResult$msf;
    centeredData = centerDataResult$outputData;
  }

  n=ncol(centeredData);
  case.sse = apply(centeredData,2,FUN=bisquare.sse)
  exonbase = c()
  for (j in 1:nrow(exonset)) {
    exonbase = c(exonbase,c(exonset[j,2]:exonset[j,3]))
  }
  goodcase = which(apply(pileupData[exonbase,],2,mean)>5)

  ## estimate g by lowess
  if (length(goodcase)>floor(0.1*n)) {
    lowess1 = lowess(case.sse[goodcase] ~ msf[goodcase],f=smoothness)

    tau = lowess1[[2]][which.min((lowess1[[1]]-1)^2)[1]]
    if (tau<0) {
      tau=min(lowess1[[2]][which(lowess1[[2]]>0)])
    }
    ## Rescale case.sse
    case.sse.adj = sqrt(case.sse/tau)
    x = lowess1[[1]][which(lowess1[[2]]>0)]
    y = sqrt(lowess1[[2]][which(lowess1[[2]]>0)]/tau)

    lowess2 = as.list(NULL);
    lowess2$x = msf;
    lowess2$y = interpolate.hy(x=x,y=y,newdata=msf)

    g = lowess2$y;
    g[which(g<1)] = 1;

    if (makePlot) {
      colmat=rep("darkgrey",n)
      colmat[goodcase[order(msf[goodcase],decreasing=TRUE)]] = rainbow(n=length(goodcase),start=0,end=0.8)

      yaxis = yaxis.hy(case.sse.adj)
      plot(msf,case.sse.adj,col=colmat,
           xlab="median scale factor",ylab="variability",
           ylim=c(yaxis[1],yaxis[2]+0.3*(yaxis[2]-yaxis[1])),...)
      lines(lowess2$x[order(lowess2$x)],lowess2$y[order(lowess2$x)],lwd=3,col="black")
      lines(msf[order(msf)],g[order(msf)],lwd=3,col="red")
      abline(v=1,h=1,lty=2)
      # points(msf[order(msf)],g[order(msf)],type="l",lwd=3,col="red")

      legend("topleft",bty="n",
             legend=c(paste("# cases involved =",length(goodcase)),
                      paste("gene length =",length(exonbase)),
                      paste("smoothness  =",round(smoothness,digits=1))))
    }
  } else {
    # No adjustment applied
    g = rep(1,n); case.sse.adj=case.sse;
    lowess2=NULL;

    if (makePlot) {
      colmat=rep("darkgrey",n)
      colmat[goodcase[order(msf[goodcase],decreasing=TRUE)]] = rainbow(n=length(goodcase),start=0,end=0.8)

      yaxis = yaxis.hy(case.sse.adj)
      plot(msf,case.sse.adj,col=colmat,
           xlab="median scale factor",ylab="variability",
           ylim=c(yaxis[1],yaxis[2]+0.3*(yaxis[2]-yaxis[1])), ...)
      abline(v=1,h=1,lty=2)
      # points(msf[order(msf)],g[order(msf)],type="l",lwd=3,col="red")

      legend("topleft",bty="n",
             legend=c(paste("No adjustment applied"),
                      paste("# cases involved =",length(goodcase)),
                      paste("gene length =",length(exonbase))))
    }

  }

  return(list(g=g,msf=msf,sse=case.sse.adj,goodcase=goodcase,
              lowess.fit=lowess2,smoothness=smoothness))
}


#' @references https://github.com/hyochoi/SCISSOR/blob/master/R/interpolate.hy.R
#' @noRd

interpolate.hy <- function(x,y,newdata) {

  inner.fn = function(x,y,newx) {
    x.sorted = x[order(x)]
    y.sorted = y[order(x)]

    # find two values that are closest to newdata and cover newdata between them
    x.diff = x.sorted-newx
    if (max(x.diff)<0) {
      # this case is when newx is beyond the maximum x
      newy = y.sorted[length(y.sorted)]
    } else if (min(x.diff)>0) {
      # this case is when newx is beyond the minimum x
      newy = y.sorted[1]
    } else {
      i = which.min(x.diff^2)
      if (x.diff[i]==0) {
        newy = y.sorted[i]
      } else {
        if (x.diff[i]<0) {
          il = i; iu = i+1;
        } else {
          il = i-1; iu = i;
        }
        if (abs(x.sorted[iu]-x.sorted[il])>1e-5) {
          b = (y.sorted[iu]-y.sorted[il])/(x.sorted[iu]-x.sorted[il])
          newy = y.sorted[il]+b*(newx-x.sorted[il])
        } else {
          newy = y.sorted[il]
        }
      }
    }
    return(newy)
  }
  return(apply(matrix(newdata,ncol=1),1,FUN=function(z){inner.fn(x=x,y=y,newx=z)}))
}


#' @references https://github.com/hyochoi/SCISSOR/blob/master/R/yaxis.hy.R
#' @noRd

yaxis.hy <- function(mat){
  #  mat : d by n matrix
  tempmax <- max(mat) ;
  tempmin <- min(mat) ;
  templen <- tempmax-tempmin ;
  return(c(tempmin-0.002*templen, tempmax+0.002*templen)) ;
}


#' @references https://github.com/hyochoi/SCISSOR/blob/master/R/plot_offset.R
#' @noRd

plot_offset <- function(offset.obj,draw.legend=TRUE,colmat=NULL,...) {

  g = offset.obj$g; msf = offset.obj$msf; n = length(g);
  goodcase = offset.obj$goodcase;
  sse = offset.obj$sse; lowess.fit = offset.obj$lowess.fit; smoothness=offset.obj$smoothness;

  if (is.null(colmat)) {
    colmat=rep("darkgrey",n)
    colmat[goodcase[order(msf[goodcase],decreasing=TRUE)]] = rainbow(n=length(goodcase),start=0,end=0.8)
  } else {
    colmat[-goodcase] = "darkgrey"
  }

  yaxis = yaxis.hy(sse)
  if (draw.legend) {
    plot(msf,sse,col=colmat,ylim=c(yaxis[1],yaxis[2]+0.2*(yaxis[2]-yaxis[1])),
         xlab="median scale factor",ylab="variation level",...)
  } else {
    plot(msf,sse,col=colmat,ylim=c(yaxis),
         xlab="median scale factor",ylab="variation level",...)
  }
  lines(lowess.fit$x[order(lowess.fit$x)],lowess.fit$y[order(lowess.fit$x)],lwd=2,col="black")
  lines(msf[order(msf)],g[order(msf)],lwd=3,col="red")
  abline(v=1,h=1,lty=2)
  # points(msf[order(msf)],g[order(msf)],type="l",lwd=3,col="red")

  if (draw.legend) {
    legend("topleft",bty="n",
           legend=c(paste("# cases involved =",length(goodcase)),
                    paste("smoothness  =",round(smoothness,digits=1))))
  }

}


#' @references https://github.com/hyochoi/SCISSOR/blob/master/R/decay.rate.hy.R
#' @noRd

decay.rate.hy <- function(Data,dist2end=NULL,robust=FALSE){
  # Calculate slopes & decay rates for every data point
  # Data = after applying center_data
  # dist2end = distance from the 3' end
  #            if null, use full transcript
  if (!is.null(dist2end)) {
    d = dim(data)[1];
    dist2end = min(d,dist2end);
    Data = Data[c((d-dist2end+1):d),];
    rm(d)
  }
  d = dim(Data)[1]; n = dim(Data)[2];
  lmcoef = matrix(0,ncol=2,nrow=n);
  x = 1:d ;
  for (case in 1:n){
    y = Data[,case];
    lmres = lm(y~x);
    lmcoef[case,] = coef(lmres);
  }
  rate = lmcoef[,2]*d
  return(list(slope=lmcoef[,2],rate=rate));
}


#' @references https://github.com/hyochoi/SCISSOR/blob/master/R/pd.rate.hy.R
#' @noRd

pd.rate.hy <- function(x,qrsc=FALSE) {
  # projection depth
  m = mad(x)
  if (m<1e-5) {
    return(rep(0,length(x)))
  } else {
    if (qrsc) {
      rsc=compScales(x)
      y=rep(0,length(x))

      above.ind=which((x-rsc$med)>=0)
      below.ind=which((x-rsc$med)<0)
      if (rsc$sa>1e-5) {
        y[above.ind]=(x[above.ind]-rsc$med)/rsc$sa
      }
      if (rsc$sb>1e-5) {
        y[below.ind]=(x[below.ind]-rsc$med)/rsc$sb
      }
      return(y)
    } else {
      return((x-median(x))/m)
    }
  }
}


#' @references https://github.com/hyochoi/SCISSOR/blob/master/R/compScales.R
#' @noRd

compScales <- function(x,
                      rmZeroes=FALSE,maxRatio=NULL,precScale=1e-10){
  # Computes the scales sa and sb (above and below the median).
  # Assumes that x is an array of numbers.
  #
  x = x[!is.na(x)] # we always take out NAs
  temp = fastSplitSample(x)
  xa   = temp$xa
  xb   = temp$xb
  med  = temp$med
  sall = scale1StepM((x-med),precScale=precScale)
  if(rmZeroes){ # reduces breakdown value but yields fewer implosions
    xa = xa[xa > precScale]
    xb = xb[xb > precScale]
  }
  sa = scale1StepM(xa,precScale=precScale)
  sb = scale1StepM(xb,precScale=precScale)
  if(!is.null(maxRatio)){
    if(maxRatio < 2) stop("maxRatio must be at least 2")
    sa = min(c(max(sa,sall/maxRatio,na.rm = TRUE),sall*maxRatio),
             na.rm = TRUE)
    sb = min(c(max(sb,sall/maxRatio,na.rm=TRUE),sall*maxRatio),
             na.rm = TRUE)
  }
  return(list(sa=sa,sb=sb,med=med))
}


#' @references https://github.com/hyochoi/SCISSOR/blob/master/R/fastSplitSample.R
#' @noRd

fastSplitSample <- function(x){
  # Centers sample by median, and divides in 2 equal halves.
  # Assumes that NAs have already been removed.
  # This function has time complexity O(n).
  #
  med = median(x) # takes only O(n) time
  x = x - med # centering
  n = length(x)
  h = n %/% 2   # = integer part of n/2
  xa = x[x > 0] # length(xa) <= h
  xb = x[x < 0] # length(xa) <= h
  xa = c(rep(0,(n - h - length(xa))),xa)
  xb = c(rep(0,(n - h - length(xb))),abs(xb)) # abs() !
  return(list(xa=xa,xb=xb,med=med))
}


#' @references https://github.com/hyochoi/SCISSOR/blob/master/R/scale1StepM.R
#' @noRd

scale1StepM <- function(x,precScale=1e-10) {
  # Computes the first step of an algorithm for
  # a scale M-estimator using the given rho function.
  # The scatter is computed relative to zero.
  #
  x = x[!is.na(x)] # we always take out NAs
  n = length(x)
  if(n == 0) { return(0.0)
  } else {
    sigma0 = 1.4826*median(abs(x))
    if(sigma0 < precScale) { return(0.0)
    } else {
      rho = rhoHuber(x/sigma0)
      return(sigma0 * sqrt(sum(rho)*2/n))
    }
  }
}


#' @references https://github.com/hyochoi/SCISSOR/blob/master/R/rhoHuber.R
#' @noRd

rhoHuber <- function(x,c=2.1){
  # x is a univariate sample
  # c is the tuning constant
  # output is rho(x)
  #
  rho = (x/c)^2
  rho[rho > 1] = 1
  1.54^2*rho
}
