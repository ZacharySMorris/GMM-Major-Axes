### Iteratively subset GPA, covariates, and PC scores for each clade into own matrix ###
simple_subsetGMM <- function(X,PCData,group)
{
  # X is a data.frame of covariate data
  # PCData is a shape PCA dataset output from geomorph::gm.prcomp()
  # group is a character string that identifies the covariate to create subgroups
  
  #create lists for data to be subset
  PCList=list()
  
  group_classifier <- X[[group]]
  Groups <- levels(group_classifier)
  
  for (g in 1:length(Groups)){
    z<-NULL
    z<-as.matrix(PCData$x[grep(Groups[g],group_classifier),])
    if (ncol(z)==1) {
      z <- t(z)
      PCList[[Groups[g]]]<-z
    }else{
      PCList[[Groups[g]]]<-z
    }
  }
  out <- list(name = paste(group),
              groups = Groups,
              PCvalues = PCList
              )
  
  out
}
###

### Function for merging datasets created with simple_subsetGMM ###
mergeGMM <- function(X,Y,Xgroups = NULL,Ygroups = NULL){
  # X and Y are outputs from simple_subsetGMM
  # Xgroups and Ygroups are otpional character strings that identify subgroups to be merged
  
  if (is.null(Xgroups)){
    Xgroups = X$groups
  }
  if (is.null(Ygroups)){
    Ygroups = Y$groups
  }
  
  Zgroups <- c(Xgroups,Ygroups)
  ZPCvalues <- append(X$PCvalues[Xgroups],Y$PCvalues[Ygroups])
  Zname <- paste(X$name,Y$name,sep = "_vs_")
  
  out <- list(name = Zname,
              groups = Zgroups,
              PCvalues = ZPCvalues)
  out
}
###

### Function for merging created with simple_subsetGMM to create a single PCA dataset and classifier label ###
mergePCA <- function(X,Y,Xgroups = NULL,Ygroups = NULL){
  # X and Y are outputs from simple_subsetGMM
  # Xgroups and Ygroups are otpional character strings that identify subgroups to be merged
  
  if (is.null(Xgroups)){
    Xgroups = X$groups
  }
  if (is.null(Ygroups)){
    Ygroups = Y$groups
  }
  
  #check that all subsets have the same PCs
  if (!0 == var(c(sapply(X$PCvalues,FUN = ncol), sapply(Y$PCvalues,FUN = ncol)))){stop("The number of PCs do not match among subsets and cannot be coerced into single matrix.")
  }
  
  Xmatrix <- do.call(rbind,X$PCvalues[c(Xgroups)])
  Xclassifier <- rep(Xgroups, sapply(X$PCvalues[c(Xgroups)],FUN = nrow))
  Ymatrix <- do.call(rbind,Y$PCvalues[c(Ygroups)])
  Yclassifier <- rep(Ygroups, sapply(Y$PCvalues[c(Ygroups)],FUN = nrow))
  
  Zmatrix <- rbind(Xmatrix,Ymatrix)
  Zclassifier <- factor(c(Xclassifier,Yclassifier))
  
  out <- list(x=Zmatrix,classifier=Zclassifier)
  out
}
###

### Iteratively subset GPA, covariates, and PC scores for each clade into own matrix ###
SubsettingGMM_V2 <- function(X,A,PCData,group)
{
  # X is a data.frame of covariate data
  # A is a shape dataset output from gpagen
  # PCData is a shape PCA dataset output from geomorph::gm.prcomp()
  # group is a character string that identifies the covariate to create subgroups
  
  #create lists for data to be subset
  GPAList=list()
  CovariatesList=list()
  SpeciesList=list()
  PCList=list()
  CSList=list()
  Taxa=list()
  
  group_classifier <- X[[group]]
  Groups <- levels(group_classifier)
  
  #pb   <- txtProgressBar(1, length(levels(X[[group]])), style=3)
  for (g in 1:length(Groups)){
    x<-NULL
    x<-as.array(A$coords[,,grep(Groups[g],group_classifier)])
    y<-NULL
    y<-X[grep(Groups[g],group_classifier),]
    z<-NULL
    z<-as.matrix(PCData$x[grep(Groups[g],group_classifier),])
    {
      GPAList[[print(Groups[g])]]<-as.array(x)
      CSList[[Groups[g]]]<-A$Csize[grep(Groups[g],group_classifier)]
      CovariatesList[[Groups[g]]]<-y
      SpeciesList[[Groups[g]]]<-y$Species
      Taxa<-names(GPAList)
    }
    if (ncol(z)==1) {
      z <- t(z)
      PCList[[Groups[g]]]<-z
    }else{
      PCList[[Groups[g]]]<-z
    }
    #setTxtProgressBar(pb, i)
  }
  out <- list(name = paste(group), groups = Taxa,
              coords = GPAList, PCvalues = PCList, CSize = CSList,
              covariates = CovariatesList,
              species = SpeciesList
              )
  
  out
}
###


### Iteratively subset GPA, covariates, and PC scores for each clade into own matrix ###
SubsettingGMM <- function(X,A,PCData,W,group,print.plot=FALSE)
{

  # X is a data.frame of covariate data
  # A is a shape dataset output from gpagen
  # PCData is a shape PCA dataset output from geomorph::gm.prcomp()
  # W is a dataset of plotting values (output from PlottingValues)
  # group is a character string that identifies the covariate to create subgroups
  
  #create lists for data to be subset
  GPAList=list()
  CovariatesList=list()
  SpeciesList=list()
  PCList=list()
  CSList=list()
  Sizes=list()
  Colors=list()
  Shapes=list()
  Legendcolors=c()
  Legendshapes=c()
  Taxa=list()
  
  PC_VarList <- summary(PCData)$PC.summary[2,]

  #make x and y axis labels with % variance on them
  x_lab <- paste("Principal Component 1", " (", round(PC_VarList[1]*100, 1), "%)", sep="")
  y_lab <- paste("Principal Component 2", " (", round(PC_VarList[2]*100, 1), "%)", sep="")


  #Plot PCA alone#
  Xlim<-c(floor(min(PCData$x[,1])*10)/10,ceiling(max(PCData$x[,1])*10)/10)
  Ylim<-c(floor(min(PCData$x[,2])*10)/10,ceiling(max(PCData$x[,2])*10)/10)

  group_classifier <- X[[group]]
  Groups <- levels(group_classifier)

  #pb   <- txtProgressBar(1, length(levels(X[[group]])), style=3)
  for (g in 1:length(Groups)){
    x<-NULL
    x<-as.array(A$coords[,,grep(Groups[g],group_classifier)])
    y<-NULL
    y<-X[grep(Groups[g],group_classifier),]
    z<-NULL
    z<-as.matrix(PCData$x[grep(Groups[g],group_classifier),])
    {
      GPAList[[print(Groups[g])]]<-as.array(x)
      CSList[[Groups[g]]]<-A$Csize[grep(Groups[g],group_classifier)]
      CovariatesList[[Groups[g]]]<-y
      SpeciesList[[Groups[g]]]<-y$Species
      Shapes[[Groups[g]]]<-W$shape[grep(Groups[g],group_classifier)]
      Sizes[[Groups[g]]]<-W$size[grep(Groups[g],group_classifier)]
      Colors[[Groups[g]]] <- W$color[grep(Groups[g],group_classifier)]
      Taxa<-names(GPAList)
    }
    if (ncol(z)==1) {
      z <- t(z)
      PCList[[Groups[g]]]<-z
    }else{
      PCList[[Groups[g]]]<-z
    }
    #setTxtProgressBar(pb, i)

    if (print.plot == TRUE){
      pdf(file=paste("PC1vPC2_DorsalSkull_", Taxa[[g]], ".pdf", sep = ""), width = 11, height = 8.5, useDingbats = FALSE)

      plot(0, 0, type = "n",
           xlim = Xlim,
           ylim = Ylim,
           xlab = x_lab,
           ylab = y_lab,
           main = paste(Taxa[[g]], "from Combined Morphospace", sep = " "),
           axes = FALSE,
           frame.plot = FALSE,
           asp=F)

      axis(1, round(seq(Xlim[1],Xlim[2],by=0.1),1), pos=Ylim[1])
      axis(2, round(seq(Ylim[1],Ylim[2],by=0.1),1), pos=Xlim[1])
      clip(Xlim[1],Xlim[2],Ylim[1],Ylim[2])
      abline(h=0, lty=3)
      abline(v=0, lty=3)

      points(z[,1], z[,2],
             pch=Shapes[[Groups[g]]],
             cex=Sizes[[Groups[g]]],
             bg=alpha(Colors[[Groups[g]]], 0.75),
             xlim=Xlim, ylim = Ylim,
             xlab=x_lab, ylab = y_lab, asp=F)

      dev.off()
    }


    LegendColors <- Colors
    Legendshapes=list()

  }
  out <- list(name = paste(group), groups = Taxa,
              coords = GPAList, PCvalues = PCList, CSize = CSList,
              covariates = CovariatesList,
              species = SpeciesList,
              color = Colors, size = Sizes, shape = Shapes)

  out
}
###
