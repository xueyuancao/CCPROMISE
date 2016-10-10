PROMISE2<- function (exprSet,            
                     exprSet2,           
                     geneSet = NULL,     
                     promise.pattern,    
                     strat.var = NULL,   
                     nbperm=FALSE,       
                     max.ntail=100,
                     nperms = 10000,      
                     seed = 13) 
{
    ph.pattern <- promise.pattern
    endpt.data <- pData(phenoData(exprSet))
    array.data <- exprs(exprSet)
    array.data2 <- exprs(exprSet2)
    if (is.data.frame(array.data)) 
        array.data <- as.matrix(array.data)
    if (!is.numeric(array.data)) {
        print("ERROR: The array data is not numeric")
        return()
    }
    array.names <- dimnames(array.data)[[2]]
    endpt.names <- dimnames(endpt.data)[[1]]
    array.inc <- is.element(array.names, endpt.names)
    endpt.inc <- is.element(endpt.names, array.names)
    nsamp.inc <- sum(array.inc)
    message(paste("No. of Matching array ids in endpt data and array data:", 
        nsamp.inc))
    if (nsamp.inc <= 4) {
        message(paste("Array IDs in endpt data:"))
        message(endpt.names)
        message("Array IDs in array data:")
        message(array.names)
        stop("Error: Unable to match array identifiers in endpt data and array data.")
    }
    Y <- array.data[, array.inc]
    Y2 <- array.data2[, array.inc]
    X <- endpt.data[endpt.inc, ]
    array.names <- array.names[array.inc]
    endpt.names <- endpt.names[endpt.inc]
    array.ord <- order(array.names)
    endpt.ord <- order(endpt.names)
    Y <- Y[, array.ord]
    Y2 <- Y2[, array.ord]
    X <- X[endpt.ord, ]
    Y <- as.matrix(Y)
    Y2 <- as.matrix(Y2)
    
    array.names <- array.names[array.ord]
    endpt.names <- endpt.names[endpt.ord]
    if (any(array.names != endpt.names)) {
        name.mtx <- cbind(endpt.arrayid = endpt.names, array.arrayid = array.names)
        message("----------------- IDs ------------------")
        message(name.mtx)
        stop("Error: Unable to align IDs in array data and endpt data")
    }
    endpt.endptcol <- is.element(dimnames(X)[[2]], 
        unique(unlist(strsplit(c(as.character(ph.pattern$endpt.vars), 
        strat.var), ","))))
    phtype <- X[, endpt.endptcol]
    names(phtype) <- dimnames(endpt.data)[[2]][endpt.endptcol]
    if (!is.numeric(phtype[, dimnames(phtype)[[2]] != strat.var]) & 
        !is.data.frame(phtype[, dimnames(phtype)[[2]] != strat.var])) {
        print("ERROR: The clinical data is not numeric or data frame")
        return()
    }
    if (is.vector(phtype)) 
        phtype <- matrix(phtype, length(phtype), 1)
    probes <- dimnames(array.data)[[1]]
    nprobes <- length(probes)
    if (!is.null(geneSet)) {
       print("ERROR: geneSet should be null")
        return()
    }
    
    if(is.null(geneSet)) {
        message(paste("Computing Stats for Observed Data:", date()))
        message(paste("Computing gene-specific statistics:", 
            date()))
        if (is.null(strat.var)) 
            strat <- NULL
        else strat <- phtype[, strat.var]
        gres1 <- PROMISE::promise.genestat(Y, phtype, ph.pattern, strat)
        gres2 <- PROMISE::promise.genestat(Y2, phtype, ph.pattern, strat)
        dimnames(gres1)[[2]]<-paste(dimnames(gres1)[[2]], "_1", sep="")
        dimnames(gres2)[[2]]<-paste(dimnames(gres2)[[2]], "_2", sep="")
        gres12<-cbind(gres1, gres2)
        prstat12<-cbind(gres1[, ncol(gres1)], gres2[, ncol(gres2)])
        gene.resf<-cbind(gres12, 
                    PROMISE.stat=apply(prstat12, 1, function(prb){
                                   return(mean(prb, na.rm=TRUE))}))
        
        m.nph <- dim(gene.resf)
        gene.pvals <- matrix(0, m.nph[1], m.nph[2])
        message(paste("Set seed: seed = ", seed))
        set.seed(seed)
        if(!nbperm){
           perms <- matrix(NA, nperms, nsamp.inc)
           if (is.null(strat.var)) 
               for (i in 1:nperms) perms[i, ] <- sample(nsamp.inc, replace = FALSE)
           else {
               strat <- endpt.data[, strat.var]
               ustrat <- sort(unique(strat))
               nstrat <- length(ustrat)
               ind <- length(strat)
               for (j in 1:nstrat) {
                   this.strat <- (strat == ustrat[j])
                   for (i in 1:nperms) perms[i, this.strat] <- 
                                       sample((1:ind)[this.strat], replace = FALSE)
               }
           }
           message(paste("Computing Stats for Permuted Data: ", date()))
           for (i in 1:nperms) {
               perm.ind <- unlist(perms[i, ])
               perm.phtype <- phtype[perm.ind, ]
               if (is.null(strat.var)) 
                   perm.strat <- NULL
               else perm.strat <- perm.phtype[, strat.var]
               gtemp1 <- PROMISE::promise.genestat(Y, perm.phtype, ph.pattern, perm.strat)
               gtemp2 <- PROMISE::promise.genestat(Y2, perm.phtype, ph.pattern, perm.strat)
                              dimnames(gtemp1)[[2]]<-paste(dimnames(gtemp1)[[2]], "_1", sep="")
               dimnames(gtemp2)[[2]]<-paste(dimnames(gtemp2)[[2]], "_2", sep="")
               gtemp12<-cbind(gtemp1, gtemp2)
               prtemp12<-cbind(gtemp1[, ncol(gtemp1)], gtemp2[, ncol(gtemp2)])
               gene.tempf<-cbind(gtemp12, 
                        PROMISE.stat=apply(prtemp12,1, function(prb){
                                     return(mean(prb, na.rm=TRUE))}))
                                                   
               gene.pvals <- gene.pvals + (abs(gene.tempf) >= abs(gene.resf))
           }
           message(paste("Generating Result List: ", date()))
           gene.pvals <- gene.pvals/nperms
           generes.tab <- cbind.data.frame(probeid = probes, gene.resf, perm.p = gene.pvals)
           message(paste("Finished at: ", date()))
           res <- list(generes = generes.tab, setres = NULL)
       }
       if (nbperm){
           #set up tracking
           probe.keep <- rep(TRUE, nprobes)
           probe.done <- rep(FALSE, nprobes)
           probe.ntail <- matrix(0, m.nph[1], m.nph[2])
           probe.apt.ntail <- rep(0, m.nph[1])
           num.perms <- rep(NA, nprobes)
           i <- 0
           while ((i < nperms) && (any(probe.keep))) {
               i <- i + 1
               if (is.null(strat.var)) 
                  perm.ind <- sample(nsamp.inc, replace = FALSE)
               else {
                  perm.ind<-1:nsamp.inc 
                  strat <- X[, strat.var]
                  ustrat <- sort(unique(strat))
                  nstrat <- length(ustrat)
                  ind <- length(strat)
                  for (j in 1:nstrat) {
                    this.strat <- (strat == ustrat[j])
                    perm.ind[this.strat] <- sample((1:ind)[this.strat], replace = FALSE)
                  }
                }
                perm.phtype <- phtype[perm.ind, ]
                if (is.null(strat.var)) 
                  perm.strat <- NULL
                else perm.strat <- perm.phtype[, strat.var]
                gtemp1 <- PROMISE::promise.genestat(Y[probe.keep, ], 
                              perm.phtype, ph.pattern, perm.strat)
                gtemp2 <- PROMISE::promise.genestat(Y2[probe.keep, ], 
                              perm.phtype, ph.pattern, perm.strat)
                dimnames(gtemp1)[[2]]<-paste(dimnames(gtemp1)[[2]], "_1", sep="")
                dimnames(gtemp2)[[2]]<-paste(dimnames(gtemp2)[[2]], "_2", sep="")
                gtemp12<-cbind(gtemp1, gtemp2)
                prtemp12<- cbind(gtemp1[, ncol(gtemp1)],gtemp2[, ncol(gtemp2)])
                gene.tempf<-cbind(gtemp12, 
                        PROMISE.stat=apply(prtemp12,  1, function(prb){
                                           return(mean(prb, na.rm=TRUE))}))
 
                probe.ntail[probe.keep, ] <- probe.ntail[probe.keep, ] + 
                           (abs(gene.tempf) >= abs(gene.resf[probe.keep, ]))
                probe.apt.ntail[probe.keep] <- probe.ntail[probe.keep, m.nph[2]]
                probe.done[probe.keep] <- probe.apt.ntail[probe.keep] ==max.ntail
                num.perms[probe.done & probe.keep] <- i
                gene.pvals[probe.done & probe.keep, ] <- 
                      probe.ntail[probe.done & probe.keep, ]/i
                probe.keep[probe.done] <- FALSE
            }
            if (any(probe.keep)) {
                gene.pvals[probe.keep, ] <- probe.ntail[probe.keep,]/nperms
                num.perms[probe.keep] <- nperms
            }
            dimnames(gene.pvals)[[2]] <- dimnames(gene.resf)[[2]]
            generes.tab <- cbind.data.frame(probeid = probes, 
                gene.resf, perm.p = gene.pvals, nperms = num.perms)
            message(paste("Finished at: ", date()))
            res <- list(generes = generes.tab, setres = NULL)
     }     
    }
    return(res)
}