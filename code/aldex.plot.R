##' plots for an ALDEx3 result object
##'
##' Provides volcano, effect, MA plots from an ALDEx3 result object.
##'
##' This method plots combinations of adjusted p-values from `object$p.val.adj`,
##' posterior estimates, standard errors and log abundance values
##' averaged across Monte Carlo samples.
##' The result is returned as a single plot with significant points colored red
##' Note-calls the aldex `summary` function internally
##' 
##' @title Plot Method for ALDEx3 Objects
##' @param output An object of class \code{aldex}
##' @param plot type of plot (default='volcano')
##' @param threshold FDR significance threshold (default=0.05)
##' @param min.diff (default=0.5) only used for MA plot type
##' @param contrast the name of the comparison
##' @return the desired plot
##' @export
##' @author Greg Gloor

aldex.plot <-function(object, plot=c("volcano", "effect", "MA"), threshold=0.05,
  min.diff=0.5, contrast=NULL){
  	# this allows partial matching for the lazy and defaults to volcano
    plot=match.arg(plot)

    # call to summary() 
    sum.output <- summary(object)
    nsamples <- length(object$data[,1])
    # get the sig features
    sig <- sum.output$p.val.adj < 0.05
    if(plot=="volcano"){
      # replace 0 with min pval/10
      p.val1 <- sum.output$p.val.adj > 0
      p.val0 <- sum.output$p.val.adj == 0
      min.p <- min(sum.output$p.val.adj[p.val1])
      sum.output$p.val.adj[p.val0] <- min.p/10
      y.val <- -(log10(sum.output$p.val.adj))
      plot(sum.output$estimate, y.val, pch=19, cex=0.5, col=rgb(0,0,0,0.3),
        xlab="estimate", ylab="-log10(p.adjust)")
      abline(h=-log10(threshold), lty=2)
      mtext(levels(object$data[[contrast]])[1], side=1, line = 2, at = min(sum.output$estimate),
            col = "grey", cex = 0.8)
      mtext(levels(object$data[[contrast]])[2], side=1, line = 2, at = max(sum.output$estimate),
            col = "grey", cex = 0.8)
      points(sum.output$estimate[sig], y.val[sig], pch=19, cex=0.5, col=rgb(1,0,0,0.5))
    }else if(plot=="effect"){
      plot(sum.output$std.error*sqrt(nsamples), sum.output$estimate,
      	pch=19, cex=0.5, col=rgb(0,0,0,0.3),xlab="std dev", ylab='estimate')
      mtext(levels(object$data[[contrast]])[1], side=2, line = 2, at = min(sum.output$estimate),
            col = "grey", cex = 0.8)
      mtext(levels(object$data[[contrast]])[2], side=2, line = 2, at = max(sum.output$estimate),
            col = "grey", cex = 0.8)
    
      points(sum.output$std.error[sig]*sqrt(nsamples), sum.output$estimate[sig],
      	pch=19, cex=0.5, col=rgb(1,0,0,0.5))
	  abline(0,1, lty=2, col='grey') 
	  abline(0,-1, lty=2, col='grey')
    }else if(plot=="MA"){
    ####
    # currently logComp and logScale can be too large to print out
    # in this case there is an error and no plot
    # example, yeast dataset with nsample > 128(ish)
    # there is a known hack in aldex()
    ####
      vals <- vector()
      for(i in 1:length(object$logComp[,1,1])){vals[i]= mean(object$logComp[i,,])}
      plot(vals,sum.output$estimate, pch=19, cex=0.5, col=rgb(0,0,0,0.3),
        xlab='log abundance', ylab='estimate')
      mtext(levels(object$data[[contrast]])[1], side=2, line = 2, at = min(sum.output$estimate),
            col = "grey", cex = 0.8)
      mtext(levels(object$data[[contrast]])[2], side=2, line = 2, at = max(sum.output$estimate),
            col = "grey", cex = 0.8)
      points(vals[sig], sum.output$estimate[sig], pch=19, cex=0.5, col=rgb(1,0,0,0.5))
      abline(h=min.diff, lty=2, col='grey')
      abline(h=-min.diff, lty=2, col='grey')
    }
  
}