# *--------------------------------------------------------------------
# | FUNCTION: visGAPath
# | Function for visualising the path of a genetic algorithmn using 
# | principal components analysis
# *--------------------------------------------------------------------
# | Version |Date      |Programmer  |Details of Change                
# |     01  |18/04/2012|Simon Raper |first version.      
# *--------------------------------------------------------------------
# | INPUTS:  func        The function to be optimised
# |          npar        The number of parameters to optimise over
# |          nruns       The number of runs of the GA to visualise
# |          parallel    If true runs on multiple cores
# |          file.path   File name and path for saving GA output
# |          domains     If needed for retricted optimisation
# |                      (see rgenoud documentation for more info)
# |          max         If true searches for maxima. Default=TRUE
# |          mut.weight  Option to upweight the mutation operator
# |                      (see rgenoud documentation for more info)
# |
# *--------------------------------------------------------------------
# | OUTPUTS: A ggplot object
# *--------------------------------------------------------------------
# | USAGE:   visGAPath(func, npar, nruns, parallel, file.path, 
# |          domains, max, mutation.weight)
# |
# *--------------------------------------------------------------------
# | DEPENDS: rgenoud, plyr, foreach, doSNOW, ggplot2
# |
# *--------------------------------------------------------------------
# | NOTES:   You may need to customise the multicore part of the code
# |          to suit your machine
# |
# *--------------------------------------------------------------------
 
 
visGAPath<-function(func, npar, nruns, parallel=TRUE, file.path='C:\\pro', domains=NULL, max=TRUE, mut.weight=50){
   
  #Function for creating data set of GA populations
   
  createGenDS<-function(gen.obj, path, run.name){
    gens<-gen.obj$generations
    pop<-gen.obj$popsize
    var.num<-length(gen.obj$par)
    project.in<-readLines(path)
    z<-2
    q<-1:2
    for (i in 1:gens){
      z<-z+pop
      t<-1:4+z
      q<-rbind(q,t(t))
      z<-z+4
    }
    numeric.only<-read.table(textConnection(project.in[-q]))
    generation<-rep(0:gens, each=pop)
    gen.df<-data.frame(generation, numeric.only)
    names(gen.df)<-c("Generation", "Ind", "Fit", paste("Var", 1:var.num, sep=""))
     
    bestInGen<-function(gen){
      best<-which.max(gen$Fit)
      best.ind<-gen[best,3:(var.num+3)]
      avg.ind<-colMeans(gen[,3:(var.num+3)])
      names(avg.ind)<-paste("avg.", names(avg.ind), sep="")
      ret<-data.frame(c(best.ind, avg.ind))
    }
     
    bestGen<-ddply(gen.df, .(Generation) ,bestInGen)
    run<-rep(run.name, (gens+1))
    bestGen<-data.frame(run, bestGen)
    bestGen
     
  }
   
  #ddply version
   
  wrap.gen<-function(arg){
    seed<-arg[1,2]
    run<-arg[1,1]
    pro.path<-paste(file.path, run, ".txt", sep="")
    gen <- genoud(func, nvars=npar,pop.size=500,max=max, project.path=pro.path, print.level=2, max.generations=100, unif.seed=seed, Domains=domains, P2=mut.weight)
    run<-createGenDS(gen, pro.path, run)
    run
  }
   
  runs<-data.frame(run=paste("run", 1:nruns, sep=" "), seed=floor(runif(nruns)*10000))
   
  if (parallel==FALSE){
     
    #Single core
    all<-ddply(runs, .(run), wrap.gen, .parallel = FALSE)
     
  } else {
    #Multi core
    cl <- makeCluster(getOption("cl.cores", min(8, nruns)))
    clusterExport(cl, objects(.GlobalEnv))
    clusterEvalQ(cl, library(plyr))
    clusterEvalQ(cl, library(mvtnorm))
    clusterEvalQ(cl, library(rgenoud))
    registerDoSNOW(cl)
    all<-ddply(runs, .(run), wrap.gen, .parallel = TRUE)
     
    stopCluster(cl)
     
  }
   
  #Visualise
   
  pc<-princomp(all[(5+npar):(5+npar*2-1)])
  scores<-data.frame(Run=factor(all$run), Fit=all$Fit, pc$scores)
  start<-scores[which(all$Generation==0),]
   
  g<-ggplot(data=scores, aes(x=Comp.1, y=Comp.2, colour=Run, group=Run))+geom_path()
  g+geom_point(data=start, size=5, aes(x=Comp.1, y=Comp.2)) + scale_x_continuous('Principal Component 1') + scale_y_continuous('Principal Component 2') 
}