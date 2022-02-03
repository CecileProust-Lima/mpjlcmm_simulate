## K = number of markers
## N = number of subjects
## ng = number of latent classes
## t = tmin1,tmax1,tecart1,trand1,...,tminK,tmaxK,tecartK,trandK
## Xbin = propX1,...,propXp1
## Xcont = meanXp1+1,sdXp1+1, ..., meanXp1+p2,sdXp1+p2
## beta_Y1 = list(I=c(I_G1,I_G2), t=c(t_G1,t_G2,etc), X1=c(X1_G1,X1_G2), tX1=c(tX1,tX1), etc)
## beta_Y2 = list(I=c(I_G1,I_G2), t=c(t_G1,t_G2,etc), X1=c(X1,X1), X2=c(X2,X2), etc)
## ... beta_YK
## Z = indicates the number of random effects for each component: no random effect (0), only intercept (1) or intercept and t (2)
## B1 = matrix of VC for random effects of marker Y1
## ... BK = matrix of VC for random effects of marker Y1
## linear1, ..., linearK = parameters for the transformations
## weib = parameters for Weibull risk of event (2 prm or 2*ng prm)
## beta_S = list(X1=c(X1_G1,X1_G2), X2=c(X2_G1,X2_G2,etc))
## age0 = for the time scale in age rather than delay _ parameters min et max 
## NB : if age0 : age=(ti+a-65)/10
## scalet : expression to go from the time scale for longi to the one for survival (ex: scalet=t*0.5/11)
## piecewise : list with $nodes, $brisq et $ph

simdata_mpjlcmm <- function(K,N,ng,t,age0=NULL,Xbin=NULL,Xcont=NULL,beta_Y1,beta_Y2=NULL,Z,B1,B2=NULL,linear1,linear2=NULL,sigma=c(1,1),prop,weib=NULL,beta_S=NULL,logscale=FALSE,rweib=TRUE,piecewise=NULL,scalet,seed)
{
    if(missing(seed)) seed <- round(abs(rnorm(1,mean=5)),5)*100000
    set.seed(seed)
    print(seed)
    if(is.vector(Xcont)) Xcont <- matrix(Xcont,length(Xcont)/2,2,byrow=TRUE)
    nx <- length(Xbin)+length(Xcont)/2
    
    argt <- t
    if(is.vector(t) & length(t)>4) t <- matrix(t,K,4,byrow=TRUE)

    beta1 <- matrix(unlist(beta_Y1),ncol=ng,byrow=TRUE)
    rownames(beta1) <- names(beta_Y1)
    beta2 <- NULL
    if(!is.null(beta_Y2))
    {
        beta2 <- matrix(unlist(beta_Y2),ncol=ng,byrow=TRUE)
        rownames(beta2) <- names(beta_Y2)
    }
    beta <- list(beta1,beta2)
    
    linear <- rbind(linear1,linear2)

    B <- list(B1,B2)

    if(!is.null(weib))
    {
        if(length(weib)==2) weib <- rep(weib,ng)
        if(logscale==FALSE) weib <- weib^2 #logscale=FALSE by defaut
        if(logscale==TRUE) weib <- exp(weib)
    }
    
    ## simul data for n subjects
    res <- NULL
    i <- 1
    nb <- 0
    while(i<=N & nb<10*N)
    {
        nb <- nb+1
        ## simuler age0
        if(length(age0))
        {
            if(length(age0)==2) a <- runif(1,age0[1],age0[2])
            if(length(age0)==1) a <- rnorm(1,Xcont[age0,1],Xcont[age0,2])
            if(!(length(age0) %in% 1:2)) stop("mauvais age0")
        }
        
        ## simuler les temps si ce sont les memes dans chaque composante
        tik <- NULL
        if(length(t)==4)
        {
            tik <- t[1]
            nbt <- 1
            while(tik[nbt]<t[2])
            {
                tik <- c(tik,t[3]*nbt+runif(1,-t[4],t[4]))
                nbt <- nbt+1
            }
            if(tik[nbt]>t[2]) tik <- tik[-nbt]
                        
            if(length(age0))
            {
                tik <- (tik+a-65)/10
            }           
            
            Qik <- tik*tik
        }

        
        ## simuler les variables explicatives binaires
        if(!is.null(Xbin))
        {
            varexpli <- sapply(Xbin,function(p) rbinom(1,size=1,prob=p))
        }
        
        ## simuler les variables continues
        if(!is.null(Xcont))
        {
            varexpli <- c(varexpli,apply(Xcont,1,function(prm) rnorm(1,mean=prm[1],sd=prm[2])))
            if(length(age0))
            {
                varexpli[length(Xbin)+age0] <- (a-65)/10
            }
        }
        
        names(varexpli) <- paste("X",1:length(varexpli),sep="")
        #if(!all(names(varexpli) %in% vars)) stop("Probleme dans les variables")

        ## fonctions pour trouver une variable
        funx <- function(x,vars){
            if(paste("X",x,sep="") %in% vars) return(which(vars==paste("X",x,sep="")))
            else {
                if(paste("x",x,sep="") %in% vars) return(which(vars==paste("x",x,sep="")))
                else return(0)
            }
        }
        
        funtx <- function(x,vars){
            if(paste("tX",x,sep="") %in% vars) return(which(vars==paste("tX",x,sep="")))
            else {
                if(paste("tx",x,sep="") %in% vars) return(which(vars==paste("tx",x,sep="")))
                else return(0)
            }
        }

        funt2x <- function(x,vars){
            if(paste("t2X",x,sep="") %in% vars) return(which(vars==paste("t2X",x,sep="")))
            else {
                if(paste("t2x",x,sep="") %in% vars) return(which(vars==paste("t2x",x,sep="")))
                else return(0)
            }
        }
        
        ## proba d'appartenance aux classes
        cx <- sapply(1:nx,funx,vars=names(prop))
        xprob <- c(1,rep(0,length(prop)-1))
        if(length(prop)>1)
        {
            for(j in 1:nx)
            {
                if(cx[j]!=0) xprob[1+cx[j]] <- varexpli[j]
            }
        }
        betap <- matrix(unlist(prop),ncol=ng,byrow=TRUE)
        rownames(betap) <- names(prop)
        proba <- as.vector(exp(t(xprob)%*%betap))
        proba <- proba/sum(proba)
        pcl <- runif(1,0,1)
        classe <- as.numeric(cut(pcl,breaks=c(0,cumsum(proba))))
                                        #if(i<=10) cat("i=",i," proba=",proba," pcl=",pcl, " classe=",classe ,"betap=",betap,"xprob=",xprob,"\n")

        

        ## simuler le temps de survie
        Tsurv <- NA
        Dsurv <- NA
        if(!is.null(weib) | !is.null(piecewise))
        {
            if(!is.null(beta_S))
            {
                cx <- sapply(1:nx,funx,vars=names(beta_S))
                xsurv <- rep(0,length(beta_S))
                if(length(beta_S)>0)
                {
                    for(j in 1:nx)
                    {
                        if(cx[j]!=0) xsurv[cx[j]] <- varexpli[j]
                    }
                }
                bsurv <- matrix(unlist(beta_S),ncol=ng,byrow=TRUE)
                names(xsurv) <- names(beta_S)
            }
            else
            {
                xsurv <- 0
                bsurv <- matrix(0,ncol=ng)
            }
            pevt <- runif(1,0,1)

            if(!is.null(piecewise))
            {
                zl <- piecewise$nodes
                bl <- piecewise$brisq
                ph <- c(piecewise$ph,0)

                if(logscale==TRUE) bl <- exp(bl)
                else bl <- bl^2
 
                surv <- function(t,zl,bl,ph=1,xb,p)
                {
                    j <- which.max(zl[which(zl<=t)])
                    if(j==1) som <- 0 
                    else som <- sum(bl[1:(j-1)]*(zl[2:j]-zl[1:(j-1)]))
                    
                    if(j<length(zl)) surv <- exp(-(som+bl[j]*(t-zl[j]))*exp(ph)*exp(xb))
                    else surv <- exp(-som*exp(ph)*exp(xb))
                    
                    return(surv-p)
                }


                zero <- try(uniroot(surv,interval=c(zl[1],zl[length(zl)]),
                                    zl=zl,bl=bl,ph=ph[classe],
                                    xb=t(xsurv)%*%bsurv[,classe],p=pevt),
                            silent=TRUE)
                if(class(zero)=="try-error") Tevt <- Inf
                else Tevt <- zero$root

            }
            else
            {
                if(rweib!=TRUE)
                {
                    if(logscale==TRUE)
                    {
                        Tevt <- (-log(pevt)/(weib[2*(classe-1)+1]*exp(t(xsurv)%*%bsurv[,classe])))^(1/weib[2*(classe-1)+2]) # pour logscale=TRUE
                    }
                    else
                    {
                        Tevt <- 1/(weib[2*(classe-1)+1])*(-log(pevt)/(exp(t(xsurv)%*%bsurv[,classe])))^(1/weib[2*(classe-1)+2])
                    }
                }
                else
                {
                    if(logscale==TRUE)
                    {
                        Tevt <- rweibull(1,weib[2*(classe-1)+2],(weib[2*(classe-1)+1]*exp(-t(xsurv)%*%bsurv[,classe]))^(1/weib[2*(classe-1)+2]))
                                        # shape=b scale=(a*exp(-xbeta))^(1/b)
                    }
                    else
                    {
                        Tevt <- rweibull(1,weib[2*(classe-1)+2],exp(-t(xsurv)%*%bsurv[,classe])^(1/weib[2*(classe-1)+2])*1/weib[2*(classe-1)+1])
                                        # shape=b scale=(exp(-xbeta)^(1/b))/a
                    }
                }
            }
#if(i<10) cat("i=",i," a=",weib[2*(classe-1)+1]," b=",weib[2*(classe-1)+2]," tevt=",Tevt,"classe=",classe,"xsurv=",xsurv,"bsirv=",bsurv,"pevt=",pevt,"\n")

            if(!is.null(scalet))
            {
                tmin_scale <- eval(parse(text=gsub("t","t[1]",scalet)))
                tmax_scale <- eval(parse(text=gsub("t","t[2]",scalet)))
            }
            else
            {
                tmin_scale <- t[1]
                tmax_scale <- t[2]
            }
            
            Tcens <- runif(1,tmin_scale,2*tmax_scale)
            if(Tcens>tmax_scale) Tcens <- tmax_scale

            Tsurv <- min(Tevt,Tcens)
            Dsurv <- ifelse(Tsurv==Tevt,1,0)
        }
        enlevi <- FALSE
#if(i<10) cat("i=",i," a=",a," tevt=",Tevt, " Tcens=",Tcens, " Tsurv=",Tsurv,"\n")
          
        ## creer X et Z
        for(k in 1:K)
        {
            ## variables
            vars <- rownames(beta[[k]])
            cx <- sapply(1:nx,funx,vars=vars)
            ctx <- sapply(1:nx,funtx,vars=vars)
            ct2x <- sapply(1:nx,funt2x,vars=vars)
            ct <- which(vars=="t")
            ct2 <- which(vars=="t2")

            quad <- FALSE
            if("t2" %in% vars) quad <- TRUE
            
            ## simuler les temps
            if(length(t)>4)
            {
                tik <- t[k,1]
                nbt <- 1
                while(tik[nbt]<t[k,2])
                {
                    tik <- c(tik,t[k,3]*nbt+runif(1,-t[k,4],t[k,4]))
                    nbt <- nbt+1
                }
                if(tik[nbt]>t[2]) tik <- tik[-nbt]
                
                if(length(age0))
                {
                    tik <- (tik+a-65)/10
                }
                
                Qik <- tik*tik
            }

            ## arreter les mesures à Tsurv
            if(!is.na(Tsurv))
            {
                if(!is.null(scalet))
                {
                    tik_scale <- eval(parse(text=gsub("t","tik",scalet)))
                }
                else
                {
                    tik_scale <- tik
                }
                if(length(age0))
                {
                    tik <- tik[which(tik_scale<=(Tsurv+a-65)/10)]
                }
                else
                {
                    tik <- tik[which(tik_scale<=Tsurv)]
                }
                Qik <- tik*tik

                if(length(tik)==0)
                {
                    enlevi <- TRUE
                    break
                }
                
            }

            ## ajouter le numero de visite
            visite <- c(0:(length(tik)-1))

            ## creer xik
            xik <- matrix(0,length(tik),length(vars))
            xik[,1] <- 1
            xik[,ct] <- tik
            if(quad==TRUE)
            {
                xik[,ct2] <- Qik
            }
            for(j in 1:nx)
            {
                if(cx[j]!=0) xik[,cx[j]] <- varexpli[j]
                if(ctx[j]!=0) xik[,ctx[j]] <- tik*varexpli[j]
                if(ct2x[j]!=0) xik[,ct2x[j]] <- Qik*varexpli[j]
            }
            
            ## creer zik
            zik <- NULL
            if(Z[k]==1)
            {
                zik <- matrix(1,nrow=length(tik),ncol=1)
            }
            
            if(Z[k]==2)
            {
                zik <- cbind(1,tik)
            }
            
            if(Z[k]==3)
            {
                zik <- cbind(1,tik,Qik)
            }

            ## simuler les effets aleatoires
            if(!is.null(zik))
            {
                if(is.vector(B[[k]]))
                {
                    BB <- matrix(0,Z[k],Z[k])
                    BB[upper.tri(BB,diag=TRUE)] <- B[[k]]
                    BB <- t(BB)
                    BB[upper.tri(BB,diag=TRUE)] <- B[[k]]
                    B[[k]] <- BB
                }
                chB <- t(chol(B[[k]]))
                u01 <- rnorm(Z[k])
                ui <- chB %*% u01
            }
            else
            {
                zik <- matrix(0,nrow(xik),ncol(xik))
                ui <- rep(0,ncol(zik))
            }
            
            ## simuler les erreurs de mesure pour Y
            epsi <- rnorm(length(tik),sd=sigma[k])

            ## calculer Xbeta+Zu pour chaque classe et Y
            xbzu <- xik%*%beta[[k]][,classe] + zik%*%ui
            yi <- (xbzu + epsi)*linear[k,2] + linear[k,1]

            ## faire le data frame
            if(k==1)
            {
                datai <- data.frame(i=i, classe=classe,t=tik, visite, matrix(varexpli,nrow=length(tik),ncol=length(varexpli),byrow=TRUE), xbzu, zik%*%ui, epsi, yi, Tsurv, Dsurv)
                colnames(datai) <- c("i","classe","t","visite",names(varexpli),"xbzu1","zu1","eps1","Y1","Tsurv","Dsurv")
            }
            else
            {
                if(length(t)==4)
                {
                    datai <- data.frame(datai, xbzu, zik%*%ui, epsi, yi)
                    colnames(datai) <- c("i","classe","t","visite",names(varexpli),"xbzu1","zu1","eps1","Y1","Tsurv","Dsurv","xbzu2",paste(c("zu","eps","Y"),k,sep=""))
                }
                else
                {
                    datai <- data.frame(K=1,datai)
                    colnames(datai) <- c("K","i","classe","t","visite",names(varexpli),"xbzu","zu","eps","Y","Tsurv","Dsurv")
                    dataik <- data.frame(K=k, i=i, t=tik, visite=visite, matrix(varexpli,nrow=length(tik),ncol=length(varexpli),byrow=TRUE), xbzu, zik%*%ui, epsi, yi, Tsurv, Dsurv)
                    colnames(dataik) <- colnames(datai)
                    
                    datai <- rbind(datai,dataik)
                }
            }
                       
        } # fin boucle k

        if(enlevi==FALSE)
        {
            res <- rbind(res,datai)
            i <- i+1
        }
    } # fin while sujet

    return(res)
}


## pour i:
## affecter un g selon les proba p_ig
## simuler Yg (pour ce g défini)
## simuler survie en weibull :  p=exp(-aT^c*exp(xb))
## p ~ unif(0,1)
##  log(p)/exp(xb) = -aT^c =>  T = sqrt_c(-log(p)/(a*exp(xb)))
## Tc = unif()
## T* = min(T,Tc), di= indic(T*==T)

## en piecewise:

zl <- c(0,0.1533196, 0.3039014, 0.495551, 1.243 ) #hazardnodes
bl <- c( 0.2416715,0.8212857, 1.5068744,0.6197200) #brisq^2
 
surv1 <- function(t,zl,bl,ph=1,p)
{
    j <- which.max(zl[which(zl<=t)])
    if(j==1) som <- 0 
    else som <- sum(bl[1:(j-1)]*(zl[2:j]-zl[1:(j-1)]))

    if(j<length(zl)) surv <- exp(-(som+bl[j]*(t-zl[j]))*exp(ph))
    else surv <- exp(-som*exp(ph))

    return(surv-p)
}


uniroot(surv1,interval=c(0,1.243),zl=zl,bl=bl,p=0.7)
