"prim_perm_pca" <- function(structure1,structure2,struct.names="IMR90_HMEC",chrom="chr21",peel=0.1,paste=0.01,mass=0.001, pasting=F,n.perm=1000){
	
	##  structure1,structure2 == mutually aligned 3D reconstructions for a given chromosome:: genomic locus, x,y,z (non-missing coordinates)
	##  struct.names, chrom: text for labeling output graphics files
	##  peel, paste, mass, pasting :: PRIM parameters
	##  n.perm :: number of permutations
	
	require(prim)
 
	structure1 <- structure1[,2:4]
	structure2 <- structure2[,2:4]
		
	t.test.1 <- t.test.2 <- t.test.3 <- t.test.med <- as.data.frame(matrix(nrow=n.perm,ncol=3))

	names(t.test.1) <- names(t.test.2) <- names(t.test.3) <- names(t.test.med) <- 
	c("t-statistic","df","p-value")
	
	perm.y.fun <- as.data.frame(matrix(nrow=n.perm,ncol=4))
	names(perm.y.fun) <- c("box_1","box_2","box_3","box_med")
	
	euc.dist <- sqrt(rowSums((structure1 - structure2)^2))
	
	rot.structure1 <- prcomp(structure1)$x
		
	prim.orig <- prim.box(x=rot.structure1,y=euc.dist,peel.alpha=peel,paste.alpha=paste,mass.min=mass,threshold.type=1,pasting=F)
		
	y.fun.o <- prim.orig$y.fun
	n.boxes.o <- length(y.fun.o) - 1
	rev.sort.yfo <- rev(sort(y.fun.o))
					
	y.fun.1 <- rev.sort.yfo[1]		
	y.fun.2 <- rev.sort.yfo[2]
	y.fun.3 <- rev.sort.yfo[3]
	y.fun.med <- rev.sort.yfo[round(((1+n.boxes.o)/2),0)]
		
	ind.1 <- which(y.fun.o == rev.sort.yfo[1])[1]		## in case not unique
	ind.2 <- which(y.fun.o == rev.sort.yfo[2])[1]
	ind.3 <- which(y.fun.o == rev.sort.yfo[3])[1]
	ind.med <- which(y.fun.o == rev.sort.yfo[round(((1+n.boxes.o)/2),0)])[1]

				
	for(i in 1:n.perm){
		y.i <- sample(euc.dist)
		prim.i <- prim.box(x=rot.structure1,y=y.i,peel.alpha=peel,paste.alpha=paste,mass.min=mass,threshold.type=1,pasting=F)
		y.fun.i <- prim.i$y.fun
		n.boxes.i <- length(y.fun.i) - 1
			
		ind.1.i <- which(y.fun.i == rev(sort(y.fun.i))[1])[1]		## in case max not unique
		ind.2.i <- which(y.fun.i == rev(sort(y.fun.i))[2])[1]
		ind.3.i <- which(y.fun.i == rev(sort(y.fun.i))[3])[1]
		ind.med.i <- which(y.fun.i == rev(sort(y.fun.i))[round(((1+n.boxes.i)/2),0)])[1]
			
		perm.y.fun[i,1] <- max(y.fun.i)		
		perm.y.fun[i,2] <- rev(sort(y.fun.i))[2]
		perm.y.fun[i,3] <- rev(sort(y.fun.i))[3]
		perm.y.fun[i,4] <- rev(sort(y.fun.i))[round(((1+n.boxes.i)/2),0)]
						
		t.test.1.i <- t.test(prim.orig$y[[ind.1]],prim.i$y[[ind.1.i]])
		t.test.2.i <- t.test(prim.orig$y[[ind.2]],prim.i$y[[ind.2.i]])
		t.test.3.i <- t.test(prim.orig$y[[ind.3]],prim.i$y[[ind.3.i]])
		t.test.med.i <- t.test(prim.orig$y[[ind.med]],prim.i$y[[ind.med.i]])
						
		t.test.1[i,] <- c(t.test.1.i$statistic,t.test.1.i$parameter,t.test.1.i$p.value)
		t.test.2[i,] <- c(t.test.2.i$statistic,t.test.2.i$parameter,t.test.2.i$p.value)
		t.test.3[i,] <- c(t.test.3.i$statistic,t.test.3.i$parameter,t.test.3.i$p.value)
		t.test.med[i,] <- c(t.test.med.i$statistic,t.test.med.i$parameter,t.test.med.i$p.value)
		}			
		
	pdf(file=paste(paste(paste(paste("relocal-hists-structures",struct.names,sep="-"),"chrom",chrom,sep="-"),"n.perm",n.perm,sep="-"),"pdf",sep="."))
	
	par(mfrow=c(2,2))

	hist(na.omit(perm.y.fun[,1]),col=4,main="Box 1",xlab="Relocalization",xlim=c(min(perm.y.fun[,1],y.fun.1)-0.01,max(perm.y.fun[,1],y.fun.1)+0.01),breaks = 20)
	abline(v=y.fun.1,col=2,lty=1)
	
	hist(na.omit(perm.y.fun[,2]),col=4,main="Box 2",xlab="Relocalization",xlim=c(min(perm.y.fun[,2],y.fun.2)-0.01,max(perm.y.fun[,2],y.fun.2)+0.01),breaks = 20)
	abline(v=y.fun.2,col=2,lty=1)

	hist(na.omit(perm.y.fun[,3]),col=4,main="Box 3",xlab="Relocalization",xlim=c(min(perm.y.fun[,3],y.fun.1)-0.01,max(perm.y.fun[,3],y.fun.3)+0.01),breaks = 20)
	abline(v=y.fun.3,col=2,lty=1)

	hist(na.omit(perm.y.fun[,4]),col=4,main="Median Box",xlab="Relocalization",xlim=c(min(perm.y.fun[,4],y.fun.med)-0.01,max(perm.y.fun[,4],y.fun.med)+0.01),breaks = 20)
	abline(v=y.fun.med,col=2,lty=1)

	dev.off()	


	results <- list(prim.orig = prim.orig, y.fun.1 = y.fun.1, y.fun.2 = y.fun.2, y.fun.3 = y.fun.3, y.fun.med = y.fun.med, ind.1 = ind.1, ind.2 = ind.2, ind.3 = ind.3, ind.med = ind.med, perm.y.fun = perm.y.fun, t.test.1 = t.test.1, t.test.2 = t.test.2, t.test.3 = t.test.3, t.test.med = t.test.med)

	results 		
	}	
			
			
			
			