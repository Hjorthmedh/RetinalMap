#nb stepbreaks values should be strictly increasing.

run_WC <-
function(rgc.file,sc.file,op.file, k=1, 
	alpha=0.1,
	sijmin=0.00001,
	timestep=0.0001,
	dist=0.07,
	dist.bij=0.0289,
	stepbreaks=c(1,2500,5000,7500,10000,12500,15000,17500,20000),
	minnoise=0,
	maxnoise=0,
	initsij="small",
	EPS=1,
	EPSA=1,
	EPSB=1,
	chem="add"){




	 
	
	
	rt<-scan(rgc.file)

	SC<-scan(sc.file)
	
	
	
	#gradients into matrix form 
	rt<-matrix(rt,nrow=6,byrow=F)
	Nr=dim(rt)[2]
	
	SC<-matrix(SC,nrow=6,byrow=F)
	Nt=dim(SC)[2]

#############################

	#lookup table of pairwise distances 

	dist.retina<-rdist(t(rt[1:2,]),t(rt[1:2,]))
	
	q<-which(dist.retina<=0.0000000001,arr.ind=TRUE)

	dist.retina[q]=0

	dist.sc<-rdist(t(SC[1:2,]),t(SC[1:2,]))

	q<-which(dist.sc<=0.0000000001,arr.ind=TRUE)

	dist.sc[q]=0


#############################

	#get the connections data for Retina

	ret.sync.mat=mat.or.vec(Nr,Nr)
		
	for(e in 1:Nr){
		s<-which(dist.retina[e,]<=dist)
		ret.sync.mat[e,s]=1
	}


	
	max_neigh=max(rowSums(ret.sync.mat))
	
	inv_strength=rowSums(ret.sync.mat)

	connecs=mat.or.vec(max_neigh,Nr)

	for(xw in 1:Nr){
		sw=which(ret.sync.mat[xw,]==1)
		connecs[c(1:inv_strength[xw]),xw]=sw
	}

	#make it index from 0 (for C)

	connecs=connecs-1





	#same for SC

	bij.mat<-int.tec.ac(Nt,dist.sc,bij,dist.bij)

	max_n_bij=max(rowSums(bij.mat))

	inv_str_bij=rowSums(bij.mat)

	conn_bij=mat.or.vec(max_n_bij,Nt)

	for(xw in 1:Nt){
		sw=which(bij.mat[xw,]==1)
		conn_bij[c(1:inv_str_bij[xw]),xw]=sw
	}

	conn_bij=conn_bij-1

######################

	
	cij<-getcij(chem,rt[3,],SC[3,],rt[4,],SC[4,],Nr,Nt)
	

	sij<-getinitialsij(Nr,Nt,initsij,smallnumber=1,min=0.01,max=0.1)



###################
#Initiate simulation part


#running the first one outside of loop

	v=Nr*Nt

	

	if(length(stepbreaks)==1){
		if(stepbreaks<1){
	cat("0 interations specified-not valid",sep="\n")
		}
	a<-.C("run",as.integer(Nr),as.integer(connecs),as.integer(inv_strength),as.integer(max_neigh),
as.double(k),as.double(sij),as.integer(Nt),as.double(timestep),
as.double(cij),as.double(alpha),as.integer(v),as.double(maxnoise), as.double(minnoise),as.double(sijmin),as.integer(1), as.integer(stepbreaks),as.integer(conn_bij),as.integer(inv_str_bij),as.integer(max_n_bij))
	
	}else{
	
a<-.C("run",as.integer(Nr),as.integer(connecs),as.integer(inv_strength),as.integer(max_neigh),
as.double(k),as.double(sij),as.integer(Nt),as.double(timestep),
as.double(cij),as.double(alpha),as.integer(v),as.double(maxnoise), as.double(minnoise),as.double(sijmin),as.integer(1), as.integer(stepbreaks[1]),as.integer(conn_bij),as.integer(inv_str_bij),as.integer(max_n_bij))

	for(i in 1:(length(stepbreaks)-1)){


		sij<-matrix(a[[6]],nrow=Nr,byrow=F)

                cat(sprintf("%s About to run epoch %d to epoch %d (max %d)\n",
                            date(), stepbreaks[i]+1, stepbreaks[i+1],
                            stepbreaks[length(stepbreaks)]))
		a<-.C("run",as.integer(Nr),as.integer(connecs),as.integer(inv_strength),as.integer(max_neigh),as.double(k),as.double(sij),as.integer(Nt),as.double(timestep),
as.double(cij),as.double(alpha),as.integer(v),as.double(maxnoise), as.double(minnoise),as.double(sijmin),as.integer(stepbreaks[i]+1),as.integer(stepbreaks[i+1]),as.integer(conn_bij),as.integer(inv_str_bij),as.integer(max_n_bij))

		#write.table(matrix(a[[6]],nrow=Nr,byrow=F),file=paste("sijmatrix-",type,"most_recent.txt",sep=""),sep=" ")

	}

}
	sij=matrix(a[[6]],nrow=Nr,byrow=F)

	 write.table(signif(sij,4),file=op.file,sep=" ")






}

