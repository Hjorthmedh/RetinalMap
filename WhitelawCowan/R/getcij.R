getcij <-
function(chem,Ai,aj,Bi,bj,Nr,Nt){
	EPSA=1
	EPSB=1
	EPS=1
	
	cij<-mat.or.vec(Nr,Nt)
	if(chem=="add"){
		for(i in 1:Nr){
			for(j in 1:Nt){
				cij[i,j]=EPSA*Ai[i]*(max(aj)-aj[j])+EPSB*Bi[i]*bj[j]
			}
		}
	}else if(chem=="mult"){
		for(i in 1:Nr){
			for(j in 1:Nt){
				cij[i,j]=EPS*Ai[i]*(max(aj)-aj[j])*Bi[i]*bj[j]
			}			
		}
	}
	cij=1+cij
	return(cij)
}

