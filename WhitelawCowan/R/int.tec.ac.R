int.tec.ac <-
function(Nt,dist.sc,bij,dist.bij){
	mat.bij<-mat.or.vec(Nt,Nt)
	q<-which(dist.sc<=dist.bij,arr.ind=TRUE)
	mat.bij[q]=1
	return(mat.bij)
}

