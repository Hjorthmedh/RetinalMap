getinitialsij <-
function(Nr,Nt,type,smallnumber=1,min=0.00005,max=0.00011){
	switch(type,
		zero=mat.or.vec(Nr,Nt),
		small=smallnumber*(mat.or.vec(Nr,Nt)+1),
		rand=matrix(runif(Nr*Nt,min,max),nrow=Nr)
	)
}

