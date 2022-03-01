Ibind_n_dim = function(p1, p2, ...,sigma =  sigma){
	C=t(chol(sigma))
	C_inv=solve(C)
	n_combine=dim(C)[1]

	der_erf=function(x) {2/(pi^0.5)*exp(-x^2)}
	der_norm=function(x) {1/(2*2^0.5)*der_erf(x/(2^0.5))}
	der_qnorm=function(x) {1/(der_norm(qnorm(x)))}
	input_x_text=paste("x",1:n_combine,collapse=",",sep="")
	for(i in 2:n_combine){
		eval(parse(text=(paste("der_var_",i,"=function(",paste("x",1:n_combine,collapse=",",sep=""),"){",
		"return(-der_norm(",paste(C_inv[i,],paste("qnorm(1-exp(",paste("x",1:n_combine,sep=""),"))",sep=""),sep="*",collapse="+"),"))}",sep=""))))
		eval(parse(text=(paste("der_var_",i,"=Vectorize(der_var_",i,")",sep=""))))
	}

	eval(parse(text=(paste("test_f_nd=function(",input_x_text,"){",
		"if(",paste("is.na(qnorm(1-exp(",paste("x",1:n_combine,sep=""),")))",sep="",collapse=" | "),"){return(0)};",
		"return(abs(",paste("(1/",diag(C),")",collapse="*",sep=""),"*",paste("der_var_",2:n_combine,"(",input_x_text,")",sep="",collapse="*"),"*",
			paste("der_qnorm(1-exp(",paste("x",2:n_combine,sep=""),"))",sep="",collapse="*"),"*",
			paste("exp(",paste("x",1:n_combine,sep=""),")",sep="",collapse="*"),"))}",sep=""))))
	test_f_nd=Vectorize(test_f_nd)

	foreheadint_text=NULL
	for(i in 1:n_combine){
		if(!i==n_combine){
			foreheadint_text=paste(foreheadint_text,"integrate(function(x",i,") {",
			"sapply(x",i,", function(x",i,") {",sep="")
		}else{
			foreheadint_text=paste(foreheadint_text,"integrate(function(x",i,") {",
			"sapply(x",i,", function(x",i,") test_f_nd(",input_x_text,"))",sep="")
		}
	}
	tailint_text=NULL
	for(i in 1:n_combine){
		if(!i==n_combine){
			tailint_text=paste(tailint_text,"}, ",
			paste("log(p",1:n_combine,"[xx])",sep="",collapse="+"),"-",paste("x",1:(n_combine-i),sep="",collapse="-"),",0)$value})",sep="")
		}else{
			tailint_text=paste(tailint_text,"}, ",
			paste("log(p",1:n_combine,"[xx])",sep="",collapse="+"),",0)$value)",sep="")
		}
	}

	eval(parse(text=(paste("int_f_nd=function(xx){",
		"if(",paste("is.na(",paste("p",1:n_combine,sep=""),"[xx])",sep="",collapse=" | "),"){return(NA)};",
		"return(1-",foreheadint_text,tailint_text,"}",sep=""))))
		
	temp=unlist(lapply(1:length(p1), function(xx) int_f_nd(xx)))
	temp=ifelse(temp<0,0,temp)
	return(temp)
}


#example
r12 = 0.2
r13 = 0.3
r23 = 0.6
COV = matrix(c(1, r12, r13, r12, 1, r23, r13, r23, 1), 3, 3)
p1 = p2 = p3 = c(0.04, 0.05)
Ibind_n_dim(p1, p2, p3, sigma = COV)
#p_I = 0.01518269 0.02151095

r12 = 0.2
r13 = 0.3
r14 = 0.4
r23 = 0.6
r24 = 0.2
r34 = 0.1
COV = matrix(c(1, r12, r13, r14, r12, 1, r23, r24, r13, r23, 1, r34, r14, r24, r34, 1), 4, 4)
p1 = p2 = p3 = p4 = 0.05
#Ibind_n_dim(p1, p2, p3, p4, sigma=COV)
#p_I = 0.01404982
#0.01404982 < 0.02151095 (<0.05) means the larger synergy happened as # of p-value to be combined.
