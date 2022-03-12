der_erf = function(x) {2/(pi^0.5)*exp(-x^2)}
der_norm = function(x) {1/(2*2^0.5)*der_erf(x/(2^0.5))}
#Important!!! Since: der_qnorm(1-exp(y)) was exchanged to der_qnorm(exp(y)) in test_f, qnorm(x) -> qnorm(x,lower.tail=FALSE)
der_qnorm = function(x) {1/(der_norm(qnorm(x,lower.tail=FALSE)))}

test_f=function(x,y,r,w=NULL){
	rr=(1-r^2)^0.5
	if(is.null(w)){
		temp_qnorm_x=qnorm(exp(x),lower.tail=FALSE)
		temp_qnorm_y=qnorm(exp(y),lower.tail=FALSE)
	if(is.na(temp_qnorm_y) | is.na(temp_qnorm_x) | is.infinite(temp_qnorm_y) | is.infinite(temp_qnorm_x)){
		return(0)
	}
		return(abs(-der_norm((temp_qnorm_y-r*temp_qnorm_x)/rr)*(-der_qnorm(exp(y))/rr))*exp(x)*exp(y))
	}else{
		temp_qnorm_x=qnorm(exp(x/w),lower.tail=FALSE)
		temp_qnorm_y=qnorm(exp(y/(1-w)),lower.tail=FALSE)
	if(is.na(temp_qnorm_y) | is.na(temp_qnorm_x) | is.infinite(temp_qnorm_y) | is.infinite(temp_qnorm_x)){
		return(0)
	}
		return(abs(-der_norm((temp_qnorm_y-r*temp_qnorm_x)/rr)*(-der_qnorm(exp(y/(1-w)))/rr))*exp(x/w)*exp(y/(1-w))*(1/w)*(1/(1-w)))
	}
}
test_f=Vectorize(test_f)

fisher=function(x,y){
	return(x*y-x*y*log(x*y))
}

integral_method=function(pm, r=NULL, method="Int", w=NULL, max_exp=17, rel_error=1e-5){
	if(is.null(r)){
		sell=(!is.na(p1))&(!is.na(p2))
		r=cor(p1[sell],p2[sell])
		if(is.na(r)){stop("please enter correlation r")}
	}
	test_f=Vectorize(test_f)
	
	if(r==1) stop("correlation r must be less than 1")
	
	if(method=="Int"){
		for(i in 1:2){
			eval(parse(text = (paste("p", i, "=pm[ ,", i, "]", sep = ""))))
		}
		
		int_f=function(xx){
			if(is.na(p1[xx])|is.na(p2[xx])){
				return(NA)
			}else if((p1[xx]*p2[xx])<=10^(-max_exp)){
				return(0)
			}
			else{
				return(1-integrate(function(y) {
					sapply(y, function(y) {
						integrate(function(x) {
							sapply(x, function(x) test_f(x,y,r))
						}, log(p1[xx])+log(p2[xx])-y, 0,rel.tol=rel_error)$value
					})
				}, log(p1[xx])+log(p2[xx]), 0,rel.tol=rel_error)$value)
			}
		}
		
		temp=unlist(lapply(1:length(p1), function(xx) {
			if(xx%in%round(quantile(1:length(p1),(1:10)/10))){print(paste("Progress",round(xx/length(p1)*100,2)))}
			return(int_f(xx))
			}))
		temp=ifelse(temp<0, 0, temp)
	}else if(method=="fisher"){
		temp=fisher(p1,p2)
	}else if(method=="decor"){
		rr=(1-r^2)^0.5
		temp=fisher(p1,1-pnorm((qnorm(1-p2)-r*qnorm(1-p1))/rr))
	}else if(method=="int_weighted"){
		temp=unlist(lapply(1:length(p1), function(xx)
			1-integrate(function(y) {
				sapply(y, function(y) {
					integrate(function(x) {
						sapply(x, function(x) test_f(x,y,r,w=w))
					}, w*log(p1[xx])+(1-w)*log(p2[xx])-y, 0)$value
				})
			}, w*log(p1[xx])+(1-w)*log(p2[xx]), 0)$value))
		temp=ifelse(temp<0,0,temp)
	}
	return(temp)
}

Ibind_n_dim = function(pm, sigma){
	C=t(chol(sigma))
	C_inv=solve(C)
	n_combine=dim(C)[1]
	for(i in 1:n_combine){
		eval(parse(text = (paste("p", i, "=pm[ ,", i, "]", sep = ""))))
	}
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

Pbine=function(pm, sigma = NULL, method = "Int", w = NULL){
	if(is.null(dim(pm))) {stop("pm must be a matrix.")}
	if(dim(pm)[2] == 1) {stop("pm must be a n by m matrix with m > 1.")}
	if(dim(pm)[2] == 2) {return(integral_method(pm, r = sigma, method = method, w = w))}
	if(dim(pm)[2] > 2) {return(Ibind_n_dim(pm=pm, sigma = sigma))}
}

#example
#1. The case of n > 1.
set.seed(123)
p1 = runif(20)
p2 = (p1 + runif(20))/2 # p2 is correlated with p1.
pm = cbind(p1, p2)
pn = Pbine(pm) # sigma is estimated by cor(p1, p2).
print(cbind(p1, p2, pn))

#2. The case of n = 1. sigma should be provided.
r12 = 0.2
r13 = 0.3
r23 = 0.6
COV = matrix(c(1, r12, r13, r12, 1, r23, r13, r23, 1), 3, 3)
p1 = p2 = p3 = 0.05
pm = t(c(p1, p2, p3))
pn = Pbine(pm, sigma = COV, method = "Int") # If n = 1, sigma should be provided.
print(c(p1, p2, pn))

#3. The case of m > 2.
r12 = 0.2
r13 = 0.3
r14 = 0.4
r23 = 0.6
r24 = 0.2
r34 = 0.1
COV = matrix(c(1, r12, r13, r14, r12, 1, r23, r24, r13, r23, 1, r34, r14, r24, r34, 1), 4, 4)
p1 = p2 = p3 = p4 = 0.05
pm = t(c(p1, p2, p3, p4))
Pbine(pm, sigma = COV)
# 0.01404982
Pbine(pm, method = "Fisher")
# 0.0023
# NOTE 1: In this example, the Fisher's method overestimates the significance because p-values are correlated.
# NOTE 2: Please wait a while because it needs more time when m > 3.











