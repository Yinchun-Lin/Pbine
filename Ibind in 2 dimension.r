der_erf=function(x) {2/(pi^0.5)*exp(-x^2)}
der_norm=function(x) {1/(2*2^0.5)*der_erf(x/(2^0.5))}
#Imprtant!!! Since: der_qnorm(1-exp(y)) was exchanged to der_qnorm(exp(y)) in test_f, qnorm(x) -> qnorm(x,lower.tail=FALSE)
der_qnorm=function(x) {1/(der_norm(qnorm(x,lower.tail=FALSE)))}

test_f=function(x,y,r,w=NULL){
  rr=(1-r^2)^0.5
  if(is.null(w)){
	temp_qnorm_x=qnorm(exp(x),lower.tail=FALSE)
	temp_qnorm_y=qnorm(exp(y),lower.tail=FALSE)
    if(is.na(temp_qnorm_y) | is.na(temp_qnorm_x) | is.infinite(temp_qnorm_y) | is.infinite(temp_qnorm_x)){
	  print("KK")
      return(0)
    }
	#print("YY")
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

integral_method=function(p1, p2, r=NULL, method="int", w=NULL, max_exp=17, rel_error=1e-5){
  if(is.null(r)){
	sell=(!is.na(p1))&(!is.na(p2))
    r=cor(p1[sell],p2[sell])
    if(is.na(r)){stop("please enter correlation r")}
  }
  test_f=Vectorize(test_f)
  if(r==1) stop("correlation r must be less than 1")
  if(method=="int"){
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
		return(int_f(xx))}))
    temp=ifelse(temp<0,0,temp)
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

#example
r=0.5
rr=(1-r^2)^0.5
p1=runif(10)
p2=runif(10)
p2=1-pnorm(rr*qnorm(1-p2)+r*qnorm(1-p1))
integral_method(p1,p2,method="int",r=r)


