#' Central Composite Designs (CCD) with original levels
#'
#' @param v Number of input factors, v(>2)
#' @param type Type of central composite design i.e.  ccc or cci or ccf. "ccc" is for Central Composite Circumscribed designs, "cci" is for Central Composite Inscribed designs and "ccf" is for Central Composite Face Centered designs
#' @param randomization It is for generating the randomized layout of the design. It takes either TRUE or FALSE  and by default, it is set to FALSE
#' @param variance This is for generating the moment matrix and prediction variance of the design based on a second order model. It gives unique prediction variance along with its frequencies. It takes either TRUE or FALSE  and by default, it is set to FALSE
#'@param min_L A vector of minimum levels of the factors
#'@param max_L A vector of maximum levels of the factors
#'@description This function generates Central Composite Designs (CCD) with original levels along with coded levels for a given number of input factors (v). The CCD constitute  combinations of factorial points, axial points and center points. Three types of CCDs can be generated using this function  i.e.  ccc or cci or ccf. "ccc" is for Central Composite Circumscribed designs, "cci" is for Central Composite Inscribed designs and "ccf" is for Central Composite Face Centered designs. It gives the randomized layout of the design along with the moment  matrix and prediction variance.
#' @return Central Composite Designs (CCD) for a given number of input factors (v) with original levels
#' @note Here, the factorial portion consists of 2^v (full factorial) combinations and there is no upper limit for the number of input factors,v (>2). To get a CCD with smaller runs,  one may use fractional factorial (of resolution V) in place of full factorial.
#' @examples
#'
#'library(AsymmetricSORDs)
#'CCD_original(5,'ccc',c(10,15,20,25,30),c(15,20,25,30,35),FALSE,FALSE)
#'
#'@references
#'1) G.E.P. Box and K.B. Wilson (1951)." On the experimental attainment of optimum conditions".
#'

#'2) M. Hemavathi, Shashi Shekhar, Eldho Varghese, Seema Jaggi, Bikas Sinha & Nripes Kumar Mandal (2022)<DOI: 10.1080/03610926.2021.1944213>. "Theoretical developments in response surface designs: an informative review and further thoughts".

#'@export
CCD_original<-function(v,type,min_L,max_L,randomization=FALSE,variance=FALSE){

  k<-2^v
  a<-(k)^(1/4)
  cbn<-(factorial(v))/(factorial(2)*factorial((v-2)))
  nc<-4*(2^(v/2))+4-2*v
  nc<-round(nc)
  aa<-c(-((2^v)^(1/4)), -1, 0, 1,((2^v)^(1/4)))
  a_CCC<-aa  ##########Circumscribed (CCC)#############
  a_CCI<-aa/((2^v)^(1/4))##########Inscribed (CCI)#############
  a_CCF<-c(-1, -1, 0, 1,1)##########Face Centered (CCF)#############
  ####################################CCC
  if(type=="ccc"){
    x_ccc<-matrix(,length(min_L),ncol=5)

    for (j in 1:length(min_L))
    {
      for (i in 1:5) {
        x_ccc[j,i]<-(((max_L[j]-min_L[j])/2)*a_CCC[i])+((min_L[j]+max_L[j])/2)
      }
    }
  }
  #########################CCI
  if(type=="cci"){
    x_cci<-matrix(,length(min_L),ncol=5)

    for (j in 1:length(min_L))
    {
      for (i in 1:5) {
        x_cci[j,i]<-(((max_L[j]-min_L[j])/2)*a_CCI[i])+((min_L[j]+max_L[j])/2)
      }
    }
  }
  #############################CCF
  if(type=="ccf"){
    x_ccf<-matrix(,length(min_L),ncol=5)

    for (j in 1:length(min_L))
    {
      for (i in 1:5) {
        x_ccf[j,i]<-(((max_L[j]-min_L[j])/2)*a_CCF[i])+((min_L[j]+max_L[j])/2)
      }
    }

  }
  ###########################factorial
  i=1
  matf<-matrix(,nrow=k,ncol=0)
  while(i<=v){
    x<-c((rep(-1,(2^(v-i)))),(rep(1,(2^(v-i)))))
    x1<-t(rep(x,k/length(x)))
    matf<-cbind(matf,t(x1))
    i=i+1
  }
  #######################################axial

  matA<-matrix(,nrow=0,ncol=v)
  d1<-diag(a,nrow=v,ncol=v)
  d2<-diag(-a,nrow=v,ncol=v)
  d=1
  while(d<=v){
    matA<-rbind(matA,d1[d,],d2[d,])
    d=d+1
  }
  matA<-matA*(-1)
  #######################################
  matfA<-rbind(matf,matA)
  ############################central

  matfAc<-rbind(matfA,matrix(0,nrow=nc,ncol=v))
  #################################square terms

  ###########################
  cccsd<-matfAc

  ##########################

  #######################################randomization

  mat3<-matrix(,nrow=(k+(2*v)+nc),ncol=v)
  rand<-sample(1:(k+(2*v)+nc),(k+(2*v)+nc),replace = FALSE)

  for(m in rand){
    x1<-matrix(,nrow=1,ncol=v)
    x1<-matfAc[m,]
    mat3<-rbind(mat3,t(x1))
  }
  ##########################
  matfAc<-mat3[((k+(2*v)+nc+1):(2*(k+(2*v)+nc))),]
  rownames(matfAc)<-NULL
  cccrand<-matfAc

  ####################CCC

  ######################
  ######################################x prime x
  if(type=="ccc"){
    if(randomization==FALSE){
      #print("CCC coded (Standard)")
      message("CCC coded (Standard)")
      print(cccsd)
      matfAc<-cccsd}else{
        #print("CCC coded (Randomized)")
        message("CCC coded (Randomized)")
        print(cccrand)
        matfAc<-cccrand
      }
  }

  ########################
  #######CCI
  if(type=="cci"){
    cci<-matfAc/a
    #############adding columns(sq terms + interactions)
    if(randomization==T){
      #print("CCI coded (Randomized)")
      message("CCI coded (Randomized)")
    }else{
      #print("CCI coded (Standard)")
      message("CCI coded (Standard)")
    }
    print(cci)

  }

  ############CCF
  if(type=="ccf"){
    ccf<-matfAc
    ccf[ccf==a]<-1
    ccf[ccf==-a]<- -1
    if(randomization==T){
      #print("CCF coded (Randomized)")
      message("CCF coded (Randomized)")
    }else{
      #print("CCF coded (Standard)")
      message("CCF coded (Standard)")
    }
    print(ccf)
  }

  if(type=="ccc"){
    raw<-matfAc
  }
  if(type=="cci"){
    raw<-cci
  }
  if(type=="ccf"){
    raw<-ccf
  }
  s=1

  while(s<=v){
    if(type=="ccc"){

      o<-c(x_ccc[s,])

      raw[which(raw[,s]==0),s]<-o[3]
      raw[which(raw[,s]==-1),s]<-o[2]
      raw[which(raw[,s]==-a),s]<-o[1]
      raw[which(raw[,s]==1),s]<-o[4]
      raw[which(raw[,s]==a),s]<-o[5]
      s=s+1
    }
    if(type=="cci"){

      o<-c(x_cci[s,])

      raw[which(raw[,s]==0),s]<-o[3]
      raw[which(raw[,s]==-1),s]<-o[1]
      raw[which(raw[,s]==-1/a),s]<-o[2]
      raw[which(raw[,s]==1),s]<-o[5]
      raw[which(raw[,s]==1/a),s]<-o[4]
      s=s+1
    }
    if(type=="ccf"){

      o<-c(x_ccf[s,])

      raw[which(raw[,s]==0),s]<-o[3]
      raw[which(raw[,s]==-1),s]<-o[1]
      raw[which(raw[,s]==-1),s]<-o[2]
      raw[which(raw[,s]==1),s]<-o[5]
      raw[which(raw[,s]==1),s]<-o[4]
      s=s+1
    }
  }
  if(type=="ccc"){
    if(randomization==T){
      #print("CCC original (Randomized)")
      message("CCC original (Randomized)")
    }else{
      #print("CCC original (Standard)")
      message("CCC original (Standard)")
    }
  }
  if(type=="cci"){
    if(randomization==T){
      #print("CCI original (Randomized)")
      message("CCI original (Randomized)")
    }else{
      #print("CCI original (Standard)")
      message("CCI original (Standard)")
    }
  }
  if(type=="ccf"){
    if(randomization==T){
      #print("CCF original (Randomized)")
      message("CCF original (Randomized)")
    }else{
     # print("CCF original (Standard)")
      message("CCF original (Standard)")
    }
  }
  print(raw)
  #############################

  p=1
  while(p<=v){
    x1<-matrix(,nrow=nrow(raw),ncol=0)
    x1<-(raw[,p])^2
    raw<-cbind(raw,x1)
    p=p+1
  }


  b1=1
  b2=1
  while(b1<v){
    b2=b1+1
    while(b2<=v){
      mat2<-matrix(,nrow=0,ncol=1)
      mat2<-raw[,b1]*raw[,b2]
      raw<-cbind(raw,mat2)
      b2=b2+1
    }

    b1=b1+1

  }
  x_matrix<-raw

  x_matrix<-cbind(matrix(1,nrow=nrow(x_matrix),ncol=1),x_matrix)
  colnames(x_matrix)<-NULL
  x_prime_x<-t(x_matrix)%*%x_matrix
  moment_mat<-(1/nrow(x_matrix))*x_prime_x
  if(variance==TRUE){
    #print("Moment Matrix",quote=FALSE)
    message("Moment Matrix")
    print(moment_mat)
    k1=1
    var<-c()
    while(k1<=k+2*v+nc){
      V=t(x_matrix[k1,])
      b<-t(V)
      v_y_hat<-V%*%solve(x_prime_x) %*% b
      var<-c(var,v_y_hat)
      k1<-k1+1
    }
    variance_of_esitmated_response<-round(var,digits = 3 )
    print(table(variance_of_esitmated_response))
  }

}
