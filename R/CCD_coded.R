#' Central Composite Designs (CCD) with coded levels
#' @param v Number of input factors, v(>2)
#' @param type Type of central composite design i.e.  ccc or cci or ccf. "ccc" is for Central Composite Circumscribed designs, "cci" is for Central Composite Inscribed designs and "ccf" is for Central Composite Face Centered designs
#' @param randomization It is for generating the randomized layout of the design. It takes either TRUE or FALSE  and by default, it is set to FALSE
#' @param variance This is for generating the moment matrix and prediction variance of the design based on a second order model. It gives unique prediction variance along with its frequencies. It takes either TRUE or FALSE  and by default, it is set to FALSE
#'@description This function generates Central Composite Designs (CCD) with coded levels for a given number of input factors (v). The CCD constitute  combinations of factorial points, axial points and center points. Three types of CCD can be generated using this function  i.e.  ccc or cci or ccf. "ccc" is for Central Composite Circumscribed designs, "cci" is for Central Composite Inscribed designs and "ccf" is for Central Composite Face Centered designs. It gives the randomized layout of the design along with the moment  matrix and prediction variance.
#' @return Central Composite Designs (CCD) for a given number of input factors (v) with coded levels
#' @note Here, the factorial portion consists of  2^v (full factorial) combinations and there is no upper limit for the number of input factors, v (>2). To get a CCD with smaller runs,  one may use fractional factorial (of resolution V) in place of full factorial.
#' @examples
#'
#'library(AsymmetricSORDs)
#'CCD_coded(5,'ccc',FALSE,FALSE)
#'CCD_coded(6,"cci",FALSE,FALSE)
#'
#'@references
#'1) G.E.P. Box and K.B. Wilson (1951)." On the experimental attainment of optimum conditions".
#'

#'2) M. Hemavathi, Shashi Shekhar, Eldho Varghese, Seema Jaggi, Bikas Sinha & Nripes Kumar Mandal (2022)<DOI: 10.1080/03610926.2021.1944213>. "Theoretical developments in response surface designs: an informative review and further thoughts".
#'@export
CCD_coded<-function(v,type,randomization=FALSE,variance=FALSE){
  v
  type
  randomization
  variance
  k<-2^v
  a<-(k)^(1/4)
  cbn<-(factorial(v))/(factorial(2)*factorial((v-2)))
  nc<-4*(2^(v/2))+4-2*v
  nc<-round(nc)

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
  p=1
  while(p<=v){
    x1<-matrix(,nrow=nrow(matfAc),ncol=0)
    x1<-(matfAc[,p])^2
    matfAc<-cbind(matfAc,x1)
    p=p+1
  }


  b1=1
  b2=1
  while(b1<v){
    b2=b1+1
    while(b2<=v){
      mat2<-matrix(,nrow=0,ncol=1)
      mat2<-matfAc[,b1]*matfAc[,b2]
      matfAc<-cbind(matfAc,mat2)
      b2=b2+1
    }

    b1=b1+1

  }
  ###########################

  matfAc<-cbind(matrix(1,nrow=nrow(matfAc),ncol=1),matfAc)
  colnames(matfAc)<-NULL

  #
  cccsd<-matfAc

  if(randomization==T){
    mat3<-matrix(,nrow=0,ncol=1+v+v+cbn)
    rand<-sample(1:(k+(2*v)+nc),(k+(2*v)+nc),replace = FALSE)

    for(m in rand){
      x11<-matrix(,nrow=1,ncol=1+v+cbn+v)
      x11<-matfAc[m,]
      mat3<-rbind(mat3,x11)
    }
    ##########################
    matfAc<-mat3
    rownames(matfAc)<-NULL

    cccrand<-matfAc

  }
  ######################################x prime x
  if(type=="ccc"){
    if(randomization==FALSE){
     # print("CCC (Standard)",quote=F)
      message("CCC (Standard)")
      print(cccsd[,2:(v+1)])
      matfAc<-cccsd
    }else{
     # print("CCC (Randomized)",quote=F)
      message("CCC (Randomized)")
      print(cccrand[,2:(v+1)])
      matfAc<-cccrand
    }

    x_prime_x<-t(matfAc)%*% matfAc

    #################################variance

    k1=1
    var<-c()
    while(k1<=(k+2*v+nc))
    {
      V=t(matfAc[k1,])
      b<-t(V)
      v_y_hat<-V %*%solve(x_prime_x) %*% b
      var<-c(var,v_y_hat)
      k1<-k1+1
    }
    if(variance==TRUE){

      variance_of_esitmated_response<-round(var,digits = 3 )
      v1<-table(variance_of_esitmated_response)
    }
  }
  ########################
  #######CCI
  if(type=="cci"){
    cci<-matfAc[,2:(v+1)]/a
    #############adding columns(sq terms + interactions)
    p=1
    while(p<=v){
      x1<-matrix(,nrow=nrow(cci),ncol=0)
      x1<-(cci[,p])^2
      cci<-cbind(cci,x1)
      p=p+1
    }
    ##############################
    b1=1
    b2=2
    while(b1<b2){

      while(b2<=v){
        mat2<-matrix(,nrow=k,ncol=0)
        mat2<-cci[,b1]*cci[,b2]
        cci<-cbind(cci,mat2)
        b2=b2+1
      }
      b2=b2-1
      b1=b1+1
    }
    ###########################

    cci<-cbind(matrix(1,nrow=nrow(cci),ncol=1),cci)
    colnames(cci)<-NULL





    ##########
    x_prime_x<-t(cci)%*% cci
    k1=1
    var<-c()
    while(k1<=k+2*v+nc)
    {
      V=t(cci[k1,])
      b<-t(V)
      v_y_hat<-V %*%solve(x_prime_x) %*% b
      var<-c(var,v_y_hat)
      k1<-k1+1
    }
    if(randomization==T){
      #print("CCI (Randomized)",quote=F)
      message("CCI (Randomized)")
    }else{
     #print("CCI (Standard)",quote=F)
      message("CCI (Standard)")
    }
    print(cci[,2:(v+1)])
    if(variance==TRUE){
      #if(unique==FALSE){
      variance_of_esitmated_response<-round(var,digits = 3 )
      v2<-table(variance_of_esitmated_response)
    }

  }
  ############CCF
  if(type=="ccf"){
    ccf<-matfAc[,2:(v+1)]
    ccf[ccf==a]<-1
    ccf[ccf==-a]<- -1

    #############adding columns(sq terms + interactions)
    p=1
    while(p<=v){
      x1<-matrix(,nrow=nrow(ccf),ncol=0)
      x1<-(ccf[,p])^2
      ccf<-cbind(ccf,x1)
      p=p+1
    }
    ##############################
    b1=1
    b2=2
    while(b1<b2){

      while(b2<=v){
        mat2<-matrix(,nrow=k,ncol=0)
        mat2<-ccf[,b1]*ccf[,b2]
        ccf<-cbind(ccf,mat2)
        b2=b2+1
      }
      b2=b2-1
      b1=b1+1
    }
    ###########################

    ccf<-cbind(matrix(1,nrow=nrow(ccf),ncol=1),ccf)
    colnames(ccf)<-NULL

    ##########
    x_prime_x<-t(ccf)%*% ccf
    k1=1
    var<-c()
    while(k1<=k+2*v+nc)
    {
      V=t(ccf[k1,])
      b<-t(V)
      v_y_hat<-V %*%solve(x_prime_x) %*% b
      var<-c(var,v_y_hat)
      k1<-k1+1
    }
    variance_of_esitmated_response<-round(var,digits = 3 )
    v3<-table(variance_of_esitmated_response)
    if(randomization==T){
     # print("CCF (Randomized)",quote=F)
      message("CCF (Randomized)")
    }else{
      #print("CCF (Standard)",quote=F)
      message("CCF (Standard)")
    }

    print(ccf[,2:(v+1)])
  }
    if(type=="ccf"){
    final<-ccf[,2:(v+1)]
    }
    if(type=="cci"){
    final<-cci[,2:(v+1)]
    }
    if(type=="ccc"){
    final<-matfAc[,2:(v+1)]
    }
    if(variance==TRUE){
      x_matrix<-final

      x_matrix<-cbind(matrix(1,nrow=nrow(x_matrix),ncol=1),x_matrix)
      colnames(x_matrix)<-NULL
      x_prime_x<-t(x_matrix)%*%x_matrix
      moment_mat<-(1/nrow(x_matrix))*x_prime_x

      #print("Moment Matrix",quote=FALSE)
      message("Moment Matrix")
      print(moment_mat)

      if(type=="ccf"){
        print(v3)
      }
      if(type=="cci"){
        print(v2)
      }
      if(type=="ccc"){
        print(v1)
      }


      }

}

