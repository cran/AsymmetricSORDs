#' Asymmetric Second Order Rotatable Designs
#' @param v Number of input factors, v(>2)
#' @param number_of_pairs Number of pairs of input factors for which asymmetry is required
#' @param z A vector of real number and its length equals to number_of_pairs
#' @param type Type of central composite design i.e.  ccc or cci. "ccc" is for Central Composite Circumscribed designs and "cci" is for Central Composite Inscribed designs
#' @param randomization It is for generating the randomized layout of the design. It takes either TRUE or FALSE  and by default, it is set to FALSE
#' @param variance This is for generating the moment matrix and prediction variance of the design based on a second order model. It gives unique prediction variance along with its frequencies. It takes either TRUE or FALSE  and by default, it is set to FALSE
#'@description This function generates  ASORDs through the orthogonal transformation of central composite designs as per the procedure given by J.S. Mehta and M.N. Das (1968). It would be providing two types of asymmetric designs for a given number of treatments (v). It requires four input parameters viz., v(>2); number_of_pairs(>0); z= vector of real number of length equals to number_of_pairs; type="ccc" or "cci" and randomization=TRUE or FALSE.
#' @return
#' Asymmetric Second Order Rotatable Designs (ASORDs) for a given v.
#'
#' @examples
#'
#'library(AsymmetricSORDs)
#'Asords(5,2,c(2,3),"ccc",TRUE)
#'
#'@references
#'1) J.S. Mehta and M.N. Das (1968)." Asymmetric rotatable designs and orthogonal transformations".
#'

#'2)M. Hemavathi, Eldho Varghese, Shashi Shekhar & Seema Jaggi (2020)<DOI: 10.1080/02664763.2020.1864817>." Sequential asymmetric third order rotatable designs (SATORDs)".
#'

#'3) M. Hemavathi, Shashi Shekhar, Eldho Varghese, Seema Jaggi, Bikas Sinha & Nripes Kumar Mandal (2022)<DOI: 10.1080/03610926.2021.1944213>." Theoretical developments in response surface designs: an informative review and further thoughts".

#'@export
Asords<-function(v,number_of_pairs,z,type,randomization=FALSE,variance=FALSE){

  np<-number_of_pairs

  B<-matrix(0,nrow=v,ncol=v)
  s=1
  ss=1
  while(s<=np*2){
    d<-z[ss]
    D<-matrix(c(1/(sqrt(1+d^2)),d/(sqrt(1+d^2)),-d/(sqrt(1+d^2)),1/(sqrt(1+d^2))),ncol=2,byrow=T)
    B[s:(s+1),s:(s+1)]<-t(D)
    s=s+2
    ss=ss+1
  }
  nxt<-s
  while(nxt<=v){
    B[nxt,nxt]<-1
    nxt=nxt+1
  }


  #variance=F
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

  ccc<-matfAc
  cci<-matfAc
  ccf<-matfAc

  ##########################

  #######################################randomization
  if(randomization==T){
    mat3<-matrix(,nrow=(k+(2*v)+nc),ncol=v)
    rand<-sample(1:(k+(2*v)+nc),(k+(2*v)+nc),replace = FALSE)

    for(m in rand){
      x11<-matrix(,nrow=1,ncol=v)
      x11<-matfAc[m,]
      mat3<-rbind(mat3,(x11))
    }
    ##########################
    matfAc<-mat3
    rownames(matfAc)<-NULL

    #if(randomization==TRUE){
    ccc<-matfAc[((k+(2*v)+nc+1)):(2*(k+(2*v)+nc)),]
    cci<-matfAc[((k+(2*v)+nc+1)):(2*(k+(2*v)+nc)),]
    ccf<-matfAc[((k+(2*v)+nc+1)):(2*(k+(2*v)+nc)),]
    #final<-matfAc
  }

  # }
  ####################CCC

  ######################
  ######################################x prime x
  if(type=="ccc"){
    final<-ccc
  }


  ########################
  #######CCI
  if(type=="cci"){

    cci<-cci/a
    final<-cci

  }
  ############CCF
  if(type=="ccf"){

    #ccf<-matfAc
    ccf[ccf==a]<-1
    ccf[ccf==-a]<- -1
    final<-ccf
  }
  #############adding columns(sq terms + interactions)

  ###########################


  #print("Asymmetric Rotatable Design",quote=F)
  message("Asymmetric Rotatable Design")
  final<-final%*%B
  print(final)

  ###############
  p=1
  while(p<=v){
    x1<-matrix(,nrow=nrow(final),ncol=0)
    x1<-(final[,p])^2
    final<-cbind(final,x1)
    p=p+1
  }


  b1=1
  b2=1
  while(b1<v){
    b2=b1+1
    while(b2<=v){
      mat2<-matrix(,nrow=0,ncol=1)
      mat2<-final[,b1]*final[,b2]
      final<-cbind(final,mat2)
      b2=b2+1
    }

    b1=b1+1

  }

  ##########
  x_prime_x<-t(final)%*% final
  moment_mat<-(1/nrow(final))*x_prime_x
  rownames(moment_mat)<-NULL
  colnames(moment_mat)<-NULL
  k1=1
  var<-c()
  while(k1<=k+2*v+nc){
    V=t(final[k1,])
    b<-t(V)
    v_y_hat<-V %*%solve(x_prime_x) %*% b
    var<-c(var,v_y_hat)
    k1<-k1+1
  }

  if(variance==TRUE){
    #print("Moment Matrix",quote=FALSE)
    message("Moment Matrix")
    print(moment_mat)
    variance_of_esitmated_response<-round(var,digits = 3 )
    print(table(variance_of_esitmated_response))
  }
}
