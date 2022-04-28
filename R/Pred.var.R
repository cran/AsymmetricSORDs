
#' Function for generating the moment matrix and variance of the predicted response
#'
#' @param matrix Design matrix with the coefficients of the corresponding input factors
#'@description This function generates the moment matrix and variance of the predicted response for a given design based on a second-order model, for measuring the rotatability of the design. The input should be the specified form of a design matrix with the coefficients of the corresponding input factors. A minimum number of centre points is to be used to ensure the non-singularity of X`X.
#' @return The moment matrix and the prediction variance for a given design based on a second-order model
#' It gives unique prediction variance along with its frequencies.
#' @examples
#'  \dontrun{
#'library(AsymmetricSORDs)
#'Pred.var(matrix)
#'}
#'@references
#'1) G.E.P. Box and K.B. Wilson (1951).' On the experimental attainment of optimum conditions'.
#'

#'2) M. Hemavathi, Shashi Shekhar, Eldho Varghese, Seema Jaggi, Bikas Sinha & Nripes Kumar Mandal (2022)<DOI: 10.1080/03610926.2021.1944213>.' Theoretical developments in response surface designs: an informative review and further thoughts'.
#'@export
Pred.var<-function(matrix){
  mat<-matrix

  v<-ncol(mat)
  p=1
  while(p<=v){
    x1<-matrix(,nrow=nrow(mat),ncol=0)
    x1<-(mat[,p])^2
    mat<-cbind(mat,x1)
    p=p+1
  }


  b1=1
  b2=1
  while(b1<v){
    b2=b1+1
    while(b2<=v){
      mat2<-matrix(,nrow=0,ncol=1)
      mat2<-mat[,b1]*mat[,b2]
      mat<-cbind(mat,mat2)
      b2=b2+1
    }

    b1=b1+1

  }
  x_matrix<-mat

  x_matrix<-cbind(matrix(1,nrow=nrow(x_matrix),ncol=1),x_matrix)
  colnames(x_matrix)<-NULL
  x_prime_x<-t(x_matrix)%*%x_matrix
  moment_mat<-(1/nrow(x_matrix))*x_prime_x

   # print("Moment Matrix",quote=FALSE)
  message("Moment Matrix")
    print(moment_mat)



  k1=1
  var<-c()
  while(k1<=nrow(x_matrix))
  {
    V=t(x_matrix[k1,])
    b<-t(V)
    v_y_hat<-V %*%solve(x_prime_x) %*% b
    var<-c(var,v_y_hat)
    k1<-k1+1
  }
  variance_of_esitmated_response<-round(var,digits = 4 )
  print(table(variance_of_esitmated_response))
}
