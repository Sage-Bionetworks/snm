downweight.outliers <-
function(mat){ 
  mns <- rowMeans(mat)
  mat.c <- mat - mns
  mat.c2 <- 0 * mat.c
  mat.c2 <- apply(mat.c,2,function(x){
    sd <- mad(x)
    md <- median(x)
    x[which(abs(x - md) > (5*sd))] <- md;
    x
  })
  mns + mat.c2		
}

