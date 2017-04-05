ctmcdlogLik <-
function(gm,tmabs,te){
  P=expm(gm*te)
  ll=0
  for(i in 1:nrow(P)){
    for(j in 1:ncol(P)){
      if(P[i,j]>0){
        ll=ll+tmabs[i,j]*log(P[i,j])
      }
    }
  }
  return(ll)
}
