Vk <- function(l,i){ # l: state, i: position
  sequence <- measures.sequence
  if (i==0)
    return(1)
  else{
    e.li = e[l, sequence[i]]
    
    vk.max = 0
    for (k in states){
      vk.temp <- Vk (k, i-1) * a[k,l]
      if (vk.temp > vk.max){
        vk.max <- vk.temp
        path[i]<<- k
      }
    }
    return(e[l, sequence[i]]*vk.max)
  }
}