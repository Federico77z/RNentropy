RN_pmi <- function(Results)
{
  if(is.null(Results$selected))
    Results <- RN_select(Results)
  
  Results$pmi <- matrix(nrow = ncol(Results$design), 
                        ncol = ncol(Results$design))
  
  colnames(Results$pmi) <- colnames(Results$design)
  rownames(Results$pmi) <- colnames(Results$design)
  
  Results$npmi <- Results$pmi
  
  colshift <- ncol(Results$selected) - ncol(Results$design) 
  
  for(x in 1:nrow(Results$pmi))
  {
    for(y in 1:nrow(Results$pmi))
    {
      if(x > y)
      {
        Results$pmi[x,y] <- Results$pmi[y,x]
        Results$npmi[x,y] <- Results$npmi[y,x]
        next
      } else
      {
        sum_x <- sum(Results$selected[,x+colshift] == 1, na.rm = TRUE)
        sum_y <- sum(Results$selected[,y+colshift] == 1, na.rm = TRUE)
        sum_xy <- sum(Results$selected[,x+colshift] == 1 & 
                        Results$selected[,y+colshift] == 1, na.rm = TRUE)
        freq_x <- sum_x / nrow(Results$selected)
        freq_y <- sum_y / nrow(Results$selected)
        freq_xy <- sum_xy / nrow(Results$selected)
        
        h_xy <- log2(1/freq_xy)
        
        Results$pmi[x,y] <- log2(freq_xy / (freq_x * freq_y))
        Results$npmi[x,y] <- Results$pmi[x,y] / h_xy
      }
    }
  }
  
  return (Results)
  }