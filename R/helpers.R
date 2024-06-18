

mf<-function(x,digits=2){
  s<-format(round(x,digits=digits),
            digits=digits,drop0trailing=T,scientific=F,nsmall=digits)
  s<-sub("(?<![0-9])0+(?=\\.)", "",s, perl = TRUE)
  return(s)
}

