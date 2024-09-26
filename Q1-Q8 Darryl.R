setwd("C:/Edinburgh/Studying/Extended Statistical Programming") ## comment out of submitted
a <- scan("4300-0.txt",what="character",skip=73,nlines=32858-73,
          fileEncoding="UTF-8")
a <- gsub("_(","",a,fixed=TRUE) ## remove "_("
length(a)

s = "An omnishambles, in a headless chicken factory."
split_punct = function(ss){
  #ss = strsplit(s," ")[[1]]
  ns = grep('[,.:;!?]$',ss)
  l = length(ss)
  lns = length(ns)
  xs <- rep("",l+lns)
  iis <- ns+0:(lns-1)
  iis2 = ns+1:lns
  xs[iis2] <- substr(ss[ns],nchar(ss[ns]),nchar(ss[ns]))
  xs[iis] = substr(ss[ns],1,nchar(ss[ns])-1)
  xs[-append(iis,iis2)] = substr(ss[-ns], 1,100)
  xs
}
a = split_punct(a)

# QUESTION 6
b = unique(tolower(a))
c = match(tolower(a),b)
tab = tabulate(c, nbins = length(c))
d = sort(tab, decreasing = TRUE, index.return = TRUE)$ix[1:1000]
b = b[d] #m most common words

# QUESTION 7
indices = match(tolower(a),b)
mlag = 4
m = matrix(rep(0,(length(a)-mlag)*(mlag+1)),nrow = length(a)-mlag)
for (i in (1:(length(indices)-mlag))){
  m[i,] = indices[i:(i+mlag)]
}

# QUESTION 8
nw = 5
sentence = rep("",nw)
word_index = 2
w = rep(NA,mlag-1)
first_word_index = sample(length(b),1)
sentence[1] = b[first_word_index]
w = append(w,first_word_index)
for (i in 2:nw) {
  for (j in mlag:1) if (i>j) { ## skip lags too long for current i
    possible_words = c()
    th_word = mlag - sum(is.na(w)) + 1
    nth = min(mlag,th_word) - 1
    for (k in 1:nrow(m)){
      if (setequal(m[k,(mlag-j):mlag], w[(mlag-j):mlag])){
        if (!is.na(m[k,(mlag + 1)])){
          possible_words = append(possible_words,m[k,(mlag + 1)])
        }
      }
    }
    if (length(possible_words)>0){
      x = possible_words[sample(length(possible_words),1)]
      sentence[word_index] = b[x]
      word_index = word_index + 1
      w = append(w[2:mlag],possible_words[x])
      cat(b[x], sep = '\n')
      break
    }
  }
}
