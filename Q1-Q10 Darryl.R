setwd("C:/Edinburgh/Studying/Extended Statistical Programming") ## comment out of submitted
a <- scan("4300-0.txt",what="character",skip=73,nlines=32858-73,
          fileEncoding="UTF-8")
a <- gsub("_(","",a,fixed=TRUE) ## remove "_("

txt = c("Inside","skeleton,", "of", "this", "doubt,", "but", "waiting", "out.")
p = '[,.:;!?]$'
# QUESTION 5
split_punct = function(ss,p){
  #ss = strsplit(s," ")[[1]]
  ns = grep(p,ss)
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
#split_punct(txt,p)
a = split_punct(a,p)

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
nw = 50
sentence = rep(".",nw) #initialise a sentence vector
word_index = 2 #index of first predicted word
w = rep(NA,mlag-1) #initialise window vector of length mlag-1, first word will be appended later to make length mlag
first_word_index = sample(length(b),1) #randomly sample first word
sentence[1] = b[first_word_index] #add it to the sentence vector
cat(sentence[1], sep = '\n')
w = append(w,first_word_index)
w = as.double(w) #ensure that elements in w are the same type as m
for (i in 2:nw) {
  for (j in mlag:1) if (i>j) { ## skip lags too long for current i
    possible_words = c() #initialise vector of possible words to sample from
    for (k in 1:nrow(m)){
      if (identical(m[k,(mlag-j+1):mlag], w[(mlag-j+1):mlag])){ #check if last j words of w are the same as in m[k,]
        if (!is.na(m[k,(mlag + 1)])){ 
          possible_words = append(possible_words,m[k,(mlag + 1)]) #add words that are not NA to possible_words
        }
      }
    }
    if (length(possible_words)>1){
      x = possible_words[sample(length(possible_words),1)]
      sentence[word_index] = b[x]
      word_index = word_index + 1
      w = append(w[2:mlag],x) #add x to last position of w, remove first element
      cat(b[x], sep = '\n')
      break
    }
  }
}

# QUESTION 9
s = sample(length(b),50, replace = TRUE) #not sure if sample from b or whole text
b[s]

# QUESTION 10
sz = 1000
b = paste(toupper(substr(b[1:sz],1,1)),substr(b[1:sz],2,nchar(b[1:sz])), sep = '')
