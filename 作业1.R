setwd("C:/Users/HUAWEI/Desktop/sds/extend statisting programming") ## comment out of submitted
a <- scan("4300-0.txt",what="character",skip=73,nlines=32858-73,
          fileEncoding="UTF-8")
a <- gsub("_(","",a,fixed=TRUE) ## remove "_("
#ba <- a

#4-1
split_punct <- function(words,pwords){
  new_words <- character(0)
  for (i in 1:length(words)) {
    indices <- grep(pwords[i],words,fixed = TRUE)
    no_p_words <- gsub(pwords[i],'',words,fixed = TRUE)
    #punct <- gsub("[^[:punct:]]",'',words)
    punct <- gsub(paste0("[^", pwords[i], "]"), "", words)
    new_words <- rep('0',length(words)+length(punct))
    p_indices <- indices + (1:length(punct))
    new_words[p_indices] <- punct[nzchar(punct)]
    new_words[-p_indices] <- no_p_words
    words <- new_words
  }
  return(new_words)
}


#4-2
split_punct <- function(words, pwords) {
  new_words <- c() 
  for (i in 1:length(words)) {
    no_p_word <- gsub(pwords, "", words[i], fixed = TRUE)
    new_words <- c(new_words, no_p_word)
    punct <- gsub("[^[:punct:]]", "", words[i])
    new_words <- c(new_words, punct)
  }
  new_words <- new_words[new_words != ""]
  return(new_words)
}

pwords <- c(",", ".", ";", "!", ":", "?")
a <- split_punct(a,pwords)

#5
pwords <- c(",", ".", ";", "!", ":", "?")
a <- split_punct(a,pwords)

#6
l_a <- tolower(a)
b <- unique(l_a)
vector_indices <- match(l_a, b)
words_number <- tabulate(vector_indices)
number <- sort(words_number, decreasing=TRUE)[1000]
b <- b[words_number >= number]

#7
indices <- match(l_a, b)
mlag <- 4 
n <- length(indices)
M <- matrix(NA, nrow=(n-mlag), ncol=(mlag+1))

for (i in 1:(mlag+1)) {
  M[,i] <- indices[i:(n-mlag+i-1)]
}
#8
nw <- 50  

text <- vector("character", nw)
text[1] <- sample(b, 1) 

for (i in 2:nw) {
  for (j in mlag:1) {
    if (i > j) {
      current_indices <- M[i-j, ]
      possible_words <- b[!is.na(current_indices)]#返回 current_indices 中不是 NA 的位置。
      if (length(possible_words) > 0) {
        text[i] <- sample(possible_words, 1)
        break
      }
    }
  }
}

cat(paste(generated_text, collapse=" "))

#9
sections <- sample(b, nw, replace=TRUE, prob=words_number[1:length(b)])
cat(paste(sections, collapse=" "))



