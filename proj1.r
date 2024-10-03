# Jiayi Wu-s2664441, Ruobing Dai-s2655029, Darryl Chia-s2740198
# Darryl Chia created github repositories and write question 8(the most difficult one), Ruobing Dai did question 2 to 5 and Jiayi Wu completed question 6,7,9,10. Finally, both of us check and modify our work.

# setwd("C:/Users/HUAWEI/Desktop/sds/extend statisting programming") ## comment out of submitted
a <- scan("4300-0.txt",what="character",skip=73,nlines=32858-73,
          fileEncoding="UTF-8")
a <- gsub("_(","",a,fixed=TRUE) ## remove "_("


# 4
split_punct <- function(words,pwords){ #words is word text, pwords is punctuations to split
  ps = grep(pwords,words) #get indexes of words with punctuations
  l = length(words)
  lps = length(ps)
  xs <- rep("",l+lps)
  iis <- ps+0:(lps-1) #pushing back indexes to make space for punctuations
  iis2 = ps+1:lps #indexes for punctuations
  xs[iis2] <- substr(words[ps],nchar(words[ps]),nchar(words[ps])) #insert last character of word, should be punctuation
  xs[iis] = substr(words[ps],1,nchar(words[ps])-1) #insert the word without punctuation
  xs[-append(iis,iis2)] = substr(words[-ps], 1,100) #insert the remaining words
  xs
}


# 5
pwords <- '[,.:;!?]$'
a <- split_punct(a,pwords)


# 6
txt <- a[nzchar(a) & !is.na(a)] # remove NA and "" from text
uni_txt <- unique(tolower(txt)) # find unique words
ii <- match(tolower(txt),uni_txt) # mark words in main text with unique words index
freq <- tabulate(ii, nbins = length(uni_txt)) # calculate word frequency
desc <- sort(freq, decreasing = TRUE) #sort words' frequency in descending order
# Since more than one words have same frequency of 1000th common words, it is concise to select all the words with same frequency of 1000th common words
threshold <- desc[1000] # get the threshold frequency of common words
b <- uni_txt[which(freq >= threshold)] # select m most common words, m is more than 1000


# 7
token <- match(tolower(txt),b) #replace the common words in main text with its token
mlag <- 4
M <- matrix(rep(0,(length(token)-mlag)*(mlag+1)), length(token)-mlag, mlag+1) # initialise M matrix
for (i in 1:(mlag+1)) {
  M[,i] <- token[i:(length(token)-mlag-1+i)]
}


# 8
nw <- 50
sentence <- rep(NA,nw) #initialise a sentence vector
word_index <- 2 #index of first predicted word
w <- rep(NA,mlag-1) #initialise window vector of length mlag-1, first word will be appended later to make length mlag
first_word_index <- sample(length(b),1) #randomly sample first word
sentence[1] <- b[first_word_index] #add it to the sentence vector
cat(sentence[1], sep = '\n')
w <- append(w,first_word_index)
w <- as.double(w) #ensure that elements in w are the same type as M
for (i in 2:nw) {
  for (j in mlag:1) if (i>j) { ## skip lags too long for current i
    possible_words <- c() #initialise vector of possible words to sample from
    for (k in 1:nrow(M)){
      if (identical(M[k,(mlag-j+1):mlag], w[(mlag-j+1):mlag])){ #check if last j words of w are the same as in M[k,]
        if (!is.na(M[k,(mlag + 1)])){ 
          possible_words <- append(possible_words,M[k,(mlag + 1)]) #add words that are not NA to possible_words
        }
      }
    }
    if (length(possible_words)>1){
      x <- possible_words[sample(length(possible_words),1)]
      sentence[word_index] <- b[x]
      word_index <- word_index + 1
      w <- append(w[2:mlag],x) #add x to last position of w, remove first element
      cat(b[x], sep = '\n')
      break
    }
  }
}


# 9
n = 50
token_select <- sample(token[!is.na(token)], n, replace = TRUE) # sample tokens without returning
word_select <- b[token_select] # find certain word to each token
cat(word_select, sep = '\n')


# 10
txt <- spl_punc[nzchar(spl_punc) & !is.na(spl_punc)]
position_upper <- grep('[A-Z]', txt) # find the index of words start with capital letter in text
uni_txt <- unique(tolower(txt))  #extract unique text
ii_lower <- match(txt[-position_upper],uni_txt) #extract tokens of lower case words
ii_upper <- match(tolower(txt[position_upper]),uni_txt) #extract tokens of words start with capital letter
u_seq <- tabulate(ii_lower, nbins = length(uni_txt)) > tabulate(ii_upper,nbins = length(uni_txt)) #find the words index which most often start with a capital letter in the main text
uni_txt[which(!u_seq)] <- paste(toupper(substr(uni_txt[which(!u_seq)],1,1)),substr(uni_txt[which(!u_seq)],2,nchar(uni_txt[which(!u_seq)])), sep = '') # replace these words with words start with a capital letter
ii <- match(tolower(txt),tolower(uni_txt))
freq <- tabulate(ii, nbins = length(uni_txt))
desc <- sort(freq, decreasing = TRUE)
threshold <- desc[1000]
b <- uni_txt[which(freq >= threshold)] # m most common words with certain words start with capital letter
## go back to question 7
token <- match(txt,b) #replace the common words in main text with its token
mlag <- 4
M <- matrix(rep(0,(length(token)-mlag)*(mlag+1)), length(token)-mlag, mlag+1) # initialise M matrix
for (i in 1:(mlag+1)) {
  M[,i] <- token[i:(length(token)-mlag-1+i)]
}
## Then we can go back to question 8 to simulate again.
