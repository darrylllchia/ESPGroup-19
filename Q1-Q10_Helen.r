#3.
setwd("/Users/macbook/Statistical-Programming") ## comment out of submitted
a <- scan("4300-0.txt",what="character",skip=73,nlines=32858-73,
fileEncoding="UTF-8")
a <- gsub("_(","",a,fixed=TRUE) ## remove "_("
length(a)

#4.
#######Q:if there are 2 punctions in one character??
# txt: a character vector where matches are sought
# + or an object which can be coerced by ‘as.character’ to a character vector. Long vectors are supported.
# punc:a character vector contains punctuation mark(s)(one or more)
split_punct <- function(txt,punc){
    punc_len <- length(punc)
    text <- character(0)
    for (i in 1:punc_len) {
        pos <- grep(punc[i],txt,fixed = TRUE)
        words <- gsub(punc[i],'',txt,fixed = TRUE)
        punction <- gsub(paste0("[^", punc[i], "]"), "", txt)
        text <- rep('0',length(txt)+length(pos))
        pos_punc <- pos+(1:length(pos))
        text[pos_punc] <- punction[nzchar(punction)]
        text[-pos_punc] <- words
        txt <- text
    }
    return(text)
}

# 5.
spl_punc <- split_punct(a,c(",",".", ";", "!", ":","?")) ##over

#6.
#a
txt <- spl_punc[nzchar(spl_punc) & !is.na(spl_punc)] #remove NA and "" from text
uni_txt <- unique(tolower(txt)) # find the vector of unique words
#b
ii <- match(tolower(txt),uni_txt)
#c
freq <- tabulate(ii, nbins = length(uni_txt))
#d
desc <- sort(freq, decreasing = TRUE)
# l_range <- match(desc[1000]+1,desc)
# h_range <- match(desc[1000]-1,desc)-1
threshold <- desc[1000]
#e
b <- uni_txt[which(freq >= threshold)] ##m most common words


# 7.
# a
token <- match(tolower(txt),b)
# b
mlag <- 4
M <- matrix(rep(0,(length(token)-mlag)*(mlag+1)), length(token)-mlag, mlag+1)
for (i in 1:(mlag+1)) {
    M[,i] <- token[i:(length(token)-mlag-1+i)]
}


# 8.
mlag <- 4
nw <- 50
token_sel <- rep(0,nw)
token1 <- sample(token[!is.na(token)], 1)
token_seq <- c(rep(NA,mlag-1),token1)
token_sel[1] <- token1 #1st word
for (i in 2:nw) {
    for (j in mlag:1) if (i>j) { ## skip lags too long for current i
        valid_match_row <- c()
        for (k in 1:nrow(M)){
            if(setequal(M[k, 1:j],token_seq[1:j]) == TRUE && is.na(M[k,j+1]) == FALSE){
                valid_match_row[k] <- k} 
            else {valid_match_row[k] <- NA}
    }
    if (sum(!is.na(valid_match_row)) > 1) break
}
token_sel[i] <- sample(M[which(!is.na(valid_match_row)),j+1],1)
token_seq <- c(token_seq[2:length(token_seq)], token_sel[i])
}
words_sel <- b[token_sel]
cat(words_sel, sep = '\n')


# 9.
n = 50
b_fac <- as.factor(b)
b_num <- as.numeric(b_fac)
uni_b_num <- unique(b_num)
freq <- tabulate(b_num)
prob <- freq / sum(freq)
b_index <- sample(uni_b_num,n,replace = FALSE,prob = prob)
word_select <- b[b_index]
cat(word_select, sep = '\n')


# 10.
txt <- spl_punc[nzchar(spl_punc) & !is.na(spl_punc)]
position_upper <- grep('[A-Z]', txt)
uni_txt <- unique(tolower(txt))  #unique text 不会把标点符号unique
ii_lower <- match(txt[-position_upper],uni_txt)
ii_upper <- match(tolower(txt[position_upper]),uni_txt)
u_seq <- tabulate(ii_lower, nbins = length(uni_txt)) > tabulate(ii_upper,nbins = length(uni_txt))
uni_txt[which(!u_seq)] <- paste(toupper(substr(uni_txt[which(!u_seq)],1,1)),substr(uni_txt[which(!u_seq)],2,nchar(uni_txt[which(!u_seq)])), sep = '')
ii <- match(tolower(txt),tolower(uni_txt))
freq <- tabulate(ii, nbins = length(uni_txt))
desc <- sort(freq, decreasing = TRUE)
threshold <- desc[1000]
b <- uni_txt[which(freq >= threshold)] ##m most common words
# ? space before punction? we do not have?