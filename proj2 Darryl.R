covid=read.table("C:/Edinburgh/Studying/Extended Statistical Programming/engcov.txt",header=TRUE)
covid = covid[1:150,]

deconv = function(t,deaths,n.rep=100,bs=FALSE,t0=NULL){
t = covid$julian
deaths = covid$nhs
distribution = dlnorm(c(1:80),meanlog = 3.152, sdlog = 0.451)
distribution = distribution/sum(distribution)
len = length(t)
days = 310
death_day = rep(c((1+61):(len+61)), times = deaths) #results are starting from 62nd day
n = length(death_day)
samps = sample(x = c(1:80),prob = distribution, size = n, replace = TRUE)
t0 = death_day - samps #initial incidence days from real death days
t0[t0<1] = 1 #no one is infected before day 1
death_plot = c(rep(0,days)) 
death_plot[covid$julian] = covid$nhs 
plot(death_plot, type = 'l', col = 'orange', xlab = 'Days of the year', 
     ylab = 'No. of deaths') #real deaths per day

n.rep = 100
inft = matrix(rep(0,n.rep*days),nrow = days)
P = c(rep(0,n.rep))

compute_p = function(t){
  death_sim = tabulate(t,nbins = days)
  copy = death_sim
  copy[copy<1]=1 #according to formula, denominator cannot be 0
  p = sum(((death_plot - death_sim)^2)/copy)
}

for (i in 1:n.rep){

if (i<=50){
  samp = c(-8,-4,-2,-1,1,2,4,8)
} else if (i<=75){
  samp = c(-4,-2,-1,1,2,4)
} else samp = c(-2,-1,1,2)

samps2 = sample(x = c(1:80),prob = distribution, size = n, replace = TRUE)
td = t0+samps2 #initial death days = initial incidence + incidence-to-death distribution data
td[td<1]=1

#keep track of values of p for comparison later
if (i==1) p_before = compute_p(td) else {p_before = P[i-1]
td = tf}

#tf is the vector of incidence days for optimal p
t1 = tf = c(rep(NA,n))
temp = sample(1:n, size = n) #all indexes in random order
c = sample(samp, size = n, replace = TRUE)
t1 = td + c #pushing back or adding days for death
t1[t1<1]=1
death_sim2 = tabulate(t1,nbins = days) #new deaths per day vector

for (j in 1:n){
  ind = temp[j] #index to consider
  add = death_sim2[t1[ind]] #num of deaths on date being shifted to
  sub = death_sim2[td[ind]] #num of deaths on date being shifted from (original date of death)
  
  #comparing only both dates that have their death numbers changed because other date values stay the same
  new_val = ((death_plot[t1[ind]]-add)^2)/max(1,add) + ((death_plot[td[ind]]-sub)^2)/max(1,sub)
  old_val = ((death_plot[td[ind]]-(sub+1))^2)/max(1,sub+1) + ((death_plot[t1[ind]]-(add-1))^2)/max(1,add-1)
  
  if (new_val<old_val) { #adding c is gives a lower(better) value of P
    tf[ind] = t1[ind] #tf contains optimal death days for each person each iteration
  } else {
    tf[ind] = td[ind]
    death_sim2[td[ind]] = sub+1 #updating deaths per day vector
    death_sim2[t1[ind]] = add-1
  }
  }

p_after = compute_p(tf) #plot simulated deaths per day
inc = tf-samps2[temp] #simulated incidence vector
lines(tabulate(tf,nbins = days), col = 'blue')
lines(tabulate(inc,nbins = days), col = 'red')
inft[,i] = tabulate(tf,nbins = days)
P[i] = p_after
#P[i] = min(p_before,p_after)


# PLOT STUFF
#lines(inft[,i], col = red)

}
#lines(inft[,n.rep], col = 'blue')
# ADD THE BOOTSTRAP
print(length(inc))
if (bs == TRUE){
  plot(inft[,100], type = 'l', col = 'orange')
  death_pois = rpois(len,covid$nhs[i])
  death_plot = c(rep(0,days))
  death_plot[covid$julian] = death_pois
  
  inft = matrix(rep(0,n.rep*days),nrow = days)
  P = c(rep(0,n.rep))
  for (k in 1:len){
    print(k)
    
    
 #   for (ii in 1:n.rep){
      if (k<=50){
        samp = c(-8,-4,-2,-1,1,2,4,8)
      } else if (k<=75){
        samp = c(-4,-2,-1,1,2,4)
      } else samp = c(-2,-1,1,2)
      
      #samps2 = sample(x = c(1:80),prob = distribution, size = n, replace = TRUE)
      #td = tf+samps2 #initial death days = initial incidence + incidence-to-death distribution data
      #td[td<1]=1
      
      #keep track of values of p for comparison later
      if (k==1) p_before = compute_p(td) else p_before = P[k-1]
      td = tf
      
      t1 = tf = c(rep(NA,n))
      temp = sample(1:n, size = n) #all indexes in random order
      c = sample(samp, size = n, replace = TRUE)
      t1 = td + c #pushing back or adding days for death
      t1[t1<1]=1
      death_sim2 = tabulate(t1,nbins = days) #new deaths per day vector
      
      for (jj in 1:n){
        ind = temp[jj] #index to consider
        add = death_sim2[t1[ind]] #num of deaths on date being shifted to
        sub = death_sim2[td[ind]] #num of deaths on date before shifting
        
        #comparing only both dates that have their death numbers changed because other date values stay the same
        new_val = ((death_plot[t1[ind]]-add)^2)/max(1,add) + ((death_plot[td[ind]]-sub)^2)/max(1,sub)
        old_val = ((death_plot[td[ind]]-(sub+1))^2)/max(1,sub+1) + ((death_plot[t1[ind]]-(add-1))^2)/max(1,add-1)
        
        if (new_val<old_val) { #adding c is gives a lower(better) value of P
          tf[ind] = t1[ind] #td contains optimal death days for each person each iteration
        } else {
          tf[ind] = td[ind]
          death_sim2[td[ind]] = sub+1 #updating deaths per day vector
          death_sim2[t1[ind]] = add-1
        }
      }
      lines(tabulate(tf,nbins = days),col = 'blue')
  }
  
  
}
#return (list(P, inft, t0))
}
#}

deconv(covid$julian,covid$nhs,bs = FALSE)
