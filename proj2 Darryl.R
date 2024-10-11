covid=read.table("C:/Edinburgh/Studying/Extended Statistical Programming/engcov.txt",header=TRUE)
covid = covid[1:150,]

#deconv = function(t,deaths,n.rep=100,bs=FALSE,t0=NULL){
distribution = dlnorm(c(1:80),meanlog = 3.152, sdlog = 0.451)
distribution = distribution/sum(distribution)
len = length(covid$nhs)
days = 310
deaths = covid$nhs
death_day = rep(c((1+61):(len+61)), times = deaths)
n = length(death_day)
samps = sample(x = c(1:80),prob = distribution, size = n, replace = TRUE)
t0 = death_day - samps #initial incidence days
t0[t0<1] = 1 #no one is infected before day 1
death_plot = c(rep(0,days))
death_plot[covid$julian] = covid$nhs
plot(death_plot, type = 'l', col = 'orange') #real deaths per day

n.rep = 100
inft = matrix(rep(0,n.rep*days),nrow = days)
P = c(rep(0,n.rep))

compute_p = function(t, plot = FALSE){
  death_sim = tabulate(t,nbins = days)
  if (plot == TRUE) {
    lines(death_sim, col = 'blue')
  }
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
td = t0+samps2 #initial death days
td[td<1]=1

if (i==1) p_before = compute_p(td) else p_before = P[i-1]

t1 = c(rep(NA,n))
temp = sample(1:n, size = n) #all indexes in random order
c = sample(samp, size = n, replace = TRUE)
t1[temp] = td + c #pushing back or adding days for death
t1[t1<1]=1
death_sim2 = tabulate(t1,nbins = days)

for (j in 1:n){
  ind = temp[j]
  add = death_sim2[t1[ind]]
  sub = death_sim2[td[ind]]
  new_val = ((death_plot[t1[ind]]-add)^2)/max(1,add) + ((death_plot[td[ind]]-sub)^2)/max(1,sub)
  old_val = ((death_plot[td[ind]]-(sub+1))^2)/max(1,sub+1) + ((death_plot[t1[ind]]-(add-1))^2)/max(1,add-1)
  if (new_val<old_val) {
    td[ind] = t1[ind]
  } else {
    death_sim2[td[ind]] = sub+1
    death_sim2[t1[ind]] = add-1
  }
  }

p_after = compute_p(td,TRUE)# plot simulated deaths per day
inc = td-samps2
inc = inc[inc>0]
lines(tabulate(inc,nbins = len), col = 'red')# plot incidence per day
print(i)
inft[,i] = tabulate(td,nbins = days)
P[i] = min(p_before,p_after)


# PLOT STUFF
#lines(inft[,i], col = red)

}
# ADD THE BOOTSTRAP
if (bs == TRUE){
  for (i in 1:len){
    rpois(len,P[n.rep,i])
    
  }
}
#return (list(P, inft, t0))
#}

plot(inft[,100], type = 'l', col = 'blue')
