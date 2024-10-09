covid=read.table("C:/Edinburgh/Studying/Extended Statistical Programming/engcov.txt",header=TRUE)
covid = covid[1:150,]

#deconv = function(t,deaths,n.rep=100,bs=FALSE,t0=NULL){
distribution = dlnorm(c(1:80),meanlog = 3.152, sdlog = 0.451)
distribution = distribution/sum(distribution)
days = length(covid$nhs)
deaths = covid$nhs
death_day = rep(c(1:days), times = deaths)
n = length(death_day)
samps = sample(x = c(1:80),prob = distribution, size = n, replace = TRUE)
t0 = death_day - samps
t0[t0<1] = 1 #no one is infected before day 1

n.rep = 100
inft = matrix(rep(NA,n.rep*days),nrow = days)
P = c(rep(NA,n.rep))
samps2 = sample(x = c(1:80),prob = distribution, size = n, replace = TRUE)
t0 = t0+samps2 #initial death days
t0[t0<1]=1
t0[t0>150]=150

compute_p = function(t0){
  death_sim = tabulate(t0,nbins = days)
  copy = death_sim
  copy = copy[copy<1]=1 #according to formula, denominator cannot be 0
  p = sum(((covid$nhs - death_sim)^2)/copy)
}
k = 0
for (i in 1:n.rep){

if (i==1) p_before = compute_p(t0) else p_before = P[i-1]

if (i<=50){
  samp = c(-8,-4,-2,-1,1,2,4,8)
} else if (i<=75){
  samp = c(-4,-2,-1,1,2,4)
} else samp = c(-2,-1,1,2)
t1 = c(rep(NA,n))
temp = sample(1:n, size = n) #all indexes in random order
c = sample(samp, size = n, replace = TRUE)
t1[temp] = t0[temp] + c #pushing back or adding days for death
t1[t1<1]=1
t1[t1>150]=150
#t1 = t1+61
death_sim2 = tabulate(t1,nbins = 150)

for (j in 1:n){
  ind = temp[j]
  add = death_sim2[t1[ind]]
  sub = death_sim2[t0[ind]]
  new_val = ((covid$nhs[t1[ind]]-add)^2)/max(1,add) + ((covid$nhs[t0[ind]]-sub)^2)/max(1,sub)
  old_val = ((covid$nhs[t0[ind]]-(sub+1))^2)/max(1,sub+1) + ((covid$nhs[t1[ind]]-(add-1))^2)/max(1,add-1)
  if (new_val<old_val) {
    t0[ind] = t1[ind]
  } else {
    death_sim2[t0[ind]] = sub+1
    death_sim2[t1[ind]] = add-1
  }
  }

p_after = compute_p(t0)
inft[,i] = tabulate(t0,nbins = days)
P[i] = min(p_before,p_after)
# PLOT STUFF
# CHANGE THE C
# ADD THE BOOTSTRAP
}
#return (list(P, inft, t0))
#}