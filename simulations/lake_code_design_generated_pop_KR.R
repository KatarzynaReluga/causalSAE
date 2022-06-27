rm(list = ls(all = TRUE))


library(MASS)
library(nlme)
library(sp)
#install.packages("spgwr")
#Geographically weighted regression
library(spgwr)

set.seed(1973)

lake.s = read.table("lake.north_fin.txt", header = TRUE, sep = "\t")
lake.frame = read.table("lake_pop_rand.txt", header = TRUE, sep = "\t")

N = nrow(lake.frame)
n = nrow(lake.s)

ANC = NULL
ELEV = NULL
for (i in 1:N)
{
  euc_dist = NULL
  a = which(lake.s[, 3] > 1)
  mat_dist = rbind(as.matrix(lake.frame[i, 5:6]), as.matrix(lake.s[a, 4:5]))
  euc_dist = c(as.matrix(dist(mat_dist))[-1, 1])
  aa = which(euc_dist == 0)
  euc_dist[aa] = 0.00001
  #weight=1/(euc_dist^2)
  weight = exp((-euc_dist ^ 2) / (2.91 ^ 2))
  t = sample(a, 1, prob = weight)
  lake.s[t, 3] = lake.s[t, 3] - 1
  ANC[i] = lake.s[t, 2]
  ELEV[i] = lake.s[t, 1]
  rm(weight)
  rm(euc_dist)
  print(i)
}

pop = cbind(lake.frame, ANC, ELEV)
ni = table(lake.s[, 5])
id.sample = matrix(0, 551, 500)
for (i in 1:500)
{
  s = NULL
  for (j in 1:86)
  {
    s1 <- sample(pop[, 1][pop[, 2] == ar[j]], ni[j])
    s = c(s, s1)
  }
  
  id.sample[, i] = s
  
}

write.table(id.sample, "StratSample551_ray.txt", sep = "\t")
write.table(pop, "PopSp.txt", sep = "\t")
