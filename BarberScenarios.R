library(R2OpenBUGS)
library(coda)
library(BRugs)
workdir <- "/home/sarah/Downloads/Milltown Offsets Review"

setwd(workdir)

#1) Segregated Habitat 
#Alewife - 254 acres
#Blueback Herring - 913 acres
#Carrying capacity - Limnotech report - 846 river herring per acre

#Using Barber Model 
#Number of adults - carrying capacity for herring x segregated habitat for Alewife and Blueback Herring
#Alewife
254*846 #214884 adults
#Blueback Herring
913*846 #772398 adults

#juveniles share resources across all watersheds
#80% Alewife and 20% Blueback Herring
#214884/(214884 + x) = 0.8
#80% alewife = 214884
#20 % blueback herring = 53721

#for the entire watershed
#987282 for Alewife
#80% Alewife and 20% Blueback Herring
#987282/(987282 + x) = 0.8
#80% alewife = 987282
#20% blueback herring = 2426827
#Model for alewife
#asymptote for juveniles = 9249 per acre
#CI from Devine paper 5240, 11663
#SD = sqrt(29)*((11663-5240)/3.92) = 8824
#Variance = (8824)^2 = 77877682
#Precision = 1/(77857682) = 1.28e-8


#juvenile habitat for all species is 1167
#adding an extra row to the dataset for juvenile habitat
cat("
model { 
  for (i in 1:3) {
  Egg[i] ~ dnorm(mu[i], tau)
  mu[i] <- ((Adults[i]*exp(-0.25*Mspawn)*0.5)*fecundity)
    Juveniles[i] ~ dnorm(muj[i], tauj) 
    muj[i] <- ((alpha*(mu[i]))/(1 + ((alpha*mu[i])/recruit_asymp)))
  }
fecundity <- (b*w) - c
b <- 871.72 
c <- 50.916
w <- 144
recruit_asymp ~ dnorm(9249, 11)
#recruit_asymp ~ dnorm(3283, 11)
Mspawn ~ dnorm(2.391,1)
#Mspawn <- 2.391
tau ~ dgamma(1000, 0.1)
tauj ~ dgamma (500, 0.001)
alpha ~ dnorm(0.0015, 1)
}
", file = "BarberModelEggsNoIE.txt")
modelCheck("BarberModelEggsNoIE.txt")
#### Data  ##########
cat(
"Adults,Egg,Juveniles
486500,16788683128,3282.572
214884,NA,NA
987282,NA,NA
", file="Data2.csv")
Data = read.csv("Data2.csv",header=TRUE)
attach(Data)

bug.dat=list("Adults", "Egg", "Juveniles")

init.fun=function(){
   recruit_asymp = rnorm(1, 9000, 11)
   alpha = rnorm(1,0.0015,0.0002)
   Mspawn=rnorm(1,2.391,0.01)
   tau = rgamma(1, 10000, 0.001)
   mu = rnorm(1, 10000, 1)
   #mu = rpois(1,1000000)
   #muj = rpois(1,1000)
   muj = rnorm(1, 1000, 10)
   tauj = rgamma(1, 1000, 0.01)
   Egg = rnorm(1, 1000000, 1000)
   Juveniles = rnorm(1, 1000, 0.1)
   return(list(Egg = Egg, Juveniles = Juveniles, mu = mu, 
  alpha = alpha, muj = muj
  ))}
        

params=c("mu", "muj", "recruit_asymp")
Sys.time()

EggBug=bugs(data = bug.dat, inits = init.fun, parameters.to.save = params, model.file="BarberModelEggsNoIE.txt",
    n.chains=3, n.iter=100000, n.burnin=50000,
    n.thin=100, codaPkg = TRUE
    )
Sys.time()
(EggBug)
plot(EggBug)

#checking gelman rubin diagnostic
eggbugmc<-as.mcmc.list(EggBug)
gelman.diag(eggbugmc)
#checking convergence using Coda
codaegg <-read.bugs(EggBug)
plot(codaegg)
gelman.plot(codaegg)
gelman.diag(codaegg)


summary(codaegg)

### same procedure for Blueback Herring
#except that Blueback herring adults are mzximized to be 20% of Alewife population

cat("
model { 
  for (i in 1:3) {
  Egg[i] ~ dnorm(mu[i], tau)
  mu[i] <- ((Adults[i]*exp(-0.25*Mspawn)*0.5)*fecundity)
    Juveniles[i] ~ dnorm(muj[i], tauj) 
    muj[i] <- ((alpha*(mu[i]))/(1 + ((alpha*mu[i])/recruit_asymp)))
  }
fecundity <- (b*w) - c
b <- 871.72 
c <- 50.916
w <- 144
recruit_asymp ~ dnorm(9249, 11)
#recruit_asymp ~ dnorm(3283, 11)
Mspawn ~ dnorm(2.391,1)
#Mspawn <- 2.391
tau ~ dgamma(1000, 0.1)
tauj ~ dgamma (500, 0.001)
alpha ~ dnorm(0.0015, 1)
}
", file = "BarberModelEggsNoIE.txt")
modelCheck("BarberModelEggsNoIE.txt")
#### Data  ##########
cat(
"Adults,Egg,Juveniles
486500,16788683128,3282.572
53721,NA,NA
246827,NA,NA
", file="Data2.csv")
Data = read.csv("Data2.csv",header=TRUE)
attach(Data)

bug.dat=list("Adults", "Egg", "Juveniles")

init.fun=function(){
   recruit_asymp = rnorm(1, 9000, 11)
   alpha = rnorm(1,0.0015,0.0002)
   Mspawn=rnorm(1,2.391,0.01)
   tau = rgamma(1, 10000, 0.001)
   mu = rnorm(1, 10000, 1)
   #mu = rpois(1,1000000)
   #muj = rpois(1,1000)
   muj = rnorm(1, 1000, 10)
   tauj = rgamma(1, 1000, 0.01)
   Egg = rnorm(1, 1000000, 1000)
   Juveniles = rnorm(1, 1000, 0.1)
   return(list(Egg = Egg, Juveniles = Juveniles, mu = mu, 
  alpha = alpha, muj = muj
  ))}
        

params=c("mu", "muj", "recruit_asymp")
Sys.time()

EggBug=bugs(data = bug.dat, inits = init.fun, parameters.to.save = params, model.file="BarberModelEggsNoIE.txt",
    n.chains=3, n.iter=100000, n.burnin=50000,
    n.thin=100, codaPkg = TRUE
    )
Sys.time()

codaegg <-read.bugs(EggBug)
plot(codaegg)
gelman.plot(codaegg)
gelman.diag(codaegg)


summary(codaegg)

#updating Scenario 1 
#adding stochastic carrying capacity 
#adding CIs from Devine paper to the asymptote with gamma
#precision for adult carrying capacity 
#(1186-672)/1.284 = 400 
#Variance = 400^2 = 1600
#precision = 1/1600 = 0.0006

cat("
model { 
  for (i in 1:3) {
  Adults[i] ~ dnorm(mua[i], 0.0006)
  mua[i] <- 846*Area[i]
  Egg[i] ~ dnorm(mu[i], tau)
  mu[i] <- ((Adults[i]*exp(-0.25*Mspawn)*0.5)*fecundity)
    Juveniles[i] ~ dnorm(muj[i], tauj) 
    muj[i] <- ((alpha*(mu[i]))/(1 + ((alpha*mu[i])/recruit_asymp)))

  }
fecundity <- (b*w) - c
b <- 871.72 
c <- 50.916
w <- 144
recruit_asymp ~ dnorm(9249, 0.01)
#recruit_asymp ~ dnorm(3283, 11)
Mspawn ~ dnorm(2.391,1)
#Mspawn <- 2.391
tau ~ dgamma(100000, 0.00000595)
tauj ~ dgamma (500, 0.001)
alpha ~ dnorm(0.0015, 1)
}
", file = "BarberModelwithArea.txt")
modelCheck("BarberModelwithArea.txt")
#### Data  ##########
cat(
"Area,Adults,Egg,Juveniles
1167,486500,16788683128,3282.572
254,NA,NA,NA
1167,NA,NA,NA
", file="Data2.csv")
Data = read.csv("Data2.csv",header=TRUE)
attach(Data)

bug.dat=list("Area", "Adults", "Egg", "Juveniles")

init.fun=function(){
   recruit_asymp = rnorm(1, 9000, 0.01)
   mua = rnorm(3,9000,0.0001)
   alpha = rnorm(1,0.0015,0.0002)
   Mspawn=rnorm(1,2.391,0.01)
   tau = rgamma(1, 2.82e9, 0.00000595)
   mu = rnorm(3, 10000, 1)
   #mu = rpois(1,1000000)
   #muj = rpois(1,1000)
   muj = rnorm(3, 1000, 10)
   tauj = rgamma(1, 1000, 0.01)
   Egg = rnorm(1, 1000000, 1000)
   Juveniles = rnorm(1, 1000, 100)
   Adults = rnorm(1,9000,800)
   return(list(Egg = Egg, Juveniles = Juveniles, mu = mu, 
  alpha = alpha, muj = muj, Adults = Adults
  ))}
        

params=c("Adults","mu", "muj", "recruit_asymp", "Juveniles", "Egg")
Sys.time()

EggBug=bugs(data = bug.dat, inits = init.fun, parameters.to.save = params, model.file="BarberModelwithArea.txt",
    n.chains=3, n.iter=150000, n.burnin=100000,
    n.thin=10, codaPkg = TRUE
    )
Sys.time()

codaegg <-read.bugs(EggBug)
plot(codaegg)
gelman.plot(codaegg)
gelman.diag(codaegg)

summary(codaegg)

#running Blueback herring in scenario 1. Making sure the herring does not exceed 20% of Alewife carrying capacity
cat("
model { 
  for (i in 1:3) {
  Egg[i] ~ dnorm(mu[i], tau)
  mu[i] <- ((Adults[i]*exp(-0.25*Mspawn)*0.5)*fecundity)
    Juveniles[i] ~ dnorm(muj[i], tauj) 
    muj[i] <- ((alpha*(mu[i]))/(1 + ((alpha*mu[i])/recruit_asymp)))
  }
fecundity <- (b*w) - c
b <- 871.72 
c <- 50.916
w <- 144
recruit_asymp ~ dnorm(9249, 0.01)
#recruit_asymp ~ dnorm(3283, 11)
Mspawn ~ dnorm(2.391,1)
#Mspawn <- 2.391
tau ~ dgamma(100000, 0.00000595)
tauj ~ dgamma (500, 0.001)
alpha ~ dnorm(0.0015, 1)
}
", file = "BarberModelEggsNoIE.txt")
modelCheck("BarberModelEggsNoIE.txt")
#### Data  ##########
cat(
"Adults,Egg,Juveniles
486500,16788683128,3282.572
53700,NA,NA
246750,NA,NA
", file="Data2.csv")
Data = read.csv("Data2.csv",header=TRUE)
attach(Data)

bug.dat=list("Adults", "Egg", "Juveniles")

init.fun=function(){
   recruit_asymp = rnorm(1, 9000, 0.01)
   alpha = rnorm(1,0.0015,0.0002)
   Mspawn=rnorm(1,2.391,0.01)
   tau = rgamma(1, 2.82e9, 0.00000595)
   mu = rnorm(1, 10000, 1)
   #mu = rpois(1,1000000)
   #muj = rpois(1,1000)
   muj = rnorm(1, 1000, 10)
   tauj = rgamma(1, 1000, 0.01)
   Egg = rnorm(1, 1000000, 1000)
   Juveniles = rnorm(1, 1000, 0.1)
   return(list(Egg = Egg, Juveniles = Juveniles, mu = mu, 
  alpha = alpha, muj = muj
  ))}
        

params=c("mu", "muj", "recruit_asymp")
Sys.time()

EggBug=bugs(data = bug.dat, inits = init.fun, parameters.to.save = params, model.file="BarberModelEggsNoIE.txt",
    n.chains=3, n.iter=100000, n.burnin=50000,
    n.thin=100, codaPkg = TRUE
    )
Sys.time()

codaegg <-read.bugs(EggBug)
plot(codaegg)
gelman.plot(codaegg)
gelman.diag(codaegg)


summary(codaegg)


#Nash Lake
#1600 acres total for BH and Alewife
#687 acres

#updating Scenario 1 
#adding stochastic carrying capacity 
#adding CIs from Devine paper to the asymptote with gamma
#precision for adult carrying capacity 
#(1186-672)/1.284 = 400 
#Variance = 400^2 = 1600
#precision = 1/1600 = 0.0006

cat("
model { 
  for (i in 1:3) {
  Adults[i] ~ dnorm(mua[i], 0.0006)
  mua[i] <- 846*Area[i]
  Egg[i] ~ dnorm(mu[i], tau)
  mu[i] <- ((Adults[i]*exp(-0.25*Mspawn)*0.5)*fecundity)
    Juveniles[i] ~ dnorm(muj[i], tauj) 
    muj[i] <- ((alpha*(mu[i]))/(1 + ((alpha*mu[i])/recruit_asymp)))

  }
fecundity <- (b*w) - c
b <- 871.72 
c <- 50.916
w <- 144
recruit_asymp ~ dnorm(9249, 0.01)
#recruit_asymp ~ dnorm(3283, 11)
Mspawn ~ dnorm(2.391,1)
#Mspawn <- 2.391
tau ~ dgamma(100000, 0.00000595)
tauj ~ dgamma (500, 0.001)
alpha ~ dnorm(0.0015, 1)
}
", file = "BarberModelwithArea.txt")
modelCheck("BarberModelwithArea.txt")
#### Data  ##########
cat(
"Area,Adults,Egg,Juveniles
1167,486500,16788683128,3282.572
687,NA,NA,NA
1600,NA,NA,NA
", file="Data2.csv")
Data = read.csv("Data2.csv",header=TRUE)
attach(Data)

bug.dat=list("Area", "Adults", "Egg", "Juveniles")

init.fun=function(){
   recruit_asymp = rnorm(1, 9000, 0.01)
   mua = rnorm(3,9000,0.0001)
   alpha = rnorm(1,0.0015,0.0002)
   Mspawn=rnorm(1,2.391,0.01)
   tau = rgamma(1, 2.82e9, 0.00000595)
   mu = rnorm(3, 10000, 1)
   #mu = rpois(1,1000000)
   #muj = rpois(1,1000)
   muj = rnorm(3, 1000, 10)
   tauj = rgamma(1, 1000, 0.01)
   Egg = rnorm(1, 1000000, 1000)
   Juveniles = rnorm(1, 1000, 100)
   Adults = rnorm(1,9000,800)
   return(list(Egg = Egg, Juveniles = Juveniles, mu = mu, 
  alpha = alpha, muj = muj, Adults = Adults
  ))}
        

params=c("Adults","mu", "muj", "recruit_asymp", "Juveniles", "Egg")
Sys.time()

EggBug=bugs(data = bug.dat, inits = init.fun, parameters.to.save = params, model.file="BarberModelwithArea.txt",
    n.chains=3, n.iter=150000, n.burnin=100000,
    n.thin=10, codaPkg = TRUE
    )
Sys.time()

codaegg <-read.bugs(EggBug)
plot(codaegg)
gelman.plot(codaegg)
gelman.diag(codaegg)

summary(codaegg)

#running Blueback herring in scenario 3. Making sure the herring does not exceed 20% of Alewife carrying capacity
cat("
model { 
  for (i in 1:3) {
  Egg[i] ~ dnorm(mu[i], tau)
  mu[i] <- ((Adults[i]*exp(-0.25*Mspawn)*0.5)*fecundity)
    Juveniles[i] ~ dnorm(muj[i], tauj) 
    muj[i] <- ((alpha*(mu[i]))/(1 + ((alpha*mu[i])/recruit_asymp)))
  }
fecundity <- (b*w) - c
b <- 871.72 
c <- 50.916
w <- 144
recruit_asymp ~ dnorm(9249, 0.01)
#recruit_asymp ~ dnorm(3283, 11)
Mspawn ~ dnorm(2.391,1)
#Mspawn <- 2.391
tau ~ dgamma(100000, 0.00000595)
tauj ~ dgamma (500, 0.001)
alpha ~ dnorm(0.0015, 1)
}
", file = "BarberModelEggsNoIE.txt")
modelCheck("BarberModelEggsNoIE.txt")
#### Data  ##########
cat(
"Adults,Egg,Juveniles
486500,16788683128,3282.572
145300,NA,NA
338500,NA,NA
", file="Data2.csv")
Data = read.csv("Data2.csv",header=TRUE)
attach(Data)

bug.dat=list("Adults", "Egg", "Juveniles")

init.fun=function(){
   recruit_asymp = rnorm(1, 9000, 0.01)
   alpha = rnorm(1,0.0015,0.0002)
   Mspawn=rnorm(1,2.391,0.01)
   tau = rgamma(1, 2.82e9, 0.00000595)
   mu = rnorm(1, 10000, 1)
   #mu = rpois(1,1000000)
   #muj = rpois(1,1000)
   muj = rnorm(1, 1000, 10)
   tauj = rgamma(1, 1000, 0.01)
   Egg = rnorm(1, 1000000, 1000)
   Juveniles = rnorm(1, 1000, 0.1)
   return(list(Egg = Egg, Juveniles = Juveniles, mu = mu, 
  alpha = alpha, muj = muj
  ))}
        

params=c("mu", "muj", "recruit_asymp")
Sys.time()

EggBug=bugs(data = bug.dat, inits = init.fun, parameters.to.save = params, model.file="BarberModelEggsNoIE.txt",
    n.chains=3, n.iter=100000, n.burnin=50000,
    n.thin=100, codaPkg = TRUE
    )
Sys.time()

codaegg <-read.bugs(EggBug)
plot(codaegg)
gelman.plot(codaegg)
gelman.diag(codaegg)


summary(codaegg)