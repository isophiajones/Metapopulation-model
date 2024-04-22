###########
# m=0
###########

set.seed(123)
tspan = 100
r1 = 2.1 # intrinsic growth rate of source population 1
r2 = 2.0 # intrinsic growth rate of source population 2
r3 = 2.2 # intrinsic growth rate of source population 3
r4 = -2.1 # intrinsic growth rate of sink population
K = 20 # carrying capacity
m = 0 # dispersal rate
extinction = 0.01 # extinction threshold

# memory preallocation
x1 = matrix(0, 1, tspan) 
x2 = matrix(0, 1, tspan) 
x3 = matrix(0, 1, tspan) 
x4 = matrix(0, 1, tspan)
# initial population sizes
x1[1] = runif(1);
x2[1] = runif(1);
x3[1] = runif(1);
x4[1] = runif(1);

for(t in 1:(tspan-1)){
  # density dependent growth
  x1[t+1] = x1[t]*exp(r1*(1-(x1[t]/K)))
  x2[t+1] = x2[t]*exp(r2*(1-(x2[t]/K)))
  x3[t+1] = x3[t]*exp(r3*(1-(x3[t]/K)))
  x4[t+1] = x4[t]*exp(r4*(1-(x4[t]/K)))
  # dispersal
  x1disp = m*x1[t+1] # inds. leaving patch 1
  x2disp = m*x2[t+1] # inds. leaving patch 2
  x3disp = m*x3[t+1] # inds. leaving patch 3
  x4disp = m*x4[t+1] # inds. leaving patch 4
  # new population size after growth and dispersal
  x1[t+1] = x1[t+1] + (x2disp/3) + (x3disp/3) + (x4disp/3) - x1disp
  x2[t+1] = x2[t+1] + (x1disp/3) + (x3disp/3) + (x4disp/3) - x2disp
  x3[t+1] = x3[t+1] + (x1disp/3) + (x2disp/3) + (x4disp/3) - x3disp
  x4[t+1] = x4[t+1] + (x1disp/3) + (x2disp/3) + (x3disp/3) - x4disp
  # check for extinction events
  x1[x1<extinction] = 0
  x2[x2<extinction] = 0
  x3[x3<extinction] = 0
  x4[x4<extinction] = 0
}

par(mfrow = c(2,2))
plot(1:tspan, x1, 
     type = "l",
     col = "darkmagenta",
     xlab = "Time",
     ylab = "Abundance",
     main = "Patch 1")
plot(1:tspan, x2,
     type = "l",
     col = "darkcyan",
     xlab = "Time",
     ylab = "Abundance",
     main = "Patch 2")
plot(1:tspan, x3,
     type = "l",
     col = "midnightblue",
     xlab = "Time",
     ylab = "Abundance",
     main = "Patch 3")
plot(1:tspan, x4,
     type = "l",
     col = "mediumslateblue",
     xlab = "Time",
     ylab = "Abundance",
     main = "Patch 4")

###########
# m=0.05
###########

set.seed(123)
m = 0.05 # dispersal rate

# memory preallocation
x1 = matrix(0, 1, tspan) 
x2 = matrix(0, 1, tspan) 
x3 = matrix(0, 1, tspan) 
x4 = matrix(0, 1, tspan)
# initial population sizes
x1[1] = runif(1);
x2[1] = runif(1);
x3[1] = runif(1);
x4[1] = runif(1);

for(t in 1:(tspan-1)){
  # density dependent growth
  x1[t+1] = x1[t]*exp(r1*(1-(x1[t]/K)))
  x2[t+1] = x2[t]*exp(r2*(1-(x2[t]/K)))
  x3[t+1] = x3[t]*exp(r3*(1-(x3[t]/K)))
  x4[t+1] = x4[t]*exp(r4*(1-(x4[t]/K)))
  # dispersal
  x1disp = m*x1[t+1] # inds. leaving patch 1
  x2disp = m*x2[t+1] # inds. leaving patch 2
  x3disp = m*x3[t+1] # inds. leaving patch 3
  x4disp = m*x4[t+1] # inds. leaving patch 4
  # new population size after growth and dispersal
  x1[t+1] = x1[t+1] + (x2disp/3) + (x3disp/3) + (x4disp/3) - x1disp
  x2[t+1] = x2[t+1] + (x1disp/3) + (x3disp/3) + (x4disp/3) - x2disp
  x3[t+1] = x3[t+1] + (x1disp/3) + (x2disp/3) + (x4disp/3) - x3disp
  x4[t+1] = x4[t+1] + (x1disp/3) + (x2disp/3) + (x3disp/3) - x4disp
  # check for extinction events
  x1[x1<extinction] = 0
  x2[x2<extinction] = 0
  x3[x3<extinction] = 0
  x4[x4<extinction] = 0
}

par(mfrow = c(2,2))
plot(1:tspan, x1, 
     type = "l",
     col = "darkmagenta",
     xlab = "Time",
     ylab = "Abundance",
     main = "Patch 1")
plot(1:tspan, x2,
     type = "l",
     col = "darkcyan",
     xlab = "Time",
     ylab = "Abundance",
     main = "Patch 2")
plot(1:tspan, x3,
     type = "l",
     col = "midnightblue",
     xlab = "Time",
     ylab = "Abundance",
     main = "Patch 3")
plot(1:tspan, x4,
     type = "l",
     col = "mediumslateblue",
     xlab = "Time",
     ylab = "Abundance",
     main = "Patch 4")

###########
# m=0.10
###########

set.seed(123)
m = 0.10 # dispersal rate

# memory preallocation
x1 = matrix(0, 1, tspan) 
x2 = matrix(0, 1, tspan) 
x3 = matrix(0, 1, tspan) 
x4 = matrix(0, 1, tspan)
# initial population sizes
x1[1] = runif(1);
x2[1] = runif(1);
x3[1] = runif(1);
x4[1] = runif(1);

for(t in 1:(tspan-1)){
  # density dependent growth
  x1[t+1] = x1[t]*exp(r1*(1-(x1[t]/K)))
  x2[t+1] = x2[t]*exp(r2*(1-(x2[t]/K)))
  x3[t+1] = x3[t]*exp(r3*(1-(x3[t]/K)))
  x4[t+1] = x4[t]*exp(r4*(1-(x4[t]/K)))
  # dispersal
  x1disp = m*x1[t+1] # inds. leaving patch 1
  x2disp = m*x2[t+1] # inds. leaving patch 2
  x3disp = m*x3[t+1] # inds. leaving patch 3
  x4disp = m*x4[t+1] # inds. leaving patch 4
  # new population size after growth and dispersal
  x1[t+1] = x1[t+1] + (x2disp/3) + (x3disp/3) + (x4disp/3) - x1disp
  x2[t+1] = x2[t+1] + (x1disp/3) + (x3disp/3) + (x4disp/3) - x2disp
  x3[t+1] = x3[t+1] + (x1disp/3) + (x2disp/3) + (x4disp/3) - x3disp
  x4[t+1] = x4[t+1] + (x1disp/3) + (x2disp/3) + (x3disp/3) - x4disp
  # check for extinction events
  x1[x1<extinction] = 0
  x2[x2<extinction] = 0
  x3[x3<extinction] = 0
  x4[x4<extinction] = 0
}

par(mfrow = c(2,2))
plot(1:tspan, x1, 
     type = "l",
     col = "darkmagenta",
     xlab = "Time",
     ylab = "Abundance",
     main = "Patch 1")
plot(1:tspan, x2,
     type = "l",
     col = "darkcyan",
     xlab = "Time",
     ylab = "Abundance",
     main = "Patch 2")
plot(1:tspan, x3,
     type = "l",
     col = "midnightblue",
     xlab = "Time",
     ylab = "Abundance",
     main = "Patch 3")
plot(1:tspan, x4,
     type = "l",
     col = "mediumslateblue",
     xlab = "Time",
     ylab = "Abundance",
     main = "Patch 4")

###########
# m=0.2
###########

set.seed(123)
m = 0.2 # dispersal rate

# memory preallocation
x1 = matrix(0, 1, tspan) 
x2 = matrix(0, 1, tspan) 
x3 = matrix(0, 1, tspan) 
x4 = matrix(0, 1, tspan)
# initial population sizes
x1[1] = runif(1);
x2[1] = runif(1);
x3[1] = runif(1);
x4[1] = runif(1);

for(t in 1:(tspan-1)){
  # density dependent growth
  x1[t+1] = x1[t]*exp(r1*(1-(x1[t]/K)))
  x2[t+1] = x2[t]*exp(r2*(1-(x2[t]/K)))
  x3[t+1] = x3[t]*exp(r3*(1-(x3[t]/K)))
  x4[t+1] = x4[t]*exp(r4*(1-(x4[t]/K)))
  # dispersal
  x1disp = m*x1[t+1] # inds. leaving patch 1
  x2disp = m*x2[t+1] # inds. leaving patch 2
  x3disp = m*x3[t+1] # inds. leaving patch 3
  x4disp = m*x4[t+1] # inds. leaving patch 4
  # new population size after growth and dispersal
  x1[t+1] = x1[t+1] + (x2disp/3) + (x3disp/3) + (x4disp/3) - x1disp
  x2[t+1] = x2[t+1] + (x1disp/3) + (x3disp/3) + (x4disp/3) - x2disp
  x3[t+1] = x3[t+1] + (x1disp/3) + (x2disp/3) + (x4disp/3) - x3disp
  x4[t+1] = x4[t+1] + (x1disp/3) + (x2disp/3) + (x3disp/3) - x4disp
  # check for extinction events
  x1[x1<extinction] = 0
  x2[x2<extinction] = 0
  x3[x3<extinction] = 0
  x4[x4<extinction] = 0
}

par(mfrow = c(2,2))
plot(1:tspan, x1, 
     type = "l",
     col = "darkmagenta",
     xlab = "Time",
     ylab = "Abundance",
     main = "Patch 1")
plot(1:tspan, x2,
     type = "l",
     col = "darkcyan",
     xlab = "Time",
     ylab = "Abundance",
     main = "Patch 2")
plot(1:tspan, x3,
     type = "l",
     col = "midnightblue",
     xlab = "Time",
     ylab = "Abundance",
     main = "Patch 3")
plot(1:tspan, x4,
     type = "l",
     col = "mediumslateblue",
     xlab = "Time",
     ylab = "Abundance",
     main = "Patch 4")

###########
# m=0.35
###########

set.seed(123)
m = 0.35 # dispersal rate

# memory preallocation
x1 = matrix(0, 1, tspan) 
x2 = matrix(0, 1, tspan) 
x3 = matrix(0, 1, tspan) 
x4 = matrix(0, 1, tspan)
# initial population sizes
x1[1] = runif(1);
x2[1] = runif(1);
x3[1] = runif(1);
x4[1] = runif(1);

for(t in 1:(tspan-1)){
  # density dependent growth
  x1[t+1] = x1[t]*exp(r1*(1-(x1[t]/K)))
  x2[t+1] = x2[t]*exp(r2*(1-(x2[t]/K)))
  x3[t+1] = x3[t]*exp(r3*(1-(x3[t]/K)))
  x4[t+1] = x4[t]*exp(r4*(1-(x4[t]/K)))
  # dispersal
  x1disp = m*x1[t+1] # inds. leaving patch 1
  x2disp = m*x2[t+1] # inds. leaving patch 2
  x3disp = m*x3[t+1] # inds. leaving patch 3
  x4disp = m*x4[t+1] # inds. leaving patch 4
  # new population size after growth and dispersal
  x1[t+1] = x1[t+1] + (x2disp/3) + (x3disp/3) + (x4disp/3) - x1disp
  x2[t+1] = x2[t+1] + (x1disp/3) + (x3disp/3) + (x4disp/3) - x2disp
  x3[t+1] = x3[t+1] + (x1disp/3) + (x2disp/3) + (x4disp/3) - x3disp
  x4[t+1] = x4[t+1] + (x1disp/3) + (x2disp/3) + (x3disp/3) - x4disp
  # check for extinction events
  x1[x1<extinction] = 0
  x2[x2<extinction] = 0
  x3[x3<extinction] = 0
  x4[x4<extinction] = 0
}

par(mfrow = c(2,2))
plot(1:tspan, x1, 
     type = "l",
     col = "darkmagenta",
     xlab = "Time",
     ylab = "Abundance",
     main = "Patch 1")
plot(1:tspan, x2,
     type = "l",
     col = "darkcyan",
     xlab = "Time",
     ylab = "Abundance",
     main = "Patch 2")
plot(1:tspan, x3,
     type = "l",
     col = "midnightblue",
     xlab = "Time",
     ylab = "Abundance",
     main = "Patch 3")
plot(1:tspan, x4,
     type = "l",
     col = "mediumslateblue",
     xlab = "Time",
     ylab = "Abundance",
     main = "Patch 4")