## Load deSolve package
library(deSolve)
## Create an SIR function
sir <- function(time, state, parameters) {
  
  with(as.list(c(state, parameters)), {
    
    dS <- -betaI * S * I - betaW * S * W
    dI <-  betaI * S * I + betaW * S * W- gamma * I
    dW <-  xi * I - xi * W
    dR <-                 gamma * I
    
    return(list(c(dS, dI, dW, dR)))
  })
}
### Set parameters
## Proportion in each compartment: Susceptible 0.998856, Infected 0.00114, Recovered 0, Bacteria conc. in water 0
init       <- c(S = 1-1.14e-3, I = 1.14e-3, W = 0.0, R = 0.0)
## beta: infection parameter; gamma: recovery parameter
parameters <- c(betaI = 0.75,betaW = 0.75,xi = 0.01, gamma = 0.25)
## Time frame 
times      <- seq(0, 161, by = 7)
## Solve using ode (General Solver for Ordinary Differential Equations)
out <- ode(y = init, times = times, func = sir, parms = parameters)
## change to data frame
out <- as.data.frame(out)
## Delete time variable
out$time <- NULL
## Show data
head(out, 10)
## Plot
matplot(x = times, y = out, type = "l",
        xlab = "Time", ylab = "Susceptible and Recovered", main = "SIWR Model",
        lwd = 1, lty = 1, bty = "l", col = 2:5)
## Add legend
legend(40, 0.7, c("Susceptible", "Infected","Bacteria conc. in water", "Recovered"), pch = 1, col = 2:5, bty = "n")
#Plotting data and model in single graph
data <- c(102,80,75,200,400,2184,3842,5408,6224,6041,5379,4164,3226,2613,2096,1745,1390,1080,955,805,657,546,487,397)
y <- out$I*89193
y1=cbind(data,y)
as.data.frame(y1)
matplot(x = times, y = y1, type = "l", xlab = "Time", ylab = "Data and y", main = "Data vs y", lwd = 1, lty = 1,col=2:3)
########################################
#Multiple for loops
index<-matrix(nrow = 4,ncol = 0)
e<-matrix(nrow = 1,ncol = 0)
z<-matrix(nrow = 24,ncol=0)
for (l in seq(from = 95000, to=110000, by = 1000))
{
for(i in seq(from = 0.1, to = 1, by = 0.1))
{
  for(j in seq(from = 0.1, to = 1, by = 0.05))
  {
    for(k in seq(from = 0.01, to = 0.1, by = 0.01))
    {
      parameters <- c(betaI = i,betaW = j,xi = k, gamma = 0.25) 
      out <- ode(y = init, times = times, func = sir, parms = parameters)
      z<- cbind(z, out[,3]*l)
      e<-cbind(e,sum((data-out[,3]*l)^2))
      index <- cbind(index,matrix(c(i,j,k,l),ncol = 1))
    }
  } 
}
}
which.min(e)
e[which.min(e)]
index[,which.min(e)]
###############################Trial with estimated parameters
## beta: infection parameter; gamma: recovery parameter
parameters <- c(betaI = index[,which.min(e)][1],betaW = index[,which.min(e)][2],xi = index[,which.min(e)][3], gamma = 0.25)
## Solve using ode (General Solver for Ordinary Differential Equations)
out_new <- ode(y = init, times = times, func = sir, parms = parameters)
## change to data frame
out_new <- as.data.frame(out_new)
## Delete time variable
out_new$time <- NULL
## Show data
head(out_new, 10)
## Plot
matplot(x = times, y = out_new, type = "l",
        xlab = "Time", ylab = "Susceptible and Recovered", main = "SIWR Model",
        lwd = 1, lty = 1, bty = "l", col = 2:5)
## Add legend
legend(40, 0.7, c("Susceptible", "Infected","Bacteria conc. in water", "Recovered"), pch = 1, col = 2:5, bty = "n")
#Plotting data and model in single graph
data <- c(102,80,75,200,400,2184,3842,5408,6224,6041,5379,4164,3226,2613,2096,1745,1390,1080,955,805,657,546,487,397)
y <- out_new$I*index[,which.min(e)][4]
y1=cbind(data,y)
as.data.frame(y1) 
matplot(x = times, y = y, type = "l", xlab = "Time", ylab = "Data and y", main = "Data vs y", lwd = 1, lty = 1,col=2:3)
plot(times,y,type = 'l',col=2)
points(times,data,col=4)

