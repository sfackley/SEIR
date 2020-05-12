require(deSolve) # loads the differential equation solving package
require(numDeriv)
require(ggplot2)
require(tuneR)
require(ggthemes)
require(tidyr)
require(dplyr)
require(magrittr)
rm(list=ls()) # clears the workspace


# model() is a function which defines the set of differential equations
# inputs: time, states, parameters
## time: the time(s) at which you want the derivates (a vector of times)
## states: the number of flies in each stage (a vector with the entries labeled)
## parameters: the parameters of the model (a vector with the entries labeled)
# outputs a list where the first entry is a vector derivatives and the second entry is NULL
# (this output format is required by the differential equation solver)
first.model <- function(time,states,parameters){
  with(as.list(c(states,parameters)),
       {
         beta.t <- R0*gamma
         if(time>lockdown.duration+start.time){
           beta.t <- R0*gamma*R0.reduction.off.lockdown
         }
         period <- intermittent.lockdown.duration+
           time.between.intermittent.lockdown
         first.lockdown <- time>=start.time&
           time<lockdown.duration+start.time
         int.lock.time <- (time-intermittent.lockdown.start.time) %% period
         intermittent.lockdown <- (int.lock.time<intermittent.lockdown.duration)&
           time>intermittent.lockdown.start.time
         if(first.lockdown){
           beta.t <- R0*gamma*R0.reduction
         }
         if(intermittent.lockdown){
           beta.t <- R0*gamma*R0.reduction.severe.lockdown
         }
         N  <- S + E + I + R
         dS <- -beta.t*S*I/N
         dE <- beta.t*S*I/N - alpha*E
         dI <- alpha*E - gamma*I 
         dR <- gamma*I
         list(c(dS,dE,dI,dR),NULL)
       })}

# gen.init() is a function which generates initial conditions of the model. 
gen.init <- function(parameters){
  with(as.list(parameters),
       {
         initial.conditions <- c(S=99.8,E=0.1,I=.1,R=0)	
         initial.conditions 
       })}

# do.sim() is a function which numerically solves the differential equation
# inputs: parameters of the model and times where you want to know the solution
# lsoda automatically determines the solving method and the time step to get the solution to within some tolerance at the specified times
# outputs a matrix with each row corresponding to the numbers of flies in each stages the input times
do.sim <- function(parameters,times){
  lsoda(gen.init(parameters),times,first.model,parameters)
}

# Lockdowns reduce R0 <- 
# Case finding increases gamma (decreases infectious time)
# Even off lockdown ppl are more cautious 

years <- 1
gamma <- 1/10
pars  <- c(R0=2.5,alpha=1/6,gamma=gamma,
           start.time=30,R0.reduction=1/3,lockdown.duration=60,
           R0.reduction.off.lockdown=2/3,
           R0.reduction.severe.lockdown=1/6,
           intermittent.lockdown.start.time=100,
           intermittent.lockdown.duration=30,
           time.between.intermittent.lockdown=60)
times <- 0:365*years

out.data <- do.sim(pars,times)
out.data %<>% as.data.frame
to.plot <- gather(out.data,
                  "State","Percent",-time)
to.plot %<>% filter(State=="I"|State=="R")

with(as.list(pars),{
lockdown.times <<- c(start.time,intermittent.lockdown.start.time,
                    intermittent.lockdown.start.time+
                      intermittent.lockdown.duration*1:1000+
                      time.between.intermittent.lockdown*1:1000)
lockdown.up <<- c(start.time+lockdown.duration,
                  intermittent.lockdown.start.time+
                    intermittent.lockdown.duration,
                     intermittent.lockdown.start.time+
                       intermittent.lockdown.duration*1:1000+
                       time.between.intermittent.lockdown*1:1000+
                        intermittent.lockdown.duration)
})

lockdown.times <- lockdown.times[lockdown.times<max(times)]
lockdown.up <- lockdown.up[lockdown.up<max(times)]

annotation=paste(names(pars),": ", round(pars,2),"\n",collapse=" ")
ggplot(to.plot,aes(time,Percent,group=State,color=State))+
  geom_line() + 
  geom_vline(xintercept=lockdown.times,color='firebrick',alpha=0.5)+
  geom_vline(xintercept=lockdown.up,color='forestgreen',alpha=0.5)+
  annotate("text",x=60,y=max(to.plot$Percent)*2/3, 
           label=annotation,size=2,family="Times")+
  theme_tufte()+xlab("Number of Days")



