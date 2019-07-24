# remove all variables
rm(list=ls())

set.seed(1)    # to produce reproducable results

# Parameter
lgd = 0.45    
r = 0.035     
m = 100000     # number of simulations
alpha = 0.99  # confidence level for CVaR
t = 1         # time horizon

# possible ratig classes
rc = 1:8

firmname = c("A", "B", "C") #could be anything

# loan portfolio
ead = c(1E6, 12E6, 6E6) 
rating = c(1, 4, 3)   
n = length(rating) # number of companies

# migration matrix
M = matrix(c(  0.9081,  0.0833,0.0068, 0.0006,0.0008, 0.0002,0.0001,0.0001,
               0.0070,  0.9065,0.0779, 0.0064,0.0006, 0.0013,0.0002,0.0001,
               0.0009,0.0227,0.9105, 0.0552,0.0074, 0.0026,0.0001,0.0006,
               0.0002,0.0033,0.0595,0.8593,0.0530, 0.0117,0.0112,0.0018,
               0.0003,0.0014,0.0067,0.0773, 0.8053, 0.0884,0.0100,0.0106,
               0.0001,0.0011,0.0024,0.0043,0.0648, 0.8346,0.0407,0.0520,
               0.0020,0.0001,0.0022,0.0130,0.0238, 0.1124,    0.6486,0.1979,
               0,     0,0, 0,  0, 0,  0,1.0000),
            8, 8, dimnames = list(rc, rc), byrow=T)

# correlation matrix msut have dim of number of companies
rho = matrix(c(1,   0.6, 0.3,
               0.6,1,   0.5,  
               0.3,0.5, 1),
            3, 3, dimnames=list(firmname, firmname), byrow=T)



# *** Step 1: Bounds ***

# reverse odering of migration matrix, default class is omitted
revM = M[1:7, 8:1]

# building cumulated sums and transposing matrix 
cumRevM = t(apply(revM, 1, cumsum))

# computation of bounds using the inverse of the CDF of the standard normal distribution
S = qnorm(cumRevM)

# Upper bounds set to Inf explicitly
S[, 8] = +Inf

# Setting lower bound
S = cbind(-Inf, S)


# *** Step 2: Computation of CS ***

# PDs are resorted
pd = M[1:7, 8]
# CS Formula
cs = -(log(1 - lgd*pd))


# *** STEP 3: Computation of Portfolio Values ***

# value in case of no change in rating
p = ead*exp(-(r + cs[rating])*t)
# value in case of no change in rating for portfolio of loans
p0 = sum(p)


# *** Step 4: Computation of State Space ***


cs2 = matrix(rep(cs, n), n, 7, byrow=T)
ead2 = matrix(rep(ead, 7), n, 7, byrow=F)
V = ead2*exp(-(r + cs2)*t)
# special case default
V = cbind(V, ead*(1 - lgd))


# *** Step 5: Simultion of random numbers ***

# Simulationen using Antithetic Sampling
X = matrix(rnorm(n * m / 2), n, m/2)
Z = cbind(X, -X)


# *** Step 6: Simulation using correlated random numbers ***

# Cholesky Matrix
A = t(chol(rho))
Y = A%*%Z


# *** Step 7: Value for every scenario ***

# allocation of secanrios to rating classes
simV = NULL
for (i in 1:n) {
  # bound depend on the company
  l = S[rating[i], ]
  
  simClass = 9 - as.numeric(cut(Y[i, ], l))
  #allocation of vlaues
  simV = rbind(simV, V[i, simClass])
}


# *** Step 8: Compuation of simulated portfolio values ***

# summing over loans
sim.p = colSums(simV)


# *** Step 9: Constructing P/L distribution ***


sim.pv = sim.p - p0


# *** Step 10: CVaR***


CVaR = -quantile(sim.pv, 1 - alpha)


# *** Step 11: some illustrations ***

b = seq(min(sim.pv), max(sim.pv), length=40)
hist(sim.pv, main="P/L Histogramm", xlab="Profit/Loss", ylab="Counts", breaks=b)

CVaR


