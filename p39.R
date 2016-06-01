# Complete Randomization for 2 Way ANOVA
# INPUT
#   factorNameA : Name of factor A
#   factorNameB : Name for factor B
#   nA : Number of levels for factor A
#   nB : Number of levels for factor B
#   nR : Number of replication for each combination
# RETURN : vector of randomzied combination of levels

CR2W = function(factorNameA="A", factorNameB="B", nA=2, nB=2, nR=1)
{
  A = paste0(factorNameA,1:nA)
  B = paste0(factorNameB,1:nB)
  sample(rep(outer(A, B, "paste0"),nR))
}

CR2W("A", "B", nA=3, nB=4, nR=2)

factorNameA = "A"
factorNameB = "B"
nA = 3
nB = 4
nR = 2
A = paste0(factorNameA,1:nA)
B = paste0(factorNameB,1:nB)
outer(A, B, "paste0")
as.vector(outer(A, B, "paste0"))
rep(outer(A, B, "paste0"),nR)


# Kendall's tau
Ktau = function(X, Y)
{
  m = length(X)
  n = length(Y)
  if (m != n | m < 1) { 
    print("Bad Input")
    return(0)
  }

  Tau = 0
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      Tau = Tau + sign(X[i] - X[j])*sign(Y[i] - Y[j])
    }
  }
  return(Tau)
}

X = c(86, 71, 77, 68, 91, 72, 77, 91, 70, 71, 88, 87)
Y = c(88, 77, 76, 64, 96, 72, 65, 90, 65, 80, 81, 72)

Ktau(X, Y)


X2 = sort(X)
Y2 = Y[order(X,Y)]
X2
Y2
n = length(X2)
P = vector(length=n)
Q = vector(length=n)

for (i in 1:(n-1)) {
  tP = 0
  tQ = 0
  for (j in (i+1):n) {
    tX = (X2[i] - X2[j]) * (Y2[i] - Y2[j])
    if (tX > 0) tP = tP + 1
    if (tX < 0) tQ = tQ + 1 
  }
  P[i] = tP
  Q[i] = tQ
}

cbind(P, Q)


