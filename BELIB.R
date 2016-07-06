# R Libraries for Equivalence Test

# Diletti E, Hauscheke D, Steinjans VW. Sample size determination for bioequivalence assessment by means of confidence intervals. Int J Clin Pharmacology, Therapy and Toxicology. 1992;29:S51-58
# Chow SC, Liu JP. Design and Analysis of Bioavailability and Bioequivalence Study. 2nd ed. p158. 2000, Marcel Dekker Inc

##################################
# Sample Size for BE : Difference

bpow.d.cv <- function(i, cv, delta=0, alpha=0.1, beta=0.2, thetaL=-20, thetaU=20)
{
  if (delta <= thetaL | delta >= thetaU) return(0);
  t0 <- qt(1-alpha/2, 2*i-2);
  p1 <- pt(-1*t0, 2*i-2, ncp=(thetaL-delta)/(cv/sqrt(i)))
  p2 <- pt(t0, 2*i-2, ncp=(thetaU-delta)/(cv/sqrt(i)))
  power <- p1 - p2;
  if (power < 0) power <- 0;  
  return(power);
}


bpow.d.mse <- function(i, mu.r, mse, true.d=0, alpha=0.1, beta=0.2, thetaL=-20, thetaU=20)
{
  cv <- 100*sqrt(mse) / mu.r;
  delta <- 100*true.d / mu.r;
  return(bpow.d.cv(i, cv, delta, alpha, beta, thetaL, thetaU));
}


bss.d.cv <- function(cv, delta=0, alpha=0.1, beta=0.2, thetaL=-20, thetaU=20)
{
  if (delta <= thetaL | delta >= thetaU) return(Inf);
  for(i in 2:1000) {
    power <- bpow.d.cv(i, cv, delta, alpha, beta, thetaL, thetaU);
    if (power > 1 - beta) return(i);
  }
  return(">1000");
}


bss.d.mse <- function(mu.r, mse, true.d=0, alpha=0.1, beta=0.2, thetaL=-20, thetaU=20)
{
  cv <- 100*sqrt(mse) / mu.r;
  delta <- 100*true.d / mu.r;
  return(bss.d.cv(cv, delta, alpha, beta, thetaL, thetaU));
}

# LL, UL like 85% - 115%
# N: total N, 2 * n per group
bss.d.ci <- function(N, LL, UL)
{
  pe <- (LL + UL)/2
  t0 <- qt(0.95, N-2);
  cv <- sqrt(N/2)*(UL - LL)/(2*t0);
  s1 <- bss.d.cv(cv)
  s2 <- bss.d.cv(cv, delta=(pe-100));
  sampsize <- cbind(s1, s2);
  p1 <- round(100 * bpow.d.cv(N/2, cv));
  p2 <- round(100 * bpow.d.cv(N/2, cv, delta=(pe-100)));
  power <- cbind(p1, p2);
  result <- rbind(sampsize, power);
  dimnames(result) <- list(c("80% Power Sample Size", paste("Power at N =",N)), c("True Percent=100", sprintf("True Percent=%.2f", pe)));
  return(result);
}

############################
# Sample Size for BE : Ratio

bpow.r.mse <- function(i, mse, true.r=1, alpha=0.1, beta=0.2, thetaL=0.8, thetaU=1.25)
{
  if (true.r <= thetaL | true.r >= thetaU) return(0);
  t0 <- qt(1-alpha/2,2*i-2);
  p1 <- pt(-1*t0, 2*i-2, ncp=log(thetaL/true.r)/sqrt(mse/i));
  p2 <- pt(t0, 2*i-2, ncp=log(thetaU/true.r)/sqrt(mse/i));
  power <- p1 - p2;
  if (power < 0) power <- 0;
  return(power);
}


bpow.r.cv <- function(i, cv, true.r=1, alpha=0.1, beta=0.2, thetaL=0.8, thetaU=1.25)
{
  mse <- log(1+(cv/100)^2);
  return(bpow.r.mse(i, mse, true.r, alpha, beta, thetaL, thetaU));
}


bss.r.mse <- function(mse, true.r=1, alpha=0.1, beta=0.2, thetaL=0.8, thetaU=1.25)
{
  if (true.r <= thetaL | true.r >= thetaU) return(Inf);
  for (i in 2:1000) {
    power <- bpow.r.mse(i, mse, true.r, alpha, beta, thetaL, thetaU);
    if (power > 1 - beta) return(i);
  }
  return(">1000");
}


bss.r.cv <- function(cv, true.r=1, alpha=0.1, beta=0.2, thetaL=0.8, thetaU=1.25)
{
  mse <- log(1+(cv/100)^2);
  return(bss.r.mse(mse, true.r, alpha, beta, thetaL, thetaU));
}

# LL, UL like 0.85 ~ 1.15
# N: toal N, 2 * n per group
bss.r.ci <- function(N, LL, UL)
{
  pe <- exp((log(UL)+log(LL))/2);
  t0 <- qt(0.95, N-2);
  sd <- (log(UL)-log(LL))/(2*t0);
  mse <- sd^2 * N / 2;
  s1 <- bss.r.mse(mse)
  s2 <- bss.r.mse(mse, true.r=pe);
  sampsize <- cbind(s1, s2);
  p1 <- round(100 * bpow.r.mse(N/2, mse));
  p2 <- round(100 * bpow.r.mse(N/2, mse, true.r=pe));
  power <- cbind(p1, p2);
  result <- rbind(sampsize, power);
  dimnames(result) <- list(c("80% Power Sample Size", paste("Power at N =",N)), c("True Ratio=1", sprintf("True Ratio=%.4f", pe)));
  return(result);
}

#########################
# 2x2 BE Test

assert <- function(bedata)
{
  Si11 <- bedata[bedata$GRP=="RT" & bedata$PRD==1, "SUBJ"]
  Si21 <- bedata[bedata$GRP=="RT" & bedata$PRD==2, "SUBJ"]
  Si12 <- bedata[bedata$GRP=="TR" & bedata$PRD==1, "SUBJ"]
  Si22 <- bedata[bedata$GRP=="TR" & bedata$PRD==2, "SUBJ"]

  return(identical(Si11, Si21) & identical(Si12, Si22))
}

betest <- function(bedata, var, logtransformed)
{
  if(!assert(bedata)) {
    cat("\n Bad Data Format !\n");
    return(NULL);
  }

  Yijk <- bedata[, var]
  Yi11 <- bedata[bedata$GRP=="RT" & bedata$PRD==1, var]
  Yi21 <- bedata[bedata$GRP=="RT" & bedata$PRD==2, var]
  Yi12 <- bedata[bedata$GRP=="TR" & bedata$PRD==1, var]
  Yi22 <- bedata[bedata$GRP=="TR" & bedata$PRD==2, var]

  n1 <- length(Yi11)
  n2 <- length(Yi12)

  Y... <- mean(Yijk)
  SStotal <- sum((Yijk - Y...)^2)

  Y.11 <- mean(Yi11)
  Y.21 <- mean(Yi21)
  Y.12 <- mean(Yi12)
  Y.22 <- mean(Yi22)

  Yi.1 <- (Yi11 + Yi21) / 2
  Yi.2 <- (Yi12 + Yi22) / 2

  Y..1 <- mean(Yi.1);
  Y..2 <- mean(Yi.2);

  mu.r <- (Y.11 + Y.22) / 2
  mu.t <- (Y.21 + Y.12) / 2

  SScarry   <- 2*n1*n2/(n1+n2)*(Y.12 + Y.22 - Y.11 - Y.21)^2 / 4
#  SSinter   <- (sum(Yi.1^2) + sum(Yi.2^2) - Y..1^2 * n1  - Y..2^2 * n2) * 2
  SSinter   <- (sum((Yi.1-Y..1)^2) + sum((Yi.2-Y..2)^2)) * 2
  SSbetween <- SScarry + SSinter

  SSperiod  <- 2*n1*n2/(n1+n2)*(Y.21 + Y.22 - Y.11 - Y.12)^2 / 4
  SSdrug    <- 2*n1*n2/(n1+n2)*(Y.21 + Y.12 - Y.11 - Y.22)^2 / 4
  SSintra   <- SStotal - SScarry - SSinter - SSdrug - SSperiod

  Source <- c("SUBJECT", "GROUP", "SUBJECT(GROUP)", "PERIOD", "DRUG", "ERROR", "TOTAL");
  SS     <- c(SSbetween, SScarry, SSinter, SSperiod, SSdrug, SSintra, SStotal);
  DF     <- c(n1+n2-1, 1, n1+n2-2, 1, 1, n1+n2-2, 2*n1+2*n2-1);
  MS     <- SS / DF
  mse    <- SSintra / (n1+n2-2);
  F      <- MS / c(mse, MS[3], mse, mse, mse, mse, mse);
  p1 <- 1 - pf(F[1], n1+n2-1, n1+n2-2)
  p2 <- 1 - pf(F[2], 1, n1+n2-2);
  p3 <- 1 - pf(F[3], n1+n2-2, n1+n2-2);
  p4 <- 1 - pf(F[4], 1, n1+n2-2);
  p5 <- 1 - pf(F[5], 1, n1+n2-2);
  p  <- c(p1, p2, p3, p4, p5, NA, NA)
  F[6] <- F[7] <- MS[7] <- NA

  ANOVA <- cbind(SS, DF, MS, F, p)
  dimnames(ANOVA) <- list(Source,c("SS", "DF", "MS", "F", "p"))

  pe <-  mu.t - mu.r
  sd <- sqrt(mse / 2 * (1/n1 + 1/n2))   # See pp 62-63
  t0 <- qt(0.95, n1+n2-2);
  ci0 <- cbind(pe - t0 * sd, pe, pe + t0 * sd)

  if (logtransformed == TRUE) {
    lsm <- cbind(exp(mu.r), exp(mu.t))
    dimnames(lsm) <- list("Geometric Means", cbind("Reference Drug", "Test Drug"))

    ci <- exp(ci0);
    dimnames(ci) <- list("90% CI for Ratio", c("Lower Limit", "Point Estimate", "Upper Limit"));

    sampsize1 <- bss.r.mse(mse);
    sampsize2 <- bss.r.mse(mse, true.r=exp(pe));
    ss <- cbind(sampsize1, sampsize2)
    dimnames(ss) <- list("80% Power Sample Size", c("True Ratio=1", "True Ratio=Point Estimate"));
  } else {
    lsm <- cbind(mu.r, mu.t)
    dimnames(lsm) <- list("Arithmetic Means", cbind("Reference Drug", "Test Drug"))

    ci1 <- (1 + ci0 / mu.r) * 100
    ci <- rbind(ci0, ci1)
    dimnames(ci) <- list(c("90% CI for Difference", "90% CI for Difference(%)"), c("Lower Limit", "Point Estimate", "Upper Limit"));

    sampsize1 <- bss.d.mse(mu.r, mse);
    sampsize2 <- bss.d.mse(mu.r, mse, true.d=pe);
    ss <- cbind(sampsize1, sampsize2)
    dimnames(ss) <- list("80% Power Sample Size", c("True Difference=0", "True Difference=Point Estimate"));
  }

  result <- list(ANOVA, lsm, ci, ss);
  names(result) <- c("Analysis of Variance", "Least Square Means", "90% Confidence Interval", "Sample Size")

  return(result);
}

hodges <- function(bedata, var)
{
  if(!assert(bedata)) {
    cat("\n Bad Data Format !\n");
    return(NULL);
  }

  Yi11 <- bedata[bedata$GRP=="RT" & bedata$PRD==1, var]
  Yi21 <- bedata[bedata$GRP=="RT" & bedata$PRD==2, var]
  Yi12 <- bedata[bedata$GRP=="TR" & bedata$PRD==1, var]
  Yi22 <- bedata[bedata$GRP=="TR" & bedata$PRD==2, var]

  n1 <- length(Yi11)
  n2 <- length(Yi12)

  if(n1 * n2 < 12) {
    cat("\n Too Small Sample Size for 90% Confidence Interval !\n");
    return(NULL);
  }

  mu.r <- (mean(Yi11) + mean(Yi22)) / 2;

  G1D <- (Yi21 - Yi11) / 2
  G2D <- (Yi22 - Yi12) / 2
  D <- sort(outer(G1D, G2D, "-"));

  pval <- pwilcox(min(length(D[D>0]), length(D[D<0])), n1, n2)
  w05 <- qwilcox(0.05, n1, n2)
  w95 <- qwilcox(0.95, n1, n2)

  names(pval) <- list(c("p-value"));

  est1 <- cbind(D[w05-1], median(D), D[w95])
  est2 <- (1 + est1 / mu.r ) * 100
  est.a <- rbind(est1, est2)
  dimnames(est.a) <- list(c("90% Confidence Interval", "90% Confidence Interval(%)"), c("Lower Limit", "Point Estimate", "Upper Limit"));

#  est3 <- cbind(D[w05], median(D), D[w95+1])
#  est4 <- (1 + est3 / mu.r ) * 100
#  est.b <- rbind(est3, est4)
#  dimnames(est.b) <- list(c("90% Confidence Interval", "90% Confidence Interval(%)"), c("Lower Limit", "Point Estimate", "Upper Limit"));

#  result <- list(pval, est.a, est.b);
#  names(result) <- c("Wilcoxon Signed-Rank Test", "Hodges-Lehmann Estimate", "Hodges-Lehmann Estimate Old")

  result <- list(pval, est.a);
  names(result) <- c("Wilcoxon Signed-Rank Test", "Hodges-Lehmann Estimate")
  return(result);
}



########################################
# BE Plot


drawind <- function(g1l, g1r, g2l, g2r, g1s, g2s)
{
  for (i in 1:length(g1l)) {
    x <- jitter(c(1, 2), factor=0.3)
    y <- c(g1l[i], g1r[i])
    lines(x, y, type="l", lty=1, col="red")
    text(x[1]-0.05, y[1], paste(g1s[i]), cex=0.6, col="red")
  }

  for (i in 1:length(g2l)) {
    x <- jitter(c(1, 2), factor=0.3)
    y <- c(g2l[i], g2r[i])
    lines(x, y, type="l", lty=2, col="blue")
    text(x[2]+0.05, y[2], paste(g2s[i]), cex=0.6, col="blue")
  }
}

drawmeansd <- function(ma, sa, mb, sb, mc, sc, md, sd, y.max)
{
  sft <- 0.03
  delta <- mean(ma, mc) - mean(mb, md)
  y.RT <- mean(ma, mc) + sign(delta) * y.max * 0.05
  y.TR <- mean(mb, md) - sign(delta) * y.max * 0.05

  lines(c(1-sft, 2-sft), c(ma, mc), type="l", lty=1, col="red")
  text(1.5-sft, y.RT, "RT", col="red")
  if (sa > 0) arrows(1-sft, ma-sa, 1-sft, ma+sa, length=0.1, code=3, angle=90, col="red")
  if (sc > 0) arrows(2-sft, mc-sc, 2-sft, mc+sc, length=0.1, code=3, angle=90, col="red")

  lines(c(1+sft, 2+sft), c(mb, md), type="l", lty=2, col="blue")
  text(1.5+sft, y.TR, "TR", col="blue")
  if (sb > 0) arrows(1+sft, mb-sb, 1+sft, mb+sd, length=0.1, code=3, angle=90, col="blue")
  if (sd > 0) arrows(2+sft, md-sd, 2+sft, md+sd, length=0.1, code=3, angle=90, col="blue")
}

beplot <- function(bedata, var)
{
  if(!assert(bedata)) {
    cat("\n Bad Data Format !\n");
    return(NULL);
  }

  Si11 <- bedata[bedata$GRP=="RT" & bedata$PRD==1, "SUBJ"]
  Si21 <- bedata[bedata$GRP=="RT" & bedata$PRD==2, "SUBJ"]
  Si12 <- bedata[bedata$GRP=="TR" & bedata$PRD==1, "SUBJ"]
  Si22 <- bedata[bedata$GRP=="TR" & bedata$PRD==2, "SUBJ"]

  Yi11 <- bedata[bedata$GRP=="RT" & bedata$PRD==1, var]
  Yi21 <- bedata[bedata$GRP=="RT" & bedata$PRD==2, var]
  Yi12 <- bedata[bedata$GRP=="TR" & bedata$PRD==1, var]
  Yi22 <- bedata[bedata$GRP=="TR" & bedata$PRD==2, var]

  n1 <- length(Yi11)
  n2 <- length(Yi12)

  Y.11 <- mean(Yi11)
  Y.21 <- mean(Yi21)
  Y.12 <- mean(Yi12)
  Y.22 <- mean(Yi22)

  sY.11 <- sd(Yi11)
  sY.21 <- sd(Yi21)
  sY.12 <- sd(Yi12)
  sY.22 <- sd(Yi22)

  y.max <- max(Y.11 + sY.11, Y.21 + sY.21, Y.12 + sY.12, Y.22 + sY.22, max(bedata[,var])) * 1.2

  windows()
  par(oma=c(1,1,3,1), mfrow=c(2,2))

  plot(0, 0, type="n", ylim=c(0, y.max), xlim=c(0.5, 2.5), axes=FALSE, xlab="Period",  ylab=var, main="(a) Individual Plot for Period")
  axis(2)
  axis(1, at=c(1,2))
  drawind(Yi11, Yi21, Yi12, Yi22, Si11, Si12)

  plot(0, 0, type="n", ylim=c(0, y.max), xlim=c(0.5, 2.5), axes=FALSE, xlab="Treatment",  ylab=var, main="(b) Individual Plot for Treatment")
  axis(2)
  axis(1, at=c(1,2), labels=c("Test", "Reference"))
  drawind(Yi21, Yi11, Yi12, Yi22, Si11, Si12)

  plot(0, 0, type="n", ylim=c(0, y.max), xlim=c(0.5, 2.5), axes=FALSE, xlab="Period",  ylab=var, main="(c) Mean and SD by Period")
  axis(2)
  axis(1, at=c(1,2))
  drawmeansd(Y.11, sY.11, Y.12, sY.12, Y.21, sY.21, Y.22, sY.22, y.max)

  plot(0, 0, type="n", ylim=c(0, y.max), xlim=c(0.5, 2.5), axes=FALSE, xlab="Treatment",  ylab=var, main="(d) Mean and SD by Treatment")
  axis(2)
  axis(1, at=c(1,2), labels=c("Test", "Reference"))
  drawmeansd(Y.21, sY.21, Y.12, sY.12, Y.11, sY.11, Y.22, sY.22, y.max)

  mtext(outer=T, side=3, paste("Equivalence Plot for", var), cex=1.5)

  windows()
  par(oma=c(1,1,3,1), mfrow=c(2,2))

  boxplot(Yi11, Yi21, Yi12, Yi22, names=c("PRD=1", "PRD=2", "PRD=1", "PRD=2"), main="(a) By Sequence and Period", sub="SEQ=RT           SEQ=TR")
  boxplot(c(Yi11, Yi21), c(Yi12, Yi22), names=c("Sequence=RT", "Sequence=TR"), main="(b) By Sequence")
  boxplot(c(Yi11, Yi12), c(Yi21, Yi22), names=c("Period=1", "Period=2"), main="(c) By Period")
  boxplot(c(Yi12, Yi21), c(Yi11, Yi22), names=c("Treatment=T", "Treatment=R"), main="(d) By Treatment")
  mtext(outer=T, side=3, paste("Box Plots for", var), cex=1.5)

}
# bedata <- read.csv("d:/csv/propofolbe.csv")
# windows()
# par(mfrow=c(2,2),oma=c(1,1,3,1))
# boxplot(Cmax ~ GRP + PRD, data=bedata)
# boxplot(Cmax ~ GRP, data=bedata)
# boxplot(Cmax ~ PRD, data=bedata)
# boxplot(Cmax ~ TRT, data=bedata)

# options(digits=3)

be <- function(filename)
{
  bedata <- read.csv(filename);
# File should have the following columns
# SUBJ : Subject ID, any data type
# GRP: "RT" or "TR"
# PRD: 1 or 2
# TRT: "R" or "T"
# AUClast: numeric data type
# AUCinf: numeric data type
# Cmax: numeric data type
# Tmax: numeric data type
# Other columns as you wish

  bedata <- bedata[order(bedata$GRP, bedata$PRD, bedata$SUBJ),];
  if(!assert(bedata)) {
    cat("\n Bad Data Format !\n");
    return(NULL);
  }

  beplot(bedata, "AUClast")
  beplot(bedata, "AUCinf")
  beplot(bedata, "Cmax")
  beplot(bedata, "Tmax")

  bedata$lnAUClast <- log(bedata$AUClast);
  bedata$lnAUCinf <- log(bedata$AUCinf);
  bedata$lnCmax <- log(bedata$Cmax);

  cat("\n\n[AUClast]\n\n");
  print(betest(bedata, "lnAUClast", logtransformed=T));

  cat("\n\n[AUCinf]\n\n");
  print(betest(bedata, "lnAUCinf", logtransformed=T));

  cat("\n\n[Cmax]\n\n");
  print(betest(bedata, "lnCmax", logtransformed=T));

  cat("\n\n[Tmax]\n\n");
  print(hodges(bedata, "Tmax"));
}


kbe <- function(filename)
{
  bedata <- read.csv(filename);
# File should have the following columns
# SUBJ : Subject ID, any data type
# GRP: "RT" or "TR"
# PRD: 1 or 2
# TRT: "R" or "T"
# AUClast: numeric data type
# Cmax: numeric data type
# Tmax: numeric data type
# Other columns as you wish

  bedata <- bedata[order(bedata$GRP, bedata$PRD, bedata$SUBJ),];
  if(!assert(bedata)) {
    cat("\n Bad Data Format !\n");
    return(NULL);
  }

  beplot(bedata, "AUClast")
  beplot(bedata, "Cmax")
  beplot(bedata, "Tmax")

  bedata$lnAUClast <- log(bedata$AUClast);
  bedata$lnCmax <- log(bedata$Cmax);

  cat("\n\n[AUClast]\n\n");
  print(betest(bedata, "lnAUClast", logtransformed=T));

  cat("\n\n[Cmax]\n\n");
  print(betest(bedata, "lnCmax", logtransformed=T));

  cat("\n\n[Tmax]\n\n");
  print(hodges(bedata, "Tmax"));
}


