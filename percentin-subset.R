# 2016.05.31
# Generate raw data.frame and subset using %in%
# http://www.ats.ucla.edu/stat/r/modules/subsetting.htm
A=data.frame(SUBJID=c("A","B","C"),CONC=c(5,10,7),ORDER=c(1,2,3))
SUB=c("A","C")
SEL=A[A$SUBJID%in%SUB, ]
B=subset(SEL,select=c(SUBJID,ORDER))
B
summary(B)

# sapply single function (NOT RECOMMENDED)
sapply(A,FUN=c("mean","median"))
sapply(A[,2:3], median, na.rm = TRUE)
sapply(A[,2:3], mean, na.rm = TRUE)

# sapply single function (RECOMMENDED)
# http://www.r-bloggers.com/applying-multiple-functions-to-data-frame/
multi.fun <- function(x) {
      c(min = min(x), mean = mean(x), max = max(x), sd=sd(x))
}

# Make SUBJID=NA
B=data.frame(SUBJID=NA,A[,2:3])

# sapply single function
# sapply(cars, multi.fun) # example 
STAT=sapply(B, multi.fun)

# rbind both data.frame
rbind(A,STAT)

