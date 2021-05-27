# Code by Camille Coux
# calculation of interaction preferences
# input: interaction frequencies and relative abundances. 
# frugiv_data.csv can be found in the online supplementary files, 
# or on Figshre (DOI: 10.6084/m9.figshare.12706598 ) 

library(readr)
library(plyr)
library(pscl)
library(MASS)


# Load data
df <- read.csv2("frugivory_data.csv", dec=".")

# check import went well
str(df)

# in each site, run neg binom model and exctract the count intercept and theta
ZInb_out <- ddply(df, .(Location), function(x){
  # run neutral model
  ZInb <- zeroinfl(freq ~ 1 + offset(log(row.abun)+log(col.abun)) |1, 
                   dist="negbin", link="logit", data=x)
  # run null model
  ZInb_null <- zeroinfl(freq ~ 1, dist="negbin", link="logit", data=x)
  
  return(data.frame(count.intercept = ZInb$coefficients$count, 
                    theta = ZInb$theta, 
                    zero=ZInb$coefficients$zero,
                    null.intercept=ZInb_null$coefficients$count, 
                    null_theta=ZInb_null$theta, 
                    null_zero=ZInb_null$coefficients$zero))
})


# This gets some additional information from the ZInb that is useful for plotting things
df <- ddply(df, .(Location), function(x){
  ZInb <- zeroinfl(freq ~ 1 + offset(log(row.abun)+log(col.abun)) |1, 
                   dist="negbin", link="logit",data=x) 
  x$predicted <- predict(ZInb)
  x$fittedzero <- predict(ZInb,type="zero")
  x$fittedcount <- predict(ZInb,type="count")
  x$nbzero <- dnbinom(0,mu=x$fittedcount,size=ZInb$theta)
  x$pearson <- residuals(ZInb,type="pearson")
  return(x)
})
df <- merge(df, ZInb_out)


# the values predicted by just the negbin part of the ZInb_null
# note you can also get this from predict(ZInb, type="count")
df$fittednull <- df$row.abun*df$col.abun*exp(df$null.intercept)

# the residuals compared to just the negbin part
# note that you can get this from residuals(ZInb, type="response")
df$resids.nb <- df$freq - df$fittedcount
df$resids.null <- df$freq - df$fittednull

# the deviance residuals from just the nb part
# this is where the theta values extracted earlier are used
df$dev.resids.nb <- 
  negative.binomial(theta=df$theta)$dev.resids(
    df$freq, df$fittedcount, rep(1,nrow(df)))
df$dev.resids.null <- 
  negative.binomial(theta=df$null_theta)$dev.resids(
    df$freq, df$fittednull, rep(1,nrow(df)))

# interaction deviances used as response variable at the interaction level:
df$log.dev.resids.nb <- log(df$dev.resids.nb)


# preferences and avoidances
df$prefs.and.avoids <- df$dev.resids.nb * sign(df$resids.nb)
df$pref.or.avoid <- "avoid"
df$pref.or.avoid[which(sign(df$resids.nb)>0)] <- "pref"
df$pref.or.avoid <- as.factor(df$pref.or.avoid)


# SITE LEVEL    
# calculating the total deviance for each site to get the neutrality ratio
DEVnull <- tapply(df$dev.resids.null, df$Location, sum)
DEVnb <- tapply(df$dev.resids.nb, df$Location, sum) 

n.ratio <-  1-(DEVnb/DEVnull) %>% as.data.frame
n.ratio$Location <- row.names(n.ratio)

# merge with df:
colnames(n.ratio) <- c("n.ratio", "Location")
df <- merge(df, n.ratio, by="Location")



