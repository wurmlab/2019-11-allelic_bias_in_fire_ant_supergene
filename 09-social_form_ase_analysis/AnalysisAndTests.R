# test and analysis of data

# load Figure 4 source data
load(file.choose())

# load AnalysisFunctions
source(file.choose())

# Test antilogit2 function
test1 <- all.equal(antilogit2(c(-Inf,0,1000))-c(0,0.5,1),rep(0,3))
test1

# extract logit2 values and weights
# x relative expression of locus between single and multiple-queen colonies
# y relative expression of SB/SB alleles within multiqueen colonies
# wts sample size for SB/Sb comparison as a weight

 xvals <- for_richard$lfcs_social_form
 yvals <- for_richard$lfcs_ase
 wts <- for_richard$mean_reads_ase

# antilogit transform
 transx <- antilogit2(xvals)
 transy <- antilogit2(yvals)


 mod1 <- lm(transy ~ transx, weights = wts)
 summary(mod1)
abline(mod1)

# Test deviation from (0.5,0.5)
summary(lm(I(transy-0.5) ~ I(transx - 0.5), weights = wts))


# Test for the subset with x<0.5 and x>0.5
mod2 <-  lm(transy[transx < 0.5] ~ transx[transx < 0.5], weights = wts[transx < 0.5])
summary(mod2)

# now try y vs (1-x)/x : the expected relationship if over-expression of SB is constant
mod3 <- lm(transy ~ I((1-transx)/2/transx), weights = wts)
summary(mod3)

predy <- predict(mod3, interval = 'predict', level = 0.95^(1/294))

low <- transy < predy[,'lwr']
hi <- transy > predy[,'upr']
cols <- rep('black', length(transx))
#cols[low] <- 'blue'
#cols[hi] <- 'blue'

plot(transx,transy,
     cex=sqrt(wts/mean(wts)),
     main='Allele specific expression',
     ylab=expression('Proportion of SB allele (P'[B]*')'),
     xlab=expression('Relative expression of the locus in multiple-queen colonies (P'[MQ]*')'),
     col = cols)

abline(h=0.5)
abline(v=0.5)
abline(mod1, lwd = 2)

# get sorted fitted values to plot
fy <- sort(fitted(mod3), index.return = TRUE)
# use the same sort for x values
fx <- transx[fy$ix]
# add relationship to the plot
lines(fx,fy$x, col='red', lwd = 3)

# Add the curve for Sb degeneration alone
xvals <- seq(0.2,0.9, length.out = 100)
yvals <- (xvals - 1) / -2 / xvals
lines(xvals, yvals, col='purple', lwd = 2)

# this is equivalent to x=1/(2y +1)
yvals <- seq(0.1, 0.9, length.out = 100)
xvals <- 1/(2*yvals +1)
lines(xvals, yvals, col = 'purple', lwd = 2)

mod4 <- lm(transy ~ I((transx - 1) / -2 / transx), offset = (transx - 1) / -2 / transx, weights = wts)
summary(mod4)
anova(mod4)
