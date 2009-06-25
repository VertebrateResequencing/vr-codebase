# R functions for 1000genomes data analysis
# Aylwyn Scally 2008

library(Biostrings)

read.depth = function(datafile, type, N = NULL, addcols = c(), addclasses = c()){ 
    #print(paste('Reading data from', datafile, '...'))
	if (type == 'bin'){
		cnames = c('chr', 'pos', 'dep', 'ref')
		cclasses = c('factor', 'integer', 'integer', 'numeric')
	}
	if (type == 'bingcn'){
		cnames = c('chr', 'pos', 'dep', 'ref', 'gc', 'n')
		cclasses = c('factor', 'integer', 'integer', 'numeric', 'integer', 'integer')
	}
	if (type == 'gc'){
		cnames = c('chr', 'pos', 'gc', 'n')
		cclasses = c('factor', 'integer', 'integer', 'integer')
	}
	if (type == 'samp'){
		cnames = c('chr', 'pos', 'dep', 'ref', 'gc', 'n')
		cclasses = c('factor', 'integer', 'integer', 'numeric', 'numeric', 'numeric')
	}

	if (length(addcols) > 0){
		cnames = c(cnames, addcols)
		cclasses = c(cclasses, addclasses)
	}

	if (!is.null(N)){
		dat = read.table(datafile, col.names = cnames, nrows = N, colClasses = cclasses)
		bin = NULL
	} else {
		if (length(grep('#', readLines(datafile, n=1))) > 0){
			N = as.numeric(scan(datafile, nlines = 1, what = 'character')[2])
			bin = as.numeric(scan(datafile, skip = 1, nlines = 1, what = 'character')[2])
			dat = read.table(datafile, col.names = cnames, nrows = N, colClasses = cclasses)
		}else{
			dat = read.table(datafile, col.names= cnames, colClasses = cclasses)
			bin = NULL
		}
	}

	return(list(data = dat, binsize = bin))
}

read.depdat = function(analfile, gcfile=NULL, Ndat = NULL, bin = NULL, gcmax = NULL){
	if (is.null(gcfile)){
		depdat = read.depth(analfile, 'samp', Ndat)
		gcmax = 1.0
	}else{
		depdat = read.depth(analfile, 'bin', Ndat)
		gcdat = read.depth(gcfile, 'gc', Ndat)
		depdat$data = cbind(depdat$data, gc=gcdat$data$gc, n=gcdat$data$n)
	}
	if (!is.null(bin))
		depdat$bin = bin
	if (is.null(gcmax))
		gcmax = depdat$bin
	depdat$gcmax = gcmax

	return(depdat)
}


depth.nbfit = function(depth.seq, ref.seq = c(), plot = TRUE, qtrim = 0, plotbin = 1, count.zero = FALSE, fit.method = 'mm', plotpois = TRUE, plotgam = TRUE, setpar = TRUE, annotate.fit = TRUE, ...){
	if (plot & setpar){
		npanels = as.numeric(plot) + as.numeric(plotgam)
		par(mfrow = c(1, npanels), xpd = FALSE, mar = c(4, 4, 2, 2)+0.5, mex = 0.7, font.main = 1)
	}

	if (length(depth.seq) == 0){
	    fit = list(mu = NA, size = NA, var = NA)
		if (plot){
			plot(c(0, xmax), c(0, 1), type = 'n', bty = 'n', xlab = 'mapped depth', ylab = 'Frequency', ...)
			#axis(1, at = seq(0, xmax, 5)/depscale, labels=seq(0, xmax, 5))
		}
		return(fit)
	}

	if (length(ref.seq > 0)){
		if (count.zero)
			ref.seq[depth.seq == 0] = 1 # otherwise these will be NA
		depth.seq = depth.seq[round(ref.seq) == 1]
	}

	depth.seq = na.omit(depth.seq)
	N = length(depth.seq)

	mu = mean(depth.seq)

	#lim = qpois(c(qtrim, 1 - qtrim), mu)
	lim = round(quantile(depth.seq, c(qtrim, 1 - qtrim)))
	depth.seq = depth.seq[depth.seq >= lim[1] & depth.seq <= lim[2]]

	depth.hist = hist(depth.seq, breaks = seq((lim[1] - 0.5), (lim[2] + 0.5), plotbin), plot = FALSE)

	mu = mean(depth.seq)
	dvar = var(depth.seq)

	if (plot){
		plot(depth.hist, xlab = 'mapped depth', col = 'gray', border = 'white', ...)
		#axis(1, at = seq(0, xmax, 5)/depscale, labels=seq(0, xmax, 5))
	}
		
	if (plot & plotpois)
    	lines(depth.hist$mids, N*dpois(depth.hist$mids, mu), lty = 2)
    	#lines(depth.hist$mids, dpois(round((depth.hist$mids)/poisson.scale), mu/poisson.scale)/poisson.scale, lty = 2)

	if (fit.method == 'ml'){ # numerical max likelihood fit
		est = fitdistr(depth.seq, 'negative binomial')$estimate
		fsize = est[1]
		fmu = est[2]
	}
	if (fit.method == 'mm'){ # method of moments
		fmu = mu
		fsize = mu^2 / (dvar - mu)
	}
	fit = list(mu = fmu, size = fsize, var = dvar)

	if (plot){
		lines(depth.hist$mids, N*dnbinom(depth.hist$mids, size=fsize, mu=fmu))
		if (annotate.fit)
			#legend('topright', sprintf('avg depth: %.2f', fmu), bty='n')
			#legend('topright', c(sprintf('mu: %.2f', fmu), sprintf('var/mu: %.2f', dvar/mu)), bty='n')
			legend('topright', c(sprintf('mu: %.2f', fmu), sprintf('size: %.2f', fsize)), bty='n')
	}

	if (plot & plotgam){
		x = seq(0, 2, 0.01)
		plot(x, dgamma(x, fit$size, fit$size), type = 'l', ylab = 'Frequency', main = 'Mixture Gamma')
	}

	return(fit)
}


plot.depth = function(depdat, reg = c(), chr = NA, plot.ref = TRUE, dep.max = NA, ref.max = NA, name = '', mark.reg = c(), mark.loc = c()){
	if (is.na(chr))
		chr = levels(depdat$data$chr)[1]
	dat = depdat$data[depdat$data$chr == chr, ]
	if (length(reg) > 0) # select region
		dat = dat[which(dat$pos >= reg[1] & dat$pos <= reg[2]), ]
	else{
		N = nrow(dat)
		reg = c(dat$pos[1], dat$pos[N])
	}

	nplots = 1 + plot.ref
	layout(matrix(1:nplots), heights = rep(1, nplots))
	par(mar = c(2, 4, 2, 2) + 0.5, mex = 0.7, tck = 0.02)
	par(yaxt = 's', xaxt = 's', xpd = FALSE, cex = 1.0, bty = 'o')

	if (is.na(dep.max))
		dep.max = 3 * mean(dat$dep)
	if (is.na(ref.max))
		ref.max = max(dat$ref)

	plot(reg, c(0, dep.max), main = name, xlab = '', ylab = 'mapped depth', type = 'n')
	if (length(mark.reg) > 0)
		for(i in 1:nrow(mark.reg))
			rect(mark.reg[i, 1], 0, sum(mark.reg[i, ]), dep.max, col = 'gray', border = FALSE)
	lines(dat$pos, dat$dep)
	
	if (plot.ref){
		par(mar = c(2, 4, 0, 2) + 0.5, mex = 0.7, tck = 0.02)
		plot(reg, c(0, ref.max), main = '', xlab = '', ylab = 'ref CN', type = 'n')
		if (length(mark.reg) > 0)
			for(i in 1:nrow(mark.reg))
				rect(mark.reg[i, 1], 0, sum(mark.reg[i, ]), ref.max, col = 'gray', border = FALSE)
		lines(dat$pos, dat$ref)
	}
}


gcdepth = function(depdat, sname = '', depmax = NULL, hc = FALSE, nbins = 30, quant = 0.1, qtrim = 0.01, binned = FALSE){
	if (hc){
		if (nchar(sname) > 0){
			pdf(file = paste(sname, 'gcdepth.pdf', sep = '-'), width = 6, height = 6)
			sink(paste(sname, 'gcdepth.txt', sep = '-'))
		}else{
			pdf(file = 'gcdepth.pdf', width = 6, height = 6)
			sink("gcdepth.txt")  
		}
		#png(file = 'gcdepth.png', width = 600, height = 600)
	}
	cat(sname, '\n')
	cat('#data.binsize', 'n.gcbins', 'quantile.plot', 'quantile.trim\n', sep='\t')
	cat(depdat$bin, nbins, quant, qtrim, sep='\t')
	cat('\n')
	if (length(levels(depdat$data$chr)) <= 30)
		cat('# sampled from chr:', levels(depdat$data$chr), '\n')

	depdat$data = na.omit(depdat$data)
	Nall = nrow(depdat$data)

	meantot.all = mean(depdat$data$dep)
	sstot.all = sum((depdat$data$dep - meantot.all)^2)
	msstot.all = Nall * meantot.all

# get unique sequence
	udata = depdat$data[depdat$data$ref == 1,]
	#udata = depdat$data[round(depdat$data$ref) == 1,]
	Nu = nrow(udata)
	#print(c(Nall, Nu, 100 * Nu/Nall))
	udata = udata[udata$n <= 0.1,]
	Nu = nrow(udata)
	#print(c(Nall, Nu, 100 * Nu/Nall))

	meantot = mean(udata$dep)
	sstot = sum((udata$dep - meantot)^2)
	msstot = Nu * meantot

	mu = rep(0, nbins)
	avgc = rep(0, nbins)
	var = rep(0, nbins)
	lquant = rep(0, nbins)
	uquant = rep(0, nbins)
	gfrac = rep(0, nbins)
	ssres = 0
	mssres = 0
	cat('#gc.low', 'gc.high', 'gc.avg', 'gc.percent', 'N', 'mean', 'low.quant', 'high.quant', 'var', sep='\t')
	cat('\n')
	for (ib in 1:nbins){
		binlow = depdat$gcmax*(ib - 1) / nbins
		binhigh = depdat$gcmax*ib / nbins
		if (ib == 1){
			#binlow = -1
			bindata = udata$dep[udata$gc >= 0 & udata$gc <= binhigh]
		}else{
			bindata = udata$dep[udata$gc > binlow & udata$gc <= binhigh]
		}
		N = length(bindata)
		if (ib == 1)
			gfrac[ib] = 0.5 * 100*N/Nu
		if (ib > 1)
			gfrac[ib] = gfrac[ib-1]+ 100*N/Nu
		fit = depth.nbfit(bindata, qtrim = 0.01, setpar=FALSE, plot=FALSE)
		mu[ib] = fit$mu
		avgc[ib] = mean(c(binlow, binhigh))/depdat$gcmax
		var[ib] = fit$var
		ssres = ssres + sum((bindata - mean(bindata))^2)
		if (N > 0)
			mssres = mssres + N * mean(bindata)

		lim = round(quantile(bindata, c(qtrim, 1 - qtrim)))
		bindata = bindata[bindata >= lim[1] & bindata <= lim[2]]
		lquant[ib] = quantile(bindata, quant)
		uquant[ib] = quantile(bindata, 1 - quant)

		if (depdat$gcmax == 1)
			cat(100* binlow, 100* binhigh, avgc[ib]*depdat$gcmax, gfrac[ib], N, mu[ib], lquant[ib], uquant[ib], var[ib], sep='\t')
			#cat(floor(100* binlow +1 ), floor(100* binhigh), avgc[ib]*depdat$gcmax, gfrac[ib], N, mu[ib], lquant[ib], uquant[ib], var[ib], sep='\t')
		else
			cat(floor(binlow +1 ), floor(binhigh), avgc[ib]*depdat$gcmax, gfrac[ib], N, mu[ib], lquant[ib], uquant[ib], var[ib], sep='\t')
		cat('\n')
	}

	ssdf = data.frame(mean = c(meantot.all, meantot, 0), N = c(Nall, Nu, Nu), SS = c(sstot.all, sstot, ssres), Pois.SS = c(msstot.all, msstot, mssres), row.names = c('all', 'unique', 'residual'))
	ssdf = cbind(ssdf, overdisp=ssdf$SS / ssdf$Pois.SS)
	#print(c(sstot, msstot, msstot/sstot))
	#print(c(ssres, mssres, mssres/ssres))
	#print(c(sstot, ssres, 1 - ssres/sstot))
	cat('# model summary\n')
	print(ssdf)

	pd = data.frame(gfrac, avgc, mu, var, lquant, uquant)
	pd = na.omit(pd)
	nb = nrow(pd)

	par(mfrow = c(1, 1), xpd = FALSE, mar = c(4, 4, 4, 3)+0.5, mex = 0.7, font.main = 1)

	if (is.null(depmax))
		depmax = 1.3 * max(pd$uquant)
	if (binned)
		plot(c(0,100), c(0,depmax), type='n', lty=2, xlab='percentile of unique sequence ordered by GC content', ylab=paste('mapped depth (bin size', depdat$bin, 'bp)'), main = '') 
	else
		plot(c(0,100), c(0,depmax), type='n', lty=2, xlab='percentile of unique sequence ordered by GC content', ylab=paste('mapped depth'), main = '') 
	polygon(c(pd$gfrac, rev(pd$gfrac)), c(pd$uquant, rev(pd$lquant)), col='grey',border=NA) 
# GC % lines
	gclocs = approx(pd$avgc, pd$gfrac, seq(0,1.0,0.1))$y
	abline(v = gclocs, col='black', lty = 3)
	axis(3, at = gclocs, labels = 100*seq(0,1.0,0.1))
	mtext('GC content (%)', side=3, line = 3)
	mtext(sname, side=3, line = 3, adj=0)

	lines(pd$gfrac, pd$mu)
	lines(pd$gfrac, qpois(quant, pd$mu), lty=2) 
	lines(pd$gfrac, qpois(1 - quant, pd$mu), lty=2) 

	if (hc){
		dev.off()
		sink()  
	}
}

#Ndat = 6e6
#gcmax = 100
#indiv = 'Trio-CEU/NA12878'
#analfile = '/lustre/scratch1/dmc/g1k/v1/Trio-CEU/NA12878/aln.bindepth-100'
#gcfile = '/lustre/scratch1/dmc/g1k/ref/human_b36_female.gcbin-100'
#
#depdat = read.depdat(analfile, gcfile, Ndat = 6e6, bin = 100)
#gcdepth(depdat, sname = indiv, depmax = 60, hc = TRUE, nbins = 30)
