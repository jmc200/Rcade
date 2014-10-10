##-----------------------------------------------------------
##Analysis with baySeq...

diffCountsBaySeq <- function(counts, targets, annoZones, cl = NULL, getLibsizesArgs = list(estimationType = "quantile", quantile = 0.75), getPriors.NBArgs = list(), getLikelihoods.NBArgs = list(), libsizes)
{

	##Check arguments--------------------------------------
	## Check targets is valid
	targets <- checkTargets(targets)

	## TODO Message explaining which targets have been found?

	##map known columns to targets columns
	map <- pmatch(c("fileid", "sampleid", "factor", "filepath"), tolower(colnames(targets)), nomatch = NA_integer_, duplicates.ok = FALSE)
	names(map) = c("fileid", "sampleid", "factor", "filepath")

	##Check targets matches coutns
	if(nrow(targets) != ncol(counts))
		stop("Targets file (", nrow(targets), " samples) does not match counts file (", ncol(counts), " samples).")
	##(TODO check colnames match (maybe just a warning?))

	##TODO check get...Args are sensible listy things

	##End check arguments----------------------------------

	##merge technical replicates together
	samplenames <- as.character(targets[,map["sampleid"]])

	counts.merged <- matrix(0, nrow(counts), length(unique(samplenames)))
	rownames(counts.merged) = rownames(counts); colnames(counts.merged) = unique(samplenames) #destinations

	for(i in 1:ncol(counts.merged))
	{
		sel <- which(samplenames %in% colnames(counts.merged)[i]) #get cols to be merged
		#sel.s <- as.character(targets[sel, map["sampleid"]])
		#print(sel.s)
		#counts.merged[,i] <- rowSums(counts[,sel])
		counts.merged[,i] <- apply(as.matrix(counts[,sel]), 1, sum)
	}

	##Sample specific targets list
	targets.reduced <- targets[,map[c("sampleid", "factor")]]
	targets.reduced <- unique(targets.reduced)

	##TODO check reduced targets against counts as a sanity check?

	##Read in to BaySeq format------------------------------------------------
	##construct count object
	NDE <- rep(1, nrow(targets.reduced))
	DE <- ifelse(tolower(targets.reduced[,2]) == "input", 1, 2) ##inputs 1, ChIP 2, chelsea 0
	replicates <- DE

	CD <- new("countData", data = counts.merged, replicates = replicates, groups = list(NDE = NDE, DE = DE), seglens = width(annoZones))
	CD@annotation <- data.frame(names = rownames(counts))

	##obtain library sizes
	getLibsizesArgs$cD = CD
	##NB: not sure if getLibsizes will return CD or not 
	temp <- do.call(getLibsizes, getLibsizesArgs) ##keep user arguments
	if("countData" %in% class(temp))
	{
		CD <- temp
	} else {
		libsizes(CD) <- temp
	}

	##get priors
	getPriors.NBArgs["cD"] = list(CD)
	getPriors.NBArgs["cl"] = list(cl)

	#stop(paste(c(getPriors.NBArgs, names(getPriors.NBArgs)), collapse=" DERP "))
	CD <- do.call(getPriors.NB, getPriors.NBArgs)

	##get likelihood
	getLikelihoods.NBArgs["cD"] = list(CD)
	getLikelihoods.NBArgs["cl"] = list(cl)
	CD <- do.call(getLikelihoods, getLikelihoods.NBArgs)

	##Spew out B values
	##
	getMA <- function (cD, group = 1, samplesA, samplesB)
	{
		if (!inherits(cD, what = "countData")) 
		    stop("variable 'cD' must be of or descend from class 'countData'")
		if (nrow(cD@posteriors) > 0) {
        tempLibs <- as.numeric(baySeq::libsizes(cD))
		    Adata <- colSums(t(cD@data[, samplesA])/tempLibs[samplesA])/length(samplesA)
		    Bdata <- colSums(t(cD@data[, samplesB])/tempLibs[samplesB])/length(samplesB)
#		    Azeros <- which(Adata == 0)
#		    Bzeros <- which(Bdata == 0)
#		    bexp <- log2(Bdata[Azeros] * mean(cD@libsizes[c(samplesA, 
#		        samplesB)]))
#		    aexp <- log2(Adata[Bzeros] * mean(cD@libsizes[c(samplesA, 
#		        samplesB)]))
#		    minZeros <- floor(min(aexp[aexp > -Inf], bexp[bexp > 
#		        -Inf]))
#		    maxZeros <- ceiling(max(aexp, bexp))
#		    logData <- log2(Adata/Bdata)
			M <- log(Adata/Bdata)
			A <- 0.5*log(Adata*Bdata)

			return(cbind(M = M, A = A))
	#        infRatios <- which(abs(logData) == Inf | is.na(logData))
	#        nonInfMax <- ceiling(max(abs(logData[-infRatios]))) + 
	#            4
	#        logData[Bzeros] <- aexp + nonInfMax - minZeros
	#        logData[Azeros] <- -bexp - nonInfMax + minZeros
	#        plot(x = logData, y = exp(cD@posteriors[, group]), ylim = c(-0.2, 
	#            1), xlim = c(-1, 1) * nonInfMax + c(-1, 1) * (maxZeros - 
	#            minZeros), ylab = "Posterior likelihood", xlab = "", 
	#            axes = FALSE, ...)
	#        abline(v = c(-nonInfMax + 3, nonInfMax - 3), col = "orange", 
	#            lty = 4)
	#        axis(side = 2, at = c(0:5/5))
	#        axis(side = 1, at = c(-maxZeros - nonInfMax, -nonInfMax, 
	#            -nonInfMax + 3, -round(nonInfMax/5), 0, round(nonInfMax/5), 
	#            nonInfMax - 3, nonInfMax, nonInfMax + maxZeros), 
	#            labels = c(maxZeros, minZeros, -Inf, -round(nonInfMax/5), 
	#                0, round(nonInfMax/5), Inf, minZeros, maxZeros))
	#        text(x = c(-maxZeros - nonInfMax, 0, nonInfMax + maxZeros), 
	#            y = -0.1, labels = c("log B", "log ratio", "log A"))
		}
		else stop("No posterior data found in 'cD' object!")
	}


	##tabulate output
	cond <- levels(targets$Condition)[i]
	t <- CD@groups$DE == 2

	baySeqOutput <- getMA(CD, group = 2, samplesA = which(t), samplesB = which(!t))
	baySeqOutput <- cbind(baySeqOutput, log.p = CD@posteriors[,"DE"], B = CD@posteriors[,"DE"] - CD@posteriors[,"NDE"])
	as.data.frame(baySeqOutput)

	#list(out = as.data.frame(baySeqOutput), prior = ...) #TODO bring forward priors from baySeq?
}
