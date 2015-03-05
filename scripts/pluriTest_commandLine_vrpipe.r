# RUN ON THE COMMANDLINE USING:
# /software/bin/R-2.15 --slave --args "fileName" ".RData location on lustre" < pluriTestCommandLine.r 
# e.g. /software/bin/R-2.15 --slave --args "pilot_data.txt-EDIT-7-names" "/lustre/scratch105/vrpipe/refs/pluritest.RData"

args <- commandArgs(trailingOnly=TRUE)

#GENERATE GRAPHS AND CSV FILE GIVING PLURIPOTENCY, NOVELTY AND COMBINED SCORE
#NEEDS R VERSION 2.1.5 AND BIOCONDUCTOR R PACKAGE INSTALLED
pluriTestCommandLine <- function (NewDataFileName, DataRepository) 
{
    require(lumi) # load the Bioconductor package
    require(xtable)
    require(GO.db)
    require (calibrate)

    load (file=DataRepository)
    sink(file = "pluritest_log.txt") # print to the log file
    working.lumi <- lumiR(NewDataFileName, convertNuID = FALSE, inputAnnotation = FALSE) # read in the gene expression data
    fData(working.lumi)[, 1] <- gsub("\"", "", fData(working.lumi)[, 1])
    png(file = "pluritest_image01.png") # Boxplot of microarray intensity
    plot(working.lumi, what = "boxplot")
    dev.off()
    working.lumi <- lumiT(working.lumi) # Transform the data
    
    png(file = "pluritest_image02a.png") # Clustering of vst-transformed samples (VST = variance stabilization normalisation)
    if (length(as.dist(1 - abs(cor(exprs(working.lumi[, ])))))!=1){
    hc <- hclust(as.dist(1 - abs(cor(exprs(working.lumi[, ])))))
    plot(hc, hang = -1, main = "Clustering of vst-transformed samples", 
        sub = "distance based on pearson correlations", xlab = "")
    } else {plot(0,type='n',axes=FALSE,ann=FALSE)}
    dev.off()

# RSN (robust spline normalisation) NORMALISATION OF THE DATA

    A <- fData(working.lumi)[, 1] #matches on ILMN_Ids for lumi/RSN
    B <- fData(H9targetArray)[, 1] # the well studied H9 stem cell line
    sel.match <- match(B, A) # check input IDs against pluritest IDs from the well studied H9 stem cell line

    C <- working.lumi[na.omit(sel.match), ]
    working.lumi <- lumiN(C, method = "rsn") # normalise data ### INCOMPLETE - NEED TO FIX
    #working.lumi <- lumiN(working.lumi, method = "rsn", target = H9targetArray[is.na(sel.match) == FALSE, ])

    A <- fData(working.lumi)[, 1]
    sel.match <- match(colnames(W15), A)

# CALCULATIONS

    try({
        sel <- match(rownames(W15), fData(working.lumi)[, 1])
        coef <- c(-126.7095, 0.004567437, 0.004377068, 0.001043193)

        workingLumiW15matches <- na.omit(sel) 
        wlExprs <- exprs(working.lumi[workingLumiW15matches])
        W15Selections <- W15[!is.na(sel), ]
        W12Selections <- W12[!is.na(sel), ]

        H15.new <- predictH(wlExprs, W15Selections) #H matrix prediction script runs the Lee & Seung multiplicative update for a maximum of
        H12.new <- predictH(wlExprs, W12Selections) #2000 iterations, checking convergence every 20 steps                                  
        rss.new <- apply((wlExprs - W12Selections %*% H12.new)^2, 2, sum)

        RMSE.new <- sqrt(rss.new/sum(!is.na(sel)))
        novel.new <- apply((wlExprs - W12Selections  %*% H12.new)^8, 2, sum)
        novel.new <- (novel.new/sum(!is.na(sel)))^(1/8)
        s.new <- drop(coef[1] + coef[2:4] %*% H15.new[c(1, 14, 13), ])
    })
    
    #re-order the scores by sample names
    ctrls <- grep("_CTRL",names(s.new),value=FALSE)
    if (length(ctrls)!=0){	
    s.new.ips <- s.new[-ctrls]
   	s.new <- c(s.new[ctrls],s.new.ips[order(names(s.new.ips))])
    }
    ctrls <- grep("_CTRL",names(novel.new),value=FALSE)
    if (length(ctrls)!=0){
    novel.new.ips <- novel.new[-ctrls]
    novel.new <- c(novel.new[ctrls],novel.new.ips[order(names(novel.new.ips))])
    }

    table.results <- matrix(, nrow = length(s.new), ncol = 5)
    rownames(table.results) <- names(s.new)
    colnames(table.results) <- c("pluri-raw", "pluri logit-p", 
        "novelty", "novelty logit-p", "RMSD")
    try({
        png(file = "pluritest_image02.png") # Graph of pluripotency scores
        par(mar = c(18, 4, 4, 2))
        par(xaxt = "n")
        plot(s.new, main = "pluripotency", xlab = "", ylab = "pluripotency", 
            ylim = c(-130, 70), cex=1, pch=23, bg="black")
        abline(h = 25.25, lty = "dashed", col = "red")
        abline(h = 59.95, lty = "dashed", col = "red")
        abline(h = -28.92, lty = "dashed", col = "lightblue")
        abline(h = -130, lty = "dashed", col = "lightblue")
        par(xaxt = "s")
        axis(1, at = c(1:length(s.new)), labels = names(s.new), 
            las = 2)
        dev.off()
    })
    table.results[, 1] <- round(s.new, 3)
    table.results[, 2] <- round(exp(s.new)/(1 + exp(s.new)), 
        3)
    table.results[, 3] <- round(novel.new, 3)
    table.results[, 5] <- round(RMSE.new, 3)
    try({
        png(file = "pluritest_image03.png") # Graph of pluripotency versus novelty
        color.palette = colorRampPalette(c("red", "pink1", "aliceblue", 
            "lightblue", "blue"), bias = 1)
        filled.contour2(y = c(-129:70), x = c((1:200)/50), background129_70x1_4, 
            col = colram(50), nlevels = 35, xlab = "novelty", 
            ylab = "pluripotency")

        sNames <- names(s.new)
        novelNames <- names(novel.new)

        ctrls <-grep("_CTRL",names(s.new),value=FALSE)
        if (length(ctrls)!=0){
        s.new.ctrls<-s.new[ctrls]
        s.new.stcls<-s.new[-ctrls]
        novel.new.ctrls<-novel.new[ctrls]
        novel.new.stcls<-novel.new[-ctrls]
        points(s.new.stcls ~ novel.new.stcls, cex = 1, main = "Overview", col="red", pch=20)
        points(s.new.ctrls ~ novel.new.ctrls, cex = 1, main = "Overview", pch=20)
        }else{points(s.new ~ novel.new, cex = 1, main = "Overview", col="red", pch=20)}

        #textxy(as.vector(novel.new[novelNames])-0.2, as.vector(s.new[sNames])-0.2, labs=sNames, cx=1) # MIGHT NEED TO LOWER cx BACK TO 0.25 FOR SMALLER FONT
        #points(s.new ~ novel.new, cex = 0.75, pch=21, main = "Overview", bg="black")
        dev.off()
    })
    try({
        palette(colorRampPalette(c("green", "orange", "orange", 
            "orange", "red"))(5))
        df.novelty.new <- data.frame(novelty = novel.new)
        png(file = "pluritest_image03c.png") # Graph of novelty scores
        par(mar = c(18, 4, 4, 2))
        par(xaxt = "n")
        barplot(novel.new, col = pmin(5, 10 * predict(logit.novelty, 
            type = "response", newdata = df.novelty.new) + 1), 
            names.arg = c(1:length(novel.new)), xlab = "", xlim = c(0, 
                length(novel.new)), width = 0.9, space = (1/9), 
            ylim = c(0, 4))
        title(main = "novelty")
        par(xaxt = "s")
        axis(1, at = c(1:nrow(table.results)) - 0.4, labels = names(s.new), 
            las = 2)
        table.results[, 4] <- round(predict(logit.novelty, type = "response", 
            newdata = df.novelty.new), 3)
        dev.off()
    })
    table.results[, 5] <- round(RMSE.new, 3)
    write.csv(table.results, file = "pluritest.csv")
    sink()
}

pluriTestCommandLine(args[1], args[2])
