args <- commandArgs(trailingOnly=TRUE)

pluriTestNoLabelsCommandLine <- function (NewDataFileName, DataRepository) 
{
    require(lumi)
    require(xtable)
    require(GO.db)

    load (file=DataRepository)    
    sink(file = "pluritest_log.txt")
    working.lumi <- lumiR(NewDataFileName, convertNuID = FALSE, 
        inputAnnotation = FALSE)
    fData(working.lumi)[, 1] <- gsub("\"", "", fData(working.lumi)[, 1])
    pdf(file = "pluritest_image01.pdf")
    plot(working.lumi, what = "boxplot")
    dev.off()
    working.lumi <- lumiT(working.lumi)
    hc <- hclust(as.dist(1 - abs(cor(exprs(working.lumi[, ])))))
    pdf(file = "pluritest_image02a.pdf")
    plot(hc, hang = -1, main = "Clustering of vst-transformed samples", 
        sub = "distance based on pearson correlations", xlab = "")
    dev.off()
    A <- fData(working.lumi)[, 1]
    B <- fData(H9targetArray)[, 1]

    sel.match <- match(B, A) #find where Bs are in A(pluritest) # sel.match = 41751 16775 41089

    C <- working.lumi[na.omit(sel.match), ]

    working.lumi <- lumiN(C, method = "rsn") ### INCOMPLETE - NEED TO FIX

    A <- fData(working.lumi)[, 1]
    try({
        sel <- match(rownames(W15), fData(working.lumi)[, 1]) # 16657 is an NA


        coef <- c(-126.7095, 0.004567437, 0.004377068, 0.001043193)

        workingLumiW15matches <- na.omit(sel) 
        length(workingLumiW15matches) #48743
        wlExprs <- exprs(working.lumi[workingLumiW15matches])

        W15Selections <- W15[!is.na(sel), ]
        W12Selections <- W12[!is.na(sel), ]

        H15.new <- predictH(wlExprs, W15Selections)
        H12.new <- predictH(wlExprs, W12Selections)
        rss.new <- apply((wlExprs                                   - W12Selections %*% H12.new)^2, 2, sum)
        RMSE.new <- sqrt(rss.new/sum(!is.na(sel)))
        novel.new <- apply((wlExprs                                   - W12Selections  %*% H12.new)^8, 2, sum)
        novel.new <- (novel.new/sum(!is.na(sel)))^(1/8)
        s.new <- drop(coef[1] + coef[2:4] %*% H15.new[c(1, 14, 13), ])
    })

    table.results <- matrix(, nrow = ncol(exprs(working.lumi)), ncol = 5)
    rownames(table.results) <- colnames(exprs(working.lumi))
    colnames(table.results) <- c("pluri-raw", "pluri logit-p", 
        "novelty", "novelty logit-p", "RMSD")
    try({
        print(s.new)
        pdf(file = "pluritest_image02.pdf")
        par(mar = c(12, 4, 4, 2))
        par(xaxt = "n")
        plot(s.new, main = "pluripotency", xlab = "", ylab = "pluripotency", 
            ylim = c(-130, 70))
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
        pdf(file = "pluritest_image03_no_labels.pdf")
        color.palette = colorRampPalette(c("red", "pink1", "aliceblue", 
            "lightblue", "blue"), bias = 1)
        filled.contour2(y = c(-129:70), x = c((1:200)/50), background129_70x1_4, 
            col = colram(50), nlevels = 35, xlab = "novelty", 
            ylab = "pluripotency")
        points(s.new ~ novel.new, cex = 0.4, main = "Overview")
        dev.off()
    })
    try({
        palette(colorRampPalette(c("green", "orange", "orange", 
            "orange", "red"))(5))
        df.novelty.new <- data.frame(novelty = novel.new)
        pdf(file = "pluritest_image03c.pdf")
        par(mar = c(12, 4, 4, 2))
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

pluriTestNoLabelsCommandLine(args[1], args[2])
