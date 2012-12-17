# .Rprofile -- commands to execute at the beginning of each R session
#
# You can use this file to load packages, set options, etc.
#
# NOTE: changes in this file won't be reflected until after you quit
# and start a new session
#made by ktr331
qpcr.input<-function(x,save){
#x <- 'C:/Users/neko/Desktop/2012Jan30_results.txt';
#save <- 'C:/Users/neko/Desktop/test';
qpcr.save <- as.character(save)
qpcr.table <- readLines(x)
qpcr.count <- suppressWarnings(grep("Well", qpcr.table) - 1)
abiv <- grep("^.*\\.txt$", x, ignore.case = FALSE)
if (length(abiv) != 0){
qpcr.rtx <- read.delim(x, skip = qpcr.count, blank.lines.skip = FALSE, fill = T)
#save(file="C:/Users/myat/Desktop/test/debug.txt", names(qpcr.rtx), "\n")
#names(qpcr.rtx) <- sub("CÑ.", "Ct", names(qpcr.rtx))
#it is dependent by os(windows or unix)
names(qpcr.rtx) <- gsub("Quantity", "Qty", names(qpcr.rtx))
names(qpcr.rtx) <- gsub("Target.Name", "Detector", names(qpcr.rtx))
}else{
qpcr.rtx <- read.csv(x, skip = qpcr.count, blank.lines.skip = FALSE,  fill = T)
}
#this point is needed by ABI analyse soft version 2
qpcr.rt <- qpcr.rtx[, c("Sample.Name", "Qty", "Task", "Ct", "Detector")]
qpcr.rt$Task <- gsub("UNKNOWN", "Unknown", qpcr.rt$Task)
qpcr.rt$Task <- gsub("STANDARD", "Standard", qpcr.rt$Task)
qpcr.rt$Ct <- gsub("Undetermined", NA, qpcr.rt$Ct)
qpcr.rt$Ct <- as.numeric(qpcr.rt$Ct)
qpcr.rt$Qty <- as.numeric(gsub("Undetermined", NA, qpcr.rt$Qty))
rt.cal <- qpcr.rt[qpcr.rt$Task=="Unknown", ]
#rt.cal <- qpcr.rt[qpcr.rt$Task=="NTC", ]
rt.std <- na.omit(qpcr.rt[qpcr.rt$Task=="Standard", ])
test <- unique(rt.std$Detector)
i <- length(test)
while (i >= 1){
  y <- test[i]
  rt.std.tg <-rt.std[rt.std$Detector==y, ]
  min.qty <- min(rt.std.tg$Qty)
  min.ct <- min(rt.std.tg$Ct)
  max.qty <- max(rt.std.tg$Qty)
  max.ct <- max(rt.std.tg$Ct)
  X<- as.numeric(rt.std.tg$Qty)
  Y<- as.numeric(rt.std.tg$Ct)
  res.nls <- nls(Y ~ a*log10(X)+b, start=c(a=1,b=1),trace=TRUE)
  cal.nls <- summary(res.nls)$parameters[, "Estimate"]
  pdt.nls <- predict(res.nls)
  cal.nls["a"] -> a
  cal.nls["b"] -> b
  if(qpcr.save != ""){
    pdf(file=paste(qpcr.save, "/std_", test[i], "_", paste(strsplit(as.character(Sys.time()), " |:|-")[[1]], sep="", collapse = "_"), ".pdf",sep = ""))
    plot(rt.std.tg$Qty,rt.std.tg$Ct,main=paste("Standerd Line of ", y), ylab="Ct", xlab=expression(paste(log[10], " Qty")), log="x")
    curve(a*log10(x)+b, lty="dashed", col="red", log="x", add = T)
    dev.off()   
    postscript(paste(qpcr.save, "/std_", test[i], "_", paste(strsplit(as.character(Sys.time()), " |:|-")[[1]], sep="", collapse = "_"), ".eps",sep = ""), horizontal=FALSE)
    plot(rt.std.tg$Qty,rt.std.tg$Ct,main=paste("Standerd Line of ", y), ylab="Ct", xlab=expression(paste(log[10], " Qty")), log="x")
    curve(a*log10(x)+b, lty="dashed", col="red", log="x", add = T)
    dev.off()   
    }
  rt.cal.tg <-rt.cal[rt.cal$Detector==y, ]
  rt.cal.tg <-na.omit(rt.cal.tg[,c("Sample.Name","Detector","Qty","Ct")])
  rt.cal.tg$Qty <- NA  #NULL
	rt.cal.tg$Qty <- 10^((as.numeric(rt.cal.tg$Ct)-b)/a)
   pcr.effi <- (10^(-1/a)-1)*100
  rt.cal.tg$PCR.effi <- pcr.effi
  qpcr.uq <- unique(rt.cal.tg$Sample.Name)
  j <- length(qpcr.uq)
  honmon2 <- rt.cal.tg
  while (j >= 1){
    honmon <- rt.cal.tg[rt.cal.tg$Sample.Name==qpcr.uq[j], ]
    if (j == length(qpcr.uq)){
  honmon2$Qty.Mean <- NA
  honmon2$Qty.SD <- NA
    }
  k <- nrow(honmon)
  if (k >= 2){
    honmon$Qty.Mean <- mean(honmon$Qty)
    honmon$Qty.SD <- sd(honmon$Qty)
    
    }else{
    honmon$Qty.Mean <- honmon$Qty
    honmon$Qty.SD <- "SimpleDATA"
    }
  honmon2 <-  rbind(honmon, honmon2)
  j <- j - 1
  }
rt.cal.tg <-  na.omit(honmon2)
#rt.cal.tg$Mean <- gsub("SimpleDATA", NA, rt.cal.tg$Mean)
#rt.cal.tg$SD <- gsub("SimpleDATA", NA, rt.cal.tg$SD)
##########DATA UNIQUE METHOD##########
rt.cal.tg <- rt.cal.tg[, -which (colnames(rt.cal.tg) %in% c("Qty", "Ct"))]
rt.cal.tg <- unique(rt.cal.tg[order(rt.cal.tg$Detector, rt.cal.tg$Sample.Name), ])
##########DATA UNIQUE METHOD##########
  rownames(rt.cal.tg) <- 1:nrow(rt.cal.tg)
	if (i == length(test)){
    qpcr.result <- as.list(rep(NA,i))
	  }
		qpcr.result[[i]] <- rt.cal.tg
  write.csv(rt.cal.tg, paste(qpcr.save, "/DATA_", test[i], "_", paste(strsplit(as.character(Sys.time()), " |:|-")[[1]], sep="", collapse = "_"), ".csv",sep = ""))
  write.csv(rt.std.tg, paste(qpcr.save, "/STD_", test[i], "_", paste(strsplit(as.character(Sys.time()), " |:|-")[[1]], sep="", collapse = "_"), ".csv",sep = ""))
  i <- i - 1
  }
return(qpcr.result)
}
