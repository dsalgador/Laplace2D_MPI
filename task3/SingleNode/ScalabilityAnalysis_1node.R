setwd("C:/Users/Daniel/Dropbox/Master/Parallel programming/Assignment3/EntregaMPI/task3/SingleNode/StrongScaling")
library(dplyr)
library(data.table)
datastrong <- read.table("StrongScalingSingleNode_k1R.csv", header = T, sep = ',')

datastrong <- as.data.table(datastrong)
datastrong[, speedup:= elapsed_time[1:4] / elapsed_time]
datastrong[, strongeff:= 100 * elapsed_time[1:4] / (N* elapsed_time)]

require(ggplot2)
#par(mfrow=c(1,3))


dt_agg <- datastrong[,.(speedup),by=list(n,N)]
p1 <- ggplot(dt_agg, aes(x = N, y = speedup, col = n)) + geom_point()
p1 <- p1 + labs(x = "Number of threads (N)")
p1 <- p1 + labs(y = "Speedup")
p1 <- p1 + labs(title = "Speedup")


dt_agg2 <- datastrong[,.(strongeff),by=list(n,N)]
p2 <- ggplot(dt_agg2, aes(x = N, y = strongeff, col = n)) + geom_point()
p2 <- p2 + labs(x = "Number of threads (N)")
p2 <- p2 + labs(y = expression('T'[1] ~ "/" ~ '(N · T'['N']~") x 100 (%)") )
p2 <- p2 + labs(title = "Strong Scaling Efficiency")

setwd("C:/Users/Daniel/Dropbox/Master/Parallel programming/Assignment3/EntregaMPI/task3/SingleNode/WeakScaling")
dataweak <- read.table("WeakScalingSingleNode_k1R.csv", header = T, sep = ',')

dataweak <- as.data.table(dataweak)
dataweak[, speedup:= elapsed_time[1:4] / elapsed_time]
dataweak[, weakeff:= 100 * elapsed_time[1] / (elapsed_time)]

dt_agg2 <- dataweak[,.(weakeff),by=list(n,N)]
p3 <- ggplot(dt_agg2, aes(x = N, y = weakeff, col = n)) + geom_point()
p3 <- p3 + labs(x = "Number of threads (N)")
p3 <- p3 + labs(y = expression('T'[1] ~ "/" ~ 'T'['N']~" x 100 (%)") )
p3 <- p3 + labs(title = "Weak Scaling Efficiency")


setwd("C:/Users/Daniel/Dropbox/Master/Parallel programming/Assignment3/EntregaMPI/task3/SingleNode")

png("speedup.png", width = 347, height = 274)
print(p1)
dev.off()

png("strongeff.png", width = 347, height = 274)
print(p2)
dev.off()


png("weakeff.png", width = 347, height = 274)
print(p3)
dev.off()

# 
# 
# size <- datak1$N
# tk1 <- datak1$elapsed_time
# # tk4 <- datak4$elapsed_time
# # tk20 <- datak20$elapsed_time
# 
# par(mfrow=c(1,3))
# #compile all this file and then the WeakScalabiltyAnalysis.R
# #to obtain a figure coposed by 3 subfigures
# 
# 
# factor = 100
# plot(size, (tk1[1]/(tk1 * size) ) * factor, ylim = c(0,1*factor),
#      #main = "Strong Scaling Efficiency, n = 4096, iterations = 1000",
#      main = "Strong Scaling Efficiency",
#      ylab = expression('T'[1] ~ "/" ~ '(N · T'['N']~") x 100 (%)"),
#      xlab = "Number of resources (N)",
#      cex=1.25)
# # points(size, tk4[1]/(tk4 *size) * factor, col = "red", cex=1.25)
# # points(size, tk20[1]/(tk20*size) * factor, col = "blue", cex=1.25)
# 
# legend("bottom",
#        c("k=1", "k=4", "k=20"), col = c("black", "red", "blue"),
#        bty = 'n',
#        pch = 16,
#        cex = 1)
# 
# 
# ####SPEEDUP
# sp_theory <- function(N, s=0,p=1){
#   return(1.0/(p/N+s));
# }
# 
# 
# plot(sp_theory, from = 1, to = 8, ylab = "Speedup", 
#      xlab = "Number of resources (N)",
#      #main = "Speedup, n = 4096, iterations = 1000",
#      main = "Speedup",
#      col = "green", lwd = 2.5, type = 'l')
# #curve(sp_theory(x,0,1), from = 1, to = 8, add = T, lw = 1.5)
# curve(sp_theory(x,0.25,0.75), from = 1, to = 8, add = T, lwd = 2.5, col = "blue")
# # points(size,as.numeric(tk20)[1]/as.numeric(tk20),pch=16, cex=1.25, col = "blue" )
# # points(size,as.numeric(tk4)[1]/as.numeric(tk20),pch=16, cex=1.25, col = "red" )
# # points(size,as.numeric(tk1)[1]/as.numeric(tk20),pch=16, cex=1.25, col = 1 )
# points(size, tk1[1]/(tk1), col = "black")
# # points(size, tk4[1]/(tk4), col = "red")
# # points(size, tk20[1]/(tk20), col = "blue")
# 
# legend("topleft", 
#        c("P = 0.75, S = 0.25","P = 1, S = 0"), 
#        col = c("blue","green"),lwd=2, bty ='n', 
#        cex = 0.8, 
#        seg.len = 0.5)
