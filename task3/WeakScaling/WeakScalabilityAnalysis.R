setwd("C:/Users/Daniel/Dropbox/Master/Parallel programming/Assignment3/EntregaMPI/task3/WeakScaling")

baseline <- read.table("Baseline.txt", header = T, sep = ',')

datak1 <- read.table("Wk1.txt", header = T, sep = ',')
datak4 <- read.table("Wk4.txt", header = T, sep = ',')
datak20 <- read.table("Wk20.txt", header = T, sep = ',')


datak1 <- rbind(baseline, datak1)
datak4 <- rbind(baseline, datak4)
datak20 <- rbind(baseline, datak20)


size <- datak1$N
tk1 <- datak1$elapsed_time
tk4 <- datak4$elapsed_time
tk20 <- datak20$elapsed_time

# plot(size,tk1[1]/(tk1))
# points(size, tk4[1]/(tk4), col = "red")
# points(size, tk20[1]/(tk20), col = "blue")
# 
# (tk1 - tk1[1])/c(1,4,16)


#par(mfrow=c(1,1))


factor = 100
plot(size, (tk1[1]/(tk1) ) * factor, ylim = c(0,1*factor),
     #main = "Strong Scaling Efficiency, n = 4096, iterations = 1000",
     main = "Weak Scaling Efficiency",
     ylab = "",
     xlab = "Number of resources (N)",
     cex=1.25, cex.lab = 1.5)
title(ylab=expression('T'[1] ~ "/" ~ 'T'['N']~" x 100 (%)"),
      line=2, cex.lab=1.5)
points(size, tk4[1]/(tk4) * factor, col = "red", cex=1.25)
points(size, tk20[1]/(tk20) * factor, col = "blue", cex=1.25)

legend("bottom",
       c("k=1", "k=4", "k=20"), col = c("black", "red", "blue"),
       bty = 'n',
       pch = 16,
       cex = 1)


####SPEEDUP
# sp_theory <- function(N, s=0,p=1){
#   return(1.0/(p/N+s));
# }
# 
# 
# plot(sp_theory, from = 0, to = 8, ylab = "Speedup", 
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
# points(size, tk4[1]/(tk4), col = "red")
# points(size, tk20[1]/(tk20), col = "blue")
# points(c(0,8), c(1,1), type = 'l')
# 
# legend("topleft", 
#        c("P = 0.75, S = 0.25","P = 1, S = 0"), 
#        col = c("blue","green"),lwd=2, bty ='n', 
#        cex = 0.8, 
#        seg.len = 0.5)