setwd("C:/Users/Daniel/Dropbox/Master/Parallel programming/Assignment3/EntregaMPI/task3/StrongScaling")

baseline <- read.table("Baseline.txt", header = T, sep = ',')

datak1 <- read.table("Sk1.txt", header = T, sep = ',')
datak4 <- read.table("Sk4.txt", header = T, sep = ',')
datak20 <- read.table("Sk20.txt", header = T, sep = ',')

datak1 <- rbind(baseline, datak1)
datak4 <- rbind(baseline, datak4)
datak20 <- rbind(baseline, datak20)

size <- datak1$N
tk1 <- datak1$elapsed_time
tk4 <- datak4$elapsed_time
tk20 <- datak20$elapsed_time

par(mfrow=c(1,3))
#compile all this file and then the WeakScalabiltyAnalysis.R
#to obtain a figure coposed by 3 subfigures


factor = 100
plot(size, (tk1[1]/(tk1 * size) ) * factor, ylim = c(0,1*factor),
     #main = "Strong Scaling Efficiency, n = 4096, iterations = 1000",
     main = "Strong Scaling Efficiency",
     ylab = "",
     xlab = "Number of resources (N)",
     cex=1.25 , cex.lab = 1.5 )
title(ylab=expression('T'[1] ~ "/" ~ '(N · T'['N']~") x 100 (%)"),
      line=2, cex.lab=1.5)

points(size, tk4[1]/(tk4 *size) * factor, col = "red", cex=1.25)
points(size, tk20[1]/(tk20*size) * factor, col = "blue", cex=1.25)

legend("bottom",
       c("k=1", "k=4", "k=20"), col = c("black", "red", "blue"),
       bty = 'n',
       pch = 16,
       cex = 1)


####SPEEDUP
sp_theory <- function(N, s=0,p=1){
  return(1.0/(p/N+s));
}


plot(sp_theory, from = 1, to = 8, ylab = "", 
     xlab = "Number of resources (N)",
     #main = "Speedup, n = 4096, iterations = 1000",
     main = "Speedup",
     col = "green", lwd = 2.5, type = 'l',
     cex.lab = 1.5,
     cex = 1.25)

title(ylab="Speedup",
      line=2, cex.lab=1.5)
#curve(sp_theory(x,0,1), from = 1, to = 8, add = T, lw = 1.5)
curve(sp_theory(x,0.05,0.95), from = 1, to = 8, add = T, 
      lwd = 1, col = "blue")
curve(sp_theory(x,0.03,0.97), from = 1, to = 8, add = T, 
      lwd = 1, col = "black")
curve(sp_theory(x,0.08,0.92), from = 1, to = 8, add = T, 
      lwd = 1, col = "red")

# points(size,as.numeric(tk20)[1]/as.numeric(tk20),pch=16, cex=1.25, col = "blue" )
# points(size,as.numeric(tk4)[1]/as.numeric(tk20),pch=16, cex=1.25, col = "red" )
# points(size,as.numeric(tk1)[1]/as.numeric(tk20),pch=16, cex=1.25, col = 1 )
points(size, tk1[1]/(tk1), col = "black")
points(size, tk4[1]/(tk4), col = "red")
points(size, tk20[1]/(tk20), col = "blue")

legend("topleft", 
       c("P = 0.95, S = 0.05", 
         "P = 0.97, S = 0.03", 
         "P = 0.92, S = 0.08",
         "P = 1, S = 0"), 
       col = c("blue", "black", "red","green"),lwd=2, bty ='n', 
       cex = 0.8, 
       seg.len = 0.5)
