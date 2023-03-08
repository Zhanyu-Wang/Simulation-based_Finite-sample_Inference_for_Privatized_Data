
library(graphics)
lower = function(x){
  bottom = (x+1)^2-1
  middle = x-1/4
  top = x^2
  return(bottom*(x<=-1/2)+middle*(x>-1/2 & x<1/2)+top*(x>=1/2))
}

upper = function(x){
  bottom = -x^2
  middle = x+1/4
  top = -(x-1)^2+1
  return(bottom*(x<=-1/2)+middle*(x>-1/2 & x<1/2)+top*(x>=1/2))
}

  
SSlower = function(x){
  I = x<0
  return(((x+1)^2-1)*I+(x^2)*(1-I))
}

SSupper = function(x){
  return(pmax(-x^2,1-(1-x)^2))
}
  
vals = seq(-1,1,by=.01)

SSy = c(SSupper(vals),rev(SSlower(vals)))
x = c(vals,rev(vals))

Cy = c(upper(vals),rev(lower(vals)))

pdf("sensitivityHull.pdf",width=5,height=5)
plot(x,SSy,type="l",xlab=expression(u[1]),ylab=expression(u[2]))
lines(x,Cy)
polygon(x,Cy,col=rgb(.5, .5, 1,0.5))
polygon(x,SSy,col=rgb(0,0,1,1))
dev.off()
