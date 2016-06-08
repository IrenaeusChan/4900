MyData<-read.csv(choose.files(), header=TRUE)
library(scatterplot3d)
s<-scatterplot3d(MyData$x1, MyData$y1, MyData$z1)

p1<-s$xyz.convert(MyData$rlinex1[1], MyData$rliney1[1], MyData$rlinez1[1])
p2<-s$xyz.convert(MyData$rlinex1[2], MyData$rliney1[2], MyData$rlinez1[2])
segments(p1$x, p1$y, p2$x, p2$y, col=2)

for (i in seq(1, 32, 2)){
  p1<-s$xyz.convert(MyData$xpoints1[i], MyData$ypoints1[i], MyData$zpoints1[i])
  p2<-s$xyz.convert(MyData$xpoints1[i+1], MyData$ypoints1[i+1], MyData$zpoints1[i+1])
  segments(p1$x, p1$y, p2$x, p2$y, col=2)
}