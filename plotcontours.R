loadxyz <- function(fname)
{  x = scan(file=paste(fname,"_x.txt",sep=""),what=0)
    y = scan(file=paste(fname,"_y.txt",sep=""),what=0)
    z = matrix(scan(file=paste(fname,"_z.txt",sep=""),what=0),length(x),length(y),byrow=TRUE)
    list(x=x,y=y,z=z)
}