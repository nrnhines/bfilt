loadxyz <- function(fname)
{  x = scan(file=paste(fname,"_x.txt",sep=""),what=0)
    y = scan(file=paste(fname,"_y.txt",sep=""),what=0)
    z = matrix(scan(file=paste(fname,"_z.txt",sep=""),what=0),nrow=length(x),ncol=length(y),byrow=TRUE)
    list(x=x,y=y,z=z)
}
loadxyc <- function(fname)
{  x = scan(file=paste(fname,"_x.txt",sep=""),what=0)
    y = scan(file=paste(fname,"_y.txt",sep=""),what=0)
    z = matrix(scan(file=paste(fname,"_c.txt",sep=""),what=0),nrow=length(x),ncol=length(y),byrow=TRUE)
    list(x=x,y=y,z=z)
}

