import Gnuplot

def gplotxy(x,y):
    xy = []
    for i in range(len(x)):
        xy.append([x[i],y[i]])
    g = Gnuplot.Gnuplot()
    g.plot(xy)
    raw_input('Please press return to continue...\n')
    g.reset()
