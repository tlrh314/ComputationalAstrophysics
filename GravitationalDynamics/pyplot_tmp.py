# import matplotlib
# matplotlib.use('TkAgg')
from matplotlib import pyplot
#print matplotlib.pyplot.get_backend()

import numpy

x = numpy.linspace(-numpy.pi, numpy.pi, 201)
pyplot.plot(x, numpy.sin(x))
pyplot.show()



print "hello, world"

# matplotlib = reload(matplotlib)
# print matplotlib.pyplot.get_backend()
# matplotlib.use('cairo.png')
# print matplotlib.pyplot.get_backend()
