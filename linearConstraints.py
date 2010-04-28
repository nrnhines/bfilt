import numpy

class LinearConstraints(object):
    def __init__(self):
        self.equal = []
        self.greater = []
        self.less = []

    def