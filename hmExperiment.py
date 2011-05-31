import hmm

class HME(object):
    def __init__(self,protocol):
        self.protocol = protocol

    def __len__(self,key):
        return len(self.protocol)

    def __getitem__(self,key):
        return self.protocol[key]

    def __setitem__(self,key,value):
        #make sure an hmm
        self.protocol[key] = value

    def __delitem__(self,key):
        del self.protocol[key]

    def __iter__(self):
        return self.protocol

    def __add__(self,other):
        #other hmm or hme
        try:
            return self.protocol + other.protocol
        except:
            return self.protocol + [other]

    def __iadd__(self,other)
        # other hmm or hme
        try:
            self.protocol += other.protocol
        except:
            self.protocol += [other]