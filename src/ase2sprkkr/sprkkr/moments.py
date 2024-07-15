class Moments:

    def __init__(self, qel, nos, smt, omt, hff):
        self.qel = qel
        self.nos = nos
        self.smt = smt
        self.omt = omt
        self.hff = hff

    def as_tuple(self):
        return (self.qel, self. nos, self. smt, self. omt, self. hff)
