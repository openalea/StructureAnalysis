
class Iterator(object):

    def __init__(self, begin, end):
        self.current = begin
        self.end = end

    def next(self):
        if self.current == self.end:
            raise StopIteration()
        current = self.current.__next__(0)
        return current.__mul__()