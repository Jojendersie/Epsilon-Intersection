import gdb.printing
import gdb

class EiMatrixPrinter:
    "Prints ei::Matrix<> type"

    def __init__(self, val):
        self.val = val

    # Does not work: https://stackoverflow.com/questions/26472066/gdb-pretty-printing-returning-string-from-a-childrens-iterator-but-displaye/29752860
    def children(self):
        rows = self.val.type.template_argument(1)
        cols = self.val.type.template_argument(2)
        data = self.val["m_data"]
        for y in range(rows):
            line = ''
            for x in range(cols):
                idx = x+y*cols
                yield (str(idx), data[idx])

    def to_string(self):
        return 'Matrix'
        #rows = self.val.type.template_argument(1)
        #cols = self.val.type.template_argument(2)
        #ps = '['
        #data = self.val["m_data"]
        #if rows == 1 or cols == 1:
        #    for i in range(rows * cols):
        #        ps += str(data[i])
        #        if (i+1) != rows*cols: ps += ', '
        #    if rows > 1: ps += "]'"
        #    else: ps += ']'
        #else:
        #    for y in range(rows):
        #        for x in range(cols):
        #            ps += "%12.6g" % self.val["m_data"][x+y*cols] # Multiversion line, looks good on hover, but fails in oneline views
        #            if x != cols-1: ps += ' '
        #        if y != rows-1: ps +=']\n['
        #        else: ps += ']'
        #return ps


def build_pretty_printer():
    pp = gdb.printing.RegexpCollectionPrettyPrinter(
        "Epsilon Intersection")
    pp.add_printer('ei::Matrix', '^ei::Matrix.*$', EiMatrixPrinter)
    return pp