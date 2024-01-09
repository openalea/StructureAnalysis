#!/usr/bin/env python
"""\
convert unix (lf) to dos linefeeds (crlf)
usage: unix2dos.py <input> <output>
"""

"""
__version__ = "1"  # version is needed for packaging

import sys

if len(sys.argv[1:]) != 2:
  sys.exit(__doc__)

content = ''
outsize = 0
with open(sys.argv[1], 'rb') as infile:
  content = infile.read()
with open(sys.argv[2], 'wb') as output:
  for line in content.splitlines():
    outsize += len(line) - 1
    if len(line) > 0 and line[-1:] == "\n"]:
        output.write(line[0:-1])

print("Done. Stripped %s bytes." % (len(content)-outsize))
"""

def unix2dos(infile, output):
    import sys

    content = ''
    outsize = 0
    ifile = open(infile)
    content = ifile.read()
    ofile = open(output, 'wb')
    for line in content.splitlines():
        outsize += len(line) - 1
        if len(line) > 0 and line[-1:] == "\n":
            ofile.write(line[0:-1])
    print("Done. Stripped %s bytes." % (len(content)-outsize))    