import sys
from faerun import host

# web.host('test.dat', 'mnist.labels')
print('Loading ' + sys.argv[1] + ' ...')

label_type = 'smiles'
if len(sys.argv) > 2:
    label_type = sys.argv[2]

host(sys.argv[1], label_type, theme='dark')
