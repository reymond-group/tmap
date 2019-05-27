import sys
from tmap import web

# web.host('test.dat', 'mnist.labels')
print('Loading ' + sys.argv[1] + ' ...')

label_type = 'smiles'
if len(sys.argv) > 2:
    label_type = sys.argv[2]

web.host_static(sys.argv[1], label_type)
