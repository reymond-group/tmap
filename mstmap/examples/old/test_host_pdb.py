import sys
from tmap import web


def label_formatter(label, index, name):
    label = label.lower()
    return 'https://cdn.rcsb.org/images/rutgers/' + label[1:3] + '/' + label + '/' + label + '.pdb1-500.jpg'


def link_formatter(label, index, name):
    return 'https://www.rcsb.org/structure/' + label


# web.host('test.dat', 'mnist.labels')
print('Loading ' + sys.argv[1] + ' ...')

label_type = 'urlimage'
if len(sys.argv) > 2:
    label_type = sys.argv[2]

web.host(sys.argv[1], label_type,
         label_formatter=label_formatter, link_formatter=link_formatter)
