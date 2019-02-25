import os
import sys
import cherrypy
import mstmap

import matplotlib.pyplot as plt

def index_file(path, out_path):
    indices = []
    with open(out_path, 'w+') as f_out:
        with open(path, 'r') as f_in:
            for line in iter(f_in.readline, ''):
                pos = f_in.tell()
                length = len(line)
                start = pos - length
                indices.append(start)
                indices.append(length)
                f_out.write(str(start) + ',' + str(length))

    return indices



class MstmapWeb(object):
    def __init__(self, path, labels_path = ''):
        if not os.path.isfile(path):
            print('File not found: ' + path)
            sys.exit(1)

        self.path = path
        self.labels_path = labels_path

        self.lf = mstmap.LSHForest()
        self.lf.restore(self.path)

        self.n_labels = 0 
        self.label_indices = mstmap.VectorUlong()
        self.has_labels = False
        self.labels_file = None
        self.load_label_indices()
        
        
    def load_label_indices(self):
        # Index labels
        if os.path.isfile(self.labels_path):
            indices = index_file(self.labels_path, self.labels_path + '.index')
            self.label_indices = mstmap.VectorUlong(indices)
            self.labels_file = open(self.labels_path, 'r')
            self.has_labels = True
            self.n_labels = len(indices) / 2

    @cherrypy.expose
    def index(self):
        return open(mstmap.get_asset('index.html'))

    @cherrypy.expose
    @cherrypy.tools.allow(methods=['POST'])
    @cherrypy.tools.json_out()
    @cherrypy.tools.json_in()
    def save(self):
        self.lf.store(self.path)
        # also store labels
        return True

    @cherrypy.expose
    @cherrypy.tools.allow(methods=['POST'])
    @cherrypy.tools.json_out()
    @cherrypy.tools.json_in()
    def get_coords(self):
        coords = mstmap.layout_from_lsh_forest(self.lf)
        s = 1200
        min_coords = min([min(coords[0]), min(coords[1])])
        diff_coords = min([max(coords[0]), max(coords[1])]) - min_coords

        return [
            [round(s * (x - min_coords) / diff_coords, 3) for x in coords[0]], 
            [round(s * (y - min_coords) / diff_coords, 3) for y in coords[1]]
        ]

    @cherrypy.expose
    @cherrypy.tools.allow(methods=['POST'])
    @cherrypy.tools.json_out()
    @cherrypy.tools.json_in()
    def get_nn(self):
        input_json = cherrypy.request.json
        id = input_json["id"]
        k = input_json["k"]
        return self.lf.query_linear_scan_by_id(id, k)

    @cherrypy.expose
    @cherrypy.tools.allow(methods=['POST'])
    @cherrypy.tools.json_out()
    @cherrypy.tools.json_in()
    def get_label(self):
        input_json = cherrypy.request.json
        pos = input_json["id"] * 2
        start = self.label_indices[pos]
        length = self.label_indices[pos + 1]
        self.labels_file.seek(start)
        return self.labels_file.read(length)

    @cherrypy.expose
    @cherrypy.tools.allow(methods=['POST'])
    @cherrypy.tools.json_out()
    @cherrypy.tools.json_in()
    def get_label_colors(self):
        input_json = cherrypy.request.json
        colormap = 'plasma'

        # with open(self.labels_path, 'r') as f:
        #     for line in f:

        
        # colors = np.array([plt.cm.get_cmap(colormap)(x) for x in ]])

        # colors = np.round(colors * 255.0)
        
        



def host(path, labels_path = ''):
    cherrypy.quickstart(MstmapWeb(path, labels_path))