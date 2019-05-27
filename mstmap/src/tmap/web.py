import os
import pickle
import sys

import cherrypy
import matplotlib.pyplot as plt
import numpy as np
from timeit import default_timer as timer
import tmap
import ujson


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

# The default cherrypy json encoder seems to be extremely slow...


def json_handler(*args, **kwargs):
    value = cherrypy.serving.request._json_inner_handler(*args, **kwargs)
    return ujson.dumps(value).encode('utf8')


class TmapWebStatic():
    def __init__(self, path, label_type='smiles', label_formatter=None, link_formatter=None):
        if not os.path.isfile(path):
            print('File not found: ' + path)
            sys.exit(1)

        if label_formatter is None:
            label_formatter = lambda label, index, name: label

        if link_formatter is None:
            link_formatter = lambda label, index, name: ''

        self.label_type = label_type
        self.data = pickle.load(open(path, 'rb'))
        self.ids = {}
        self.link_formatter = link_formatter
        self.label_formatter = label_formatter

        for name in self.data:
            if self.data[name]['type'] != 'scatter':
                continue

            # This is currently implemented for two values:
            # e.g. smiles and id
            # seperated by __ in the label field
            self.ids[name] = []

            if '__' in self.data[name]['labels'][0]:
                for _, label in enumerate(self.data[name]['labels']):
                    vals = label.split('__')
                    self.ids[name].append(vals[1].lower())
            else:
                for _, label in enumerate(self.data[name]['labels']):
                    self.ids[name].append(label.lower())

    @cherrypy.expose
    def index(self, **params):
        return open(tmap.get_asset('index_static.html'))

    @cherrypy.expose
    @cherrypy.tools.allow(methods=['POST'])
    @cherrypy.tools.json_out(handler=json_handler)
    @cherrypy.tools.json_in()
    def get_meta(self):
        meta = {}
        meta['label_type'] = self.label_type

        for name in self.data:
            data_type = self.data[name]['type']
            if data_type not in meta:
                meta[data_type] = {}
            meta[data_type][name] = self.data[name]['meta']
        return meta

    @cherrypy.expose
    @cherrypy.tools.allow(methods=['POST'])
    # @cherrypy.tools.json_out(handler=json_handler)
    @cherrypy.tools.json_in()
    def get_values(self):
        input_json = cherrypy.request.json
        name = input_json['name']
        coord = input_json['coord']
        dtype = input_json['dtype']

        if dtype == 'float32':
            dtype = np.float32
        elif dtype == 'uint8':
            dtype = np.uint8

        return bytes(np.array(self.data[name][coord].tolist(), dtype=dtype))

    @cherrypy.expose
    @cherrypy.tools.allow(methods=['POST'])
    @cherrypy.tools.json_out(handler=json_handler)
    @cherrypy.tools.json_in()
    def get_label(self):
        input_json = cherrypy.request.json
        index = input_json['id']
        name = input_json['name']
        label = self.data[name]['labels'][index]
        return {
            'label': self.label_formatter(label, index, name),
            'link': self.link_formatter(label, index, name)
        }

    @cherrypy.expose
    @cherrypy.tools.allow(methods=['POST'])
    @cherrypy.tools.json_out(handler=json_handler)
    @cherrypy.tools.json_in()
    def get_index(self):
        input_json = cherrypy.request.json
        labels = input_json['label']
        name = input_json['name']

        labels = str(labels).upper().strip()

        results = []
        for label in labels.split(','):
            label = label.strip().lower()
            try:
                results.append([label, self.ids[name].index(label)])
            except ValueError:
                results.append([label, -1])

        return results

    @cherrypy.expose
    def ws(self):
        cherrypy.log('Handler created: %s' % repr(cherrypy.request.ws_handler))


def host_static(path, label_type='smiles', label_formatter=None, link_formatter=None):
    cherrypy.quickstart(TmapWebStatic(path, label_type, label_formatter=label_formatter, link_formatter=link_formatter))
