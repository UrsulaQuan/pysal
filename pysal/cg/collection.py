import numpy as np
import json
import warnings

from ops import ring_bbox, ring_centroid

class Geometry(object):
    """
    Abstract Class for PySAL Geometry Types
    """
    def __init__(self, coordinates, bbox=None, rep_point=None):
        self.coordinates = coordinates
        self.set_bbox(bbox)
        self.set_rep_point(rep_point)

    def __str__(self):
        return self.type + ': ' + str(self.coordinates)

    def __repr__(self):
        return '{0}'.format(object.__repr__(self))

    def get_bbox(self):
        return self.__bbox

    def set_bbox(self, value):
        if value is None:
            value = self.build_bbox()
        self.__bbox = value
    bbox = property(get_bbox, set_bbox)

    def build_bbox(self):
        raise NotImplementedError('Must implement in subclass')

    def set_rep_point(self, value=None):
        if value is None:
            value = self.do_rep_point()
        self.__rep_point = value

    def get_rep_point(self):
        return self.__rep_point

    rep_point = property(get_rep_point, set_rep_point)

    def do_rep_point(self):
        raise NotImplementedError('Must implement in subclass')


class Polygon(Geometry):
    """
    Polygon geometry
    """
    def __init__(self, coordinates):
        """
        Arguments
        ---------
        coordinates: list of lists
                   Each list is a ring. First ring is exterior polygon ring,
                   subsequent rings (if any) are holes
        """
        super(Polygon, self).__init__(coordinates)
        self.__centroid = None

    def __str__(self):
        nrings = len(self.coordinates)
        nholes = nrings - 1
        return '{0}: {1} Rings, {2} Holes'.format(self.type, nrings, nholes)

    def build_bbox(self):
        x0 = np.Infinity
        y0 = np.Infinity
        x1 = -np.Infinity
        y1 = -np.Infinity
        self.bboxes = []
        for r, ring in enumerate(self.coordinates):
            l, b, r, t = ring_bbox(ring)
            self.bboxes.append([l, b, r, t])
            x0 = min(l, x0)
            y0 = min(b, y0)
            x1 = max(r, x1)
            y1 = max(t, y1)
        return [x0, y0, x1, y1]

    @property
    def centroid(self):
        if self.__centroid is None:
            self.__centroid = ring_centroid(self.coordinates[0])
        return self.__centroid

    def do_rep_point(self, method='bbc'):
        """
        determine representative point

        Arguments
        ---------
        method: string
                bbc: center of bounding box
                contained: point will be within exterior ring, not in hole
                centroid: center of mass
        """

        if method == 'bbc':
            # center of bbox
            l, b, t, r = self.bbox
            x = (l + r) / 2.
            y = (t + b) / 2.
            return (x, y)
        elif method.lower() == 'contained':
            # within exterior polygon, not in hole
            raise NotImplementedError
        elif method.lower() == 'centroid':
            # centroid of exterior polygon
            return self.__centroid
        else:
            raise ValueError('method not supported')


class MultiPolygon(Geometry):
    """
    MultiPolygon geometry

    """
    def __init__(self, coordinates):
        """
        Arguments
        ---------
        coordinates: list of Polygon coordinates
        """
        super(MultiPolygon, self).__init__(coordinates)
        self.type = 'MultiPolygon'

    def __str__(self):
        n_polygons = len(self.coordinates)
        return '{0}: {1} Polygon(s)'.format(self.type, n_polygons)

    def build_bbox(self):
        x0 = y0 = np.Infinity
        x1 = y1 = -np.Infinity
        self.bboxes = []
        for p, polygon in enumerate(self.coordinates):
            pbbox = []
            for r, ring in enumerate(polygon):
                l, b, r, t = ring_bbox(ring)
                pbbox.append([l, b, r, t])
                x0 = min(l, x0)
                y0 = min(b, y0)
                x1 = max(r, x1)
                y1 = max(t, y1)
            self.bboxes.append(pbbox)
        return [x0, y0, x1, y1]


class LineString(Geometry):
    def __init__(self, coordinates):
        super(LineString, self).__init__(coordinates)

    def build_bbox(self):
        return ring_bbox(self.coordinates)


class MultiLineString(Geometry):
    def __init__(self, coordinates):
        super(MultiLineString, self).__init__(coordinates)

    def build_bbox(self):
        x0 = y0 = np.Infinity
        x1 = y1 = -np.Infinity
        for line_string in self.coordinates:
            l, b, r, t = LineString(line_string).bbox
            x0 = min(l, x0)
            y0 = min(b, y0)
            x1 = max(r, x1)
            y1 = max(t, y1)
        return [x0, y0, x1, y1]


class Point(Geometry):
    def __init__(self, coordinates):
        super(Point, self).__init__(coordinates)

    def build_bbox(self):
        return self.coordinates


class MultiPoint(Geometry):
    def __init__(self, coordinates):
        super(MultiPoint, self).__init__(coordinates)

    def build_bbox(self):
        # assume all points have same number of coordinates
        m = np.array([point for point in self.coordinates])
        minc = m.min(axis=0).tolist()
        maxc = m.max(axis=0).tolist()
        return minc+maxc

geometry_dispatcher = {}
geometry_dispatcher[u'Point'] = Point
geometry_dispatcher[u'Polygon'] = Polygon
geometry_dispatcher[u'LineString'] = LineString
geometry_dispatcher[u'MultiPoint'] = MultiPoint
geometry_dispatcher[u'MultiPolygon'] = MultiPolygon
geometry_dispatcher[u'MultiLineString'] = MultiLineString


class FeatureCollection(object):
    def __init__(self, file_path, id_variable=None):
        """
        FeatureCollection extension of geojson FC for handling internal data in
        PySAL

        Arguments
        --------
        source: file path or json string
                input data

        id_variable: int
                Name of feature property that is used for spatial joins,
                sub-setting, and sorting of features

        """
        with open(file_path) as source:
            data = json.load(source)
        features = {}
        self.n_features = 0
        for i, feature in enumerate(data['features']):
            ft = feature['geometry']['type']
            fc = feature['geometry']['coordinates']
            f = {}
            f['geometry'] = geometry_dispatcher[ft](fc)
            f['properties'] = feature['properties']
            features[i] = f
            self.n_features += 1
        self.features = features
        self.id_variable = id_variable

    def get_property_types(self):
        propertyTypes = {}
        for key in self.features[0]['properties']:
            propertyTypes[key] = type(self.features[0]['properties'][key])
        return propertyTypes

    def get_properties_as_lists(self, properties):
        a = []
        for i in xrange(self.n_features):
            f = self.features[i]['properties']
            a.append([f[p] for p in properties])
        return a

    def get_properties_as_array(self, properties):
        return np.array(self.get_properties_as_lists(properties))

    def get_geometry_coordinates(self):
        fs = self.features.itervalues()
        return [f['geometry'].coordinates for f in fs]

    @property
    def id_variable(self):
        return self.__id_variable

    @id_variable.setter
    def id_variable(self, property_name):
        if property_name is None:
            self.__id_variable = property_name
        elif self.is_unique(self.get_properties_as_array([property_name])):
            self.__id_variable = property_name
        else:
            msg = "{0} is not uniquely valued. ".format(property_name)
            msg = "{0} Cannot be used as id_variable".format(msg)
            warnings.warn(msg)

    def get_id_variable_order(self):
        return self.get_properties_as_array([self.__id_variable])

    def get_feature_ids(self, args):
        """
        Find the feature ids for features with id_variables matching args

        Arguments
        ---------
        args: list
              values that id_variable matches on

        Returns
        -------
        _  :  list
              feature ids
        """
        idOrder = self.get_id_variable_order().flatten().tolist()
        ids = [idOrder.index(v) for v in args]
        return ids

    def get_features(self, args=None):
        """
        Return features matching on id_variable for arg

        Arguments
        ---------
        args: list
              values that id_variable is to be matched on. If args is None, all
              features are returned

        Returns

        _  : generator
             features matching query

        """

        if args:
            ids = self.get_feature_ids(args)
        else:
            ids = xrange(self.n_features)
        return (self.features[i] for i in ids)

    def get_geometries(self, args=None):
        if args:
            ids = self.get_feature_ids(args)
        else:
            ids = xrange(self.n_features)
        return (self.features[i]['geometry'] for i in ids)

    def is_unique(self, values):
        """
        Check if the values in the array are all unique

        Arguments
        --------
        values: array (1d)
        """
        if np.unique(values).shape[0] == values.shape[0]:
            return True
        else:
            return False

    def save_as_geojson(self, fileName="Untiled.geojson"):
        d = {}
        d["type"] = "FeatureCollection"
        d['features'] = []
        for fi in self.features:
            feature = self.features[fi]
            f = {}
            f["type"] = "Feature"
            f["geometry"] = {}
            f["geometry"]["type"] = feature["geometry"].type
            f["geometry"]["coordinates"] = feature["geometry"].coordinates
            f["properties"] = feature['properties']
            d['features'].append(f)
        with open(fileName, 'w') as out:
            json.dump(d, out)

    def get_rep_points(self):
        """
        Get representative points for all feature geometries

        Returns
        -------
        _  : array (n x k)
             each row is the representative point for a feature geometry
        """
        fit = self.features.iteritmes()
        return np.array([f['geometry'].representative_point for i, f in fit])
