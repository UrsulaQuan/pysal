
def asPolygonCollection(collection, **kwargs):
    """
    Construct a PolygonCollection from an arbitrary iterable of shapes.
    """
    ids = kwargs.pop('ids', None)
    if ids is None:
        ids = ()
    if isinstance(collection, PolygonCollection):
        out = collection
    elif isinstance(collection, dict):
        # if we have a standardized geojson handler, then we could check if this
        # is geojson FeatureCollection or GeometryCollection and pass the handler
        collection = {k:asShape(v) for k,v in diter(collection)}
        out = PolygonCollection(collection, **kwargs)
    else:
        try:
            collection = {i:asShape(v) for i,v in enumerate(collection)}
        except TypeError:
            raise TypeError("collection is not iterable")
        out = PolygonCollection(collection, **kwargs)
    return out

def asShape(obj):
    """
    Returns a pysal shape object from obj.
    obj must support the __geo_interface__.
    """
    if type(obj) in _geoJSON_type_to_Pysal_type.values():
        return obj #already pysal object
    if hasattr(obj, '__geo_interface__'):
        geo = obj.__geo_interface__
    else:
        geo = obj
    if hasattr(geo, 'type'):
        raise TypeError('%r does not appear to be a shape object' % (obj))
    geo_type = geo['type'].lower()
    #if geo_type.startswith('multi'):
    #    raise NotImplementedError, "%s are not supported at this time."%geo_type
    if geo_type in _geoJSON_type_to_Pysal_type:
        return _geoJSON_type_to_Pysal_type[geo_type].__from_geo_interface__(geo)
    else:
        raise NotImplementedError(
            "%s is not supported at this time." % geo_type)

def ring_bbox(ring):
    xs = [c[0] for c in ring]
    ys = [c[1] for c in ring]
    return min(xs), min(ys), max(xs), max(ys)


def ring_centroid(ring):
    x = [v[0] for v in ring]
    y = [v[1] for v in ring]
    n = len(x)
    a = 0.0  # area
    for i in xrange(n-1):
        a += (x[i] * y[i+1] - x[i+1] * y[i])
    a /= 2.0
    cx = cy = 0.0
    for i in xrange(n-1):
        f = (x[i] * y[i + 1] - x[i + 1] * y[i])
        cx += (x[i] + x[i + 1]) * f
        cy += (y[i] + y[i + 1]) * f
    cx = 1.0 / (6 * a) * cx
    cy = 1.0 / (6 * a) * cy
    return (cx, cy)


def bbox_overlap(A, B):
    """
    Test if two bounding boxes overlap

    Arguments
    ---------
    A, B: list of [left, bottom, right, top]
    """
    l = 0
    b = 1
    r = 2
    t = 3
    if A[b] > B[t] or A[t] < B[b] or A[r] < B[l] or A[l] > B[r]:
        return False
    else:
        return True

def _queen(ring1, ring2):
    common_vertices = [v for v in ring1 if v in ring2]
    if common_vertices:
        return True
    return False


def _rook(ring1, ring2):
    n1 = len(ring1)
    n2 = len(ring2)
    edges1 = [sorted(ring1[l:l+2]) for l in range(n1-1)]
    edges2 = [sorted(ring2[l:l+2]) for l in range(n2-1)]
    common_edges = [e for e in edges1 if e in edges2]
    if common_edges:
        return True
    return False

def contiguous(g1, g2, criterion='QUEEN'):

    if type(g1) != type(g2):
        raise ValueError('Geometries must be of same type.')
    elif not bbox_overlap(g1.bbox, g2.bbox):
        # bounding box containing all g1 rings not overlapping that of g2
        return False
    else:
        criterion = criterion.lower()
        # check if the bounding boxes for any rings in g1 overlap with any rings in g2
        overlap = False
        if g1.type == 'Polygon':
            # find all pairs where a ring from g1 and a ring from g2 have overlapping bbs
            pairs = []
            for i, bbox1 in enumerate(g1.bboxes):
                for j, bbox2 in enumerate(g2.bboxes):
                    if bbox_overlap(bbox1, bbox2):
                        overlap=True
                        pairs.append((i,j))
            if overlap:
                # search overlapping pairs and the first (if any) neighbor pair
                for pair in pairs:
                    i, j = pair
                    pi = g1.coordinates[i]
                    pj = g2.coordinates[j]
                    if _queen(pi, pj):
                        # queen is a necessary condition for rook
                        if criterion == 'queen':
                            return True
                        else:
                            return _rook(pi, pj)
                return False
            else:
                return False

        elif g1.type == 'MultiPolygon':
            pairs = []
            overlap = False
            for i, pi in enumerate(g1.bboxes):
                for ii, bi in enumerate(pi):
                    for j, pj in enumerate(g2.bboxes):
                        for jj, bj in enumerate(pj):
                            if bbox_overlap(bi, bj):
                                overlap = True
                                pair = (i,ii,j,jj)
                                pairs.append(pair)

            if overlap:
                for pair in pairs:
                    i, ii, j, jj = pair
                    pi = g1.coordinates[i][ii]
                    pj = g2.coordinates[j][jj]
                    if _queen(pi, pj):
                        if criterion == 'queen':
                            return  True
                        else:
                            return _rook(pi, pj)
                return False
            else:
                return False
        else:
            return False


