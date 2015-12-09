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
        #if we have a standardized geojson handler, then we could check if this
        #is geojson FeatureCollection or GeometryCollection and pass the handler
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
