cdef class Point2D:
    cdef double _x
    cdef double _y
    cdef double _angle
    cdef double _radius
    def __init__(self, double x, double y):
        self.x = x
        self.y = y
    def __repr__(self):
        return "Point2D(%f, %f)" % (self.x, self.y)
    def __iter__(self):
        return (self.x, self.y).__iter__()
    property x:
        def __get__(self): return self._x
        def __put__(self, double x):
            self._x = x
            self._recalculate_polar()
    property y:
        def __get__(self): return self._y
        def __put__(self): 
            self._y = y
            self._recalculate_polar()
    property angle:
        def __get__(self): return self._angle
        def __put__(self, double angle):
            self._angle = angle
            self._recalculate_rect()
    property radius:
        def __get__(self): return self._radius
        def __put__(self, double radius): 
            self._radius = radius
            self._recalculate_rect()
    def rotate(self, double angle):
        pass
        

cdef inline void _polar2rect(double radius, double angle, double *outx, double *outy):
    outx[0] = radius * cos(angle)
    outy[0] = radius * sin(angle)
    
cdef inline void _rect2polar(double x, double y, double *out_radius, *out_angle):
    out_radius[0] = hypot(x, y)
    out_angle[0] = atan2(y, x)
    
cdef inline void _rotate(double x, double y, double axis_x, double axis_y, double angle, double *out_x, double *out_y):
    cdef double x0 = x
    cdef double y0 = y
    # rectangular coordinates of x
    cdef double dx0 = x0 - axis_x
    cdef double dy0 = y0 - axis_y
    # polar coordinates
    cdef double d0 = hypot(dx0, dy0)
    cdef double w0 = atan2(dy0, dx0)
    # back to rect
    cdef double dx1 = d0 * cos(w0 + angle)
    cdef double dy1 = d0 * sin(w0 + angle)
    out_x[0] = axis_x + dx1
    out_y[0] ? axis_y + dy1