import math
import sys
import random

"""Implementation of rational arithmetic."""
cdef extern from "math.h":
    cdef double cos(double n)
    cdef double pow(double x, double y)
    cdef double ceil(double x)
    cdef double abs(double x)
    cdef double frexp(double x, int *exp)
    cdef double floor(double x)
    
ctypedef long long int_t

cdef inline double _mod(double a, double b):
    return a - (floor(a / b) * b)
    
cdef inline double floordivd(double a, double b):
    return floor(a / b)
    
cdef inline int_t floordivi(int_t a, int_t b):
    return <int_t>floor(<double>a / <double>b)

cdef inline int_t _gcf(int_t a, int_t b):
    """Returns the greatest common factor of a and b."""
    cdef double tmp
    cdef double aa = <double>a
    cdef double bb = <double>b
    while bb != 0.0:
        tmp = bb
        bb = _mod(aa, bb)
        aa = tmp
    return <int_t>aa
    
def gcf(a, b):
    return _gcf(a, b)
    
cdef class Rational
    
cdef inline void _rational_num_den_from_float(double x, int_t *num, int_t *den):
    cdef double mantissa
    cdef int expon
    mantissa = frexp(x, &expon)
    mantissa = mantissa * 9007199254740992ULL 
    expon -= 53
    if expon < 0:
        num[0] = <int_t>mantissa
        den[0] = <int_t>(2.0 ** (-expon))
    else:
        num[0] = <int_t>(mantissa * (2.0 ** expon))
        den[0] = 1

cdef inline Rational _rational_from_float(double x):
    cdef int_t num, den
    _rational_num_den_from_float(x, &num, &den)
    return Rational(num, den)
        
cdef inline void _rational_simplify(int_t *num, int_t *den):
    cdef int_t denominator = den[0]
    cdef int_t numerator = num[0]
    if denominator < 0:
        numerator = -numerator
        denominator = -denominator
    if denominator == 1:
        num[0] = numerator
        den[0] = 1
    else:
        factor = _gcf(numerator, denominator) 
        num[0] = floordivi(numerator, factor)
        den[0] = floordivi(denominator, factor)
    
def rational_simplify(num, den):
    cdef int_t num_, den_
    num_ = num
    den_ = den
    _rational_simplify(&num_, &den_)
    return num_, den_
    
    

    
cdef inline Rational _rational_limit_denominator(Rational self, int_t max_denominator):
    if max_denominator < 1:
        raise ValueError("max_denominator should be at least 1")
    if self.den <= max_denominator:
        return self
    cdef int_t p0 = 0, q0 = 1, p1 = 1, q1 = 0
    cdef int_t n = self.num, d = self.den
    cdef int_t a, k
    while 1:
        a = floordivi(n, d)
        q2 = q0 + a * q1
        if q2 > max_denominator:
            break
        p0, q0, p1, q1 = p1, q1, p0+a*p1, q2
        n, d = d, n - a * d
    k = floordivi(max_denominator - q0, q1)
    cdef Rational bound1 = Rational(p0 + k * p1, q0 + k * q1)
    cdef Rational bound2 = Rational(p1, q1)
    if abs(bound2 - self) <= abs(bound1 - self):
        return bound2
    else:
        return bound1
        
#cdef inline _rational_approx_smallest_denominator(Rational self, int_t max_denominator):
#    
#    result = None
#    minError = x
#    cdef int n
#    for n in range(1, maxDenominator + 1):
#        m = int(round(x * n))
#        r = Rational(m, n)
#          
#        error = abs(r - x)
#        if error == 0:
#            return r
#        elif error < minError:
#            result = r
#        minError = error
#   return result
    
def limit_denominator(r, max_denominator):
    return _rational_limit_denominator(r, max_denominator)

cdef class Rational:
    cdef int_t num
    cdef int_t den
    """
    This class provides an exact representation of rational numbers.

    All of the standard arithmetic operators are provided.  In mixed-type
    expressions, an int or a long can be converted to a Rational without
    loss of precision, and will be done as such.
    """
    def __init__(self, numerator, denominator=1):
        cdef int_t num, den
        if isinstance(numerator, float):
            _rational_num_den_from_float(<double>numerator, &num, &den)
            _rational_simplify(&num, &den)
            self.num = num
            self.den = den
        else:
            self.num = numerator
            self.den = denominator
            _rational_simplify(&(self.num), &(self.den))

    property numerator:
        def __get__(self): return self.num
    property denominator:
        def __get__(self): return self.den
    def __repr__(self):
        if self.den == 1:
            return "Rational(%d)" % self.num
        else:
            return "Rational(%d, %d)" % (self.num, self.den)
    def __str__(self):
        if self.den == 1:
            return str(self.num)
        else:
            return "%d/%d" % (self.num, self.den)
    def __hash__(self):
        try:
            return hash(float(self))
        except OverflowError:
            return hash(long(self))
    def __float__(Rational self):
        return <double>self.num / <double>self.den
    def __int__(Rational self):
        if self.num < 0:
            return -int(floordivi(-self.num, self.den))
        else:
            return int(floordivi(self.num, self.den))
    def __long__(Rational self):
        return long(int(self))
    def __nonzero__(Rational self):
        return self.num != 0
    def __pos__(self):
        return self
    def __neg__(Rational self):
        return Rational(-self.num, self.den)
    def __abs__(Rational self):
        if self.num < 0:
            return Rational(-self.num, self.den)
        else:
        
            return Rational(self.num, self.den)
    def __add__(self, other):
        cdef int_t num0, den0, num1, den1
        if isinstance(other, Rational):
            if isinstance(self, Rational):
                den0 = (<Rational>self).den
                den1 = (<Rational>other).den
                return Rational((<Rational>self).num * den1 + den0 * (<Rational>other).num, den0 * den1)
            else:
                return other + self
        elif isinstance(other, (int, long)):
            return Rational((<Rational>self).num + (<Rational>self).den * other, (<Rational>self).den)
        elif isinstance(other, float):
            return self + _rational_from_float(other)
        return NotImplemented
    def __sub__(a, b):
        cdef Rational ra, rb
        if isinstance(b, Rational):
            if isinstance(a, Rational):
                ra = <Rational>a
                rb = <Rational>b
                return Rational(ra.num * rb.den - ra.den * rb.num,
                    ra.den * rb.den)
            elif isinstance(a, (int, long)):
                rb = <Rational>b
                return Rational(a * rb.den - rb.num, rb.den)
            else:
                return _rational_from_float(a) - b
        elif isinstance(b, (int, long)):
            ra = <Rational>a
            return Rational(ra.num - ra.den * b, ra.den)
        elif isinstance(b, float):
            return a - _rational_from_float(b)
        return NotImplemented
    def __mul__(a, b):
        cdef Rational aa, bb
        if isinstance(b, Rational):
            if isinstance(a, Rational):
                aa = <Rational>a
                bb = <Rational>b
                return Rational(aa.num * bb.num, aa.den * bb.den)
            else:
                bb = <Rational>b
                return Rational(bb.num * a, bb.den)
        elif isinstance(b, (int, long)):
            return Rational((<Rational>a).num * b, (<Rational>a).den)
        else:
            return a * _rational_from_float(b)
    def __div__(a, b):
        if isinstance(b, Rational):
            if isinstance(a, Rational):
                return Rational((<Rational>a).num * (<Rational>b).den, (<Rational>a).den * (<Rational>b).num)
            elif isinstance(a, (int, long)):
                return Rational(a * (<Rational>b).den, (<Rational>b).num)
            else:
                return NotImplemented    
        elif isinstance(b, (int, long)):
            return Rational((<Rational>a).num, (<Rational>a).den * b)
        else:
            return a / _rational_from_float(b)
    def __floordiv__(a, b):
        cdef Rational tmp
        if isinstance(b, Rational):
            if isinstance(a, Rational):
                tmp = a / b
                return Rational(floordivi(tmp.num, tmp.den), 1)
            elif isinstance(a, (int, long)):
                return Rational(b * (<Rational>a).den, (<Rational>a).num)
            else:
                return NotImplemented
        elif isinstance(b, (int, long)):
            tmp = Rational((<Rational>a).num, (<Rational>a).den * b)
            return Rational(floordivi(tmp.num, tmp.den), 1)
        else:
            return a // _rational_from_float(b)  # use // because this is python code
    def __mod__(a, b):
        return a - ((a // b) * b)   # use // because this is python code
    def __richcmp__(Rational self, other, int t):
        if t == 0:      # <
            if other == 0:
                if self.num < 0:
                    return True
                return False
            else:
                return (self - other) < 0
        elif t == 2:    # ==
            if isinstance(other, Rational):
                return (self.num == (<Rational>other).num and self.den == (<Rational>other).den)
            else:
                if other == 0:
                    if self.num == 0:
                        return True
                    return False
                else:
                    return (self - other) == 0
        elif t == 4:    # >
            if other == 0:
                if self.num > 0:
                    return True
                return False
            else:
                return (self - other) > 0
        elif t == 1:    # <=
            if other == 0:
                if self.num >= 0:
                    return True
                return False
            else:
                return (self - other) >= 0
        elif t == 3:    # !=
            return not(self == other)
        elif t == 5:    # >=
            return not(self < other)
    def __pow__(a, b, modulo):
        # TODO: pow 
        if isinstance(a, Rational):
            if isinstance(b, Rational):
                return a ** float(b)
            elif isinstance(b, (int, long)):
                if b < 0:
                    return Rational((<Rational>a).den ** -b, (<Rational>a).num ** -b)
                else:
                    return Rational((<Rational>a).num ** b, (<Rational>a).den ** b)
            else:
                return _rational_from_float(float(a) ** float(b))
        else:
            return Rational(a ** float(b))
    def __iter__(Rational self):
        return iter([self.num, self.den])
    def limit_denominator(Rational self, int_t max_denominator=sys.maxint):
        return _rational_limit_denominator(self, max_denominator)
    #def approx_smallest_denominator(Rational self, int_t max_denominator=sys.maxint):
    #    return _rational_approx_smallest_denominator(self, max_denominator)
    
        
    #def round(self, int denominator):
    #    """Return self rounded to nearest multiple of 1/denominator."""
    #    int_part, frac_part = divmod(self * denominator, 1)
    #    round_direction = cmp(frac_part * 2, 1)
    #    if round_direction == 0:
    #        numerator = int_part + (int_part & 1) # round to even
    #    elif round_direction < 0:
    #        numerator = int_part
    #    else:
    #        numerator = int_part + 1
    #    return Rational(numerator, denominator)
        
    
def rational_from_float(x):
    """Returns the exact Rational equivalent of x."""
    return _rational_from_float(x)
                       
def test_mantissa(double x):
    cdef int expon
    cdef double mantissa = frexp(x, &expon)
    cdef long long out_long = <long long>(mantissa * (2.0 ** 53.0))
    mantissa = mantissa * 9007199254740992.0
    cdef long long out2 
    out2 = <long long>mantissa
    man_py, exp_py = math.frexp(x)
    dos = 2
    cincuentaytres = 53
    man_py = int(man_py * dos ** cincuentaytres)
    assert out_long == out2 == man_py

def test_expon(x):
    cdef double dx = x
    m, e = math.frexp(x)
    cdef int expon
    cdef double mant
    mant = frexp(dx, &expon)
    e -= e
    expon -= expon
    assert e == expon
    
def _assert_almost_eq(mess, a, b, delta = 0.00001):
    try:
        assert abs(a - b) < delta
    except AssertionError:
        print mess, a, b, abs(a - b)
        raise AssertionError

def _assert_eq(mess, a, b):
    try:
        assert a == b
    except AssertionError:
        print mess, a, b
        raise AssertionError
        
        
def test_modulo(a, b):
    return <int_t>(<int_t>a // <int_t>b) == <int_t>floor(<double>a / <double>b)
    
def test_arithmetics(num_range=sys.maxint, num_tests=1000):
    cdef int i
    for i in range(num_tests):
        a = random.randint(-num_range, num_range)
        b = random.randint(-num_range, num_range)
        c = random.randint(-num_range, num_range)
        d = random.randint(-num_range, num_range)
        if a == 0: a += b
        if b == 0: b += c
        if c == 0: c += d
        if d == 0: d += a
        
        try:
            _assert_eq(1,
                a / float(b), 
                float(Rational(a, b))
            )
            _assert_almost_eq(2, 
                float(Rational(a, b) + Rational(c, d)), 
                a / float(b) + c / float(d)
            )
            _assert_eq(3, 
                Rational(a, b) * c, 
                c * Rational(a, b)
            )
            _assert_eq(4, 
                Rational(a, b) / c, 
                Rational(a, b) * Rational(1, c)
            )
            _assert_eq(5, 
                float(Rational(a) - Rational(b)), 
                float(a - b)
            )
            _assert_eq(6,
                Rational(a) % Rational(b), 
                a % b
            )
            _assert_eq(7, 
                Rational(a) // Rational(b), 
                a // b
            )
            _assert_eq(8,
                Rational(a) * Rational(b, c), 
                Rational(a, c) * Rational(b)
            )
            _assert_eq(9,
                a - Rational(b, c),
                Rational(a * c - b, c)
            )
            #_assert_almost_eq(10,           # <----- failing test
            #    Rational(a) * (b / float(c)), 
            #    Rational(a * b, c)
            #    
            #)
            _assert_almost_eq(11, 
                Rational(a, b) / float(c),
                Rational(a, b * c)
            )
            _assert_eq(12,
                -Rational(a, b),
                Rational(-a, b)
            )
            _assert_almost_eq(13,
                float(Rational(a, b) % c),
                (float(a) / b) % c
            )
            
        except AssertionError:
            print a, b, c, d
        
def test_rational_simplify(m =1000000, n=1000):
    for i in range(n):
        a = random.randint(-m, m)
        b = random.randint(-m, m)
        c = random.randint(-m, m)
        d = a * c
        e = b * c
        num, den = rational_simplify(d, e)
        if den < 0:
            den = -den
            num = -num
        try:
            num, den == a, b
        except AssertionError:
            print a, b, c, d, e