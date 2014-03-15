#cython: embedsignature=True
cdef extern from "math.h":
    double modf (double value, double *iptr)
    
cdef extern from "time.h" nogil:
    ctypedef int time_t
    
    struct timespec:
        time_t  tv_sec
        long    tv_nsec
    int c_nanosleep "nanosleep" (timespec *req, timespec *rem)
    int EINTR
    int EINVAL
    
def nanosleep(double dt):
    """
    sleep for `dt` seconds
    """
    cdef timespec ts
    cdef timespec rem
    cdef double ns, seconds
    cdef int error
    
    ns = modf(dt, &seconds) * 1e9
    ts.tv_sec = <time_t>seconds
    ts.tv_nsec = <long>ns
    
    error = c_nanosleep(&ts, &rem)
    return error