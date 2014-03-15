#cython: boundscheck=False
#cython: embedsignature=True
import numpy
cimport numpy
import rtmidi

DTYPE = numpy.float #np.float64
ctypedef numpy.float_t DTYPE_t
ctypedef numpy.int_t INT_t
from time import sleep
NOTEON  = rtmidi.MidiMessage.noteOn
NOTEOFF = rtmidi.MidiMessage.noteOff

def send_midi(numpy.ndarray[DTYPE_t, ndim=1] notes, numpy.ndarray[DTYPE_t, ndim=1] amps, numpy.ndarray[DTYPE_t] channels, sendfunc):
    cdef int i
    cdef DTYPE_t note, amp
    cdef list noteoffs = []
    cdef object noteOn = NOTEON
    cdef object noteOff = NOTEOFF
    for i in range(amps.shape[0] - 1):
        amp = amps[i]
        if amp > 0:
            note = notes[i]
            channel = int(channels[i])
            sendfunc(noteOn(channel, note, amp))
            noteoffs.append(noteOff(channel, note))
    return noteoffs
    
def send_all(double delay, list messages, sendfunc):
    cdef int i
    if delay > 0:
        sleep(delay)
    for i in range(len(messages)):
        sendfunc(messages[i])
    
        
        