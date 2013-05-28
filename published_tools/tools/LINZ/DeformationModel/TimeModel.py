
import math
import re
from datetime import datetime, date, time
from Time import Time
from Error import ModelDefinitionError, InvalidValueError

class TimeModel(object):
    '''
    Class to calculate a component time model.  Note that this implementation ignores 
    time zones.  It assumes that all time supplied are UTC times.

    The TimeModel object has a function calcFactor for calculating the
    scale factor to apply to the deformation model component at a specified 
    date/time.
    '''

    daysinyear=365.2425

    def __init__(self,type,f0,t0,f1,t1,decay):
        '''
        Create the time model

        type - one of velocity, step, ramp, decay
        f0   - initial scale factor (time t0)
        t0   - initial date/time
        f1   - final scale factor (time t1)
        t1   - final date/time
        decay - decay rate for model
        '''
        if type not in ('velocity','step','ramp','decay'):
            raise ModelDefinitionError("Invalid temporal model type "+str(type))
        f0 = float(f0) if f0 != '' else None
        f1 = float(f1) if f1 != '' else None
        t0 = Time.Parse(t0)
        t1 = Time.Parse(t1)
        decay = float(decay) if decay != '' else None

        self._description = None
        if type == 'velocity':
            self._description = "velocity model"
            if t0 == None:
                raise ModelDefinitionError("Reference time missing for velocity time model")
            def calc(t):
                t = Time.Parse(t)
                return t.daysAfter(t0)/self.daysinyear
            self.calcFactor = calc
        
        elif type=='step':
            self._description = "step from "+str(f0)+" to "+str(f1)+" at "+str(t0)
            if t0 == None:
                raise ModelDefinitionError("Reference time missing for step time model")
            if f0 == None or f1 == None:
                raise ModelDefinitionError("Initial or final scale factor missing for step time model")
            def calc(t):
                t = Time.Parse(t)
                return f0 if t < t0 else f1
            self.calcFactor = calc

        elif type=="ramp":
            self._description = "ramp from "+str(f0)+" at "+str(t0)+" to "+str(f1)+" at "+str(t1)
            if t0 == None or t1 == None:
                raise ModelDefinitionError("Reference time missing for ramp time model")
            if t0 > t1:
                raise ModelDefinitionError("End time before start time for ramp time model")
            if f0 == None or f1 == None:
                raise ModelDefinitionError("Initial or final scale factor missing for ramp time model")
            vel = (f1-f0)/t1.daysAfter(t0) if t1 > t0 else 0.0
            def calc(t):
                t = Time.Parse(t)
                return (
                    f0 if t <= t0
                    else f1 if t >= t1
                    else f0+t.daysAfter(t0)*vel
                )
            self.calcFactor = calc

        elif type=="decay":
            self._description = "exponential decay (rate "+str(decay)+" from "+str(f0)+" at "+str(t0)+" to "+str(f1)+" at "+str(t1)
            if t0 == None:
                raise ModelDefinitionError("Reference time missing for decay time model")
            if t1 != None and t0 > t1:
                raise ModelDefinitionError("End time before start time for decay time model")
            if f0 == None or f1 == None:
                raise ModelDefinitionError("Initial or final scale factor missing for decay time model")
            if decay == None or decay <= 0:
                raise ModelDefinitionError("Decay rate missing or not greater than 0 for decay time model")

            fdif = f1-f0
            if t1 == None:
                def calc(t):
                    t = Time.Parse(t)
                    return (f0 if t <= t0 else
                            f0+fdif*(1-math.exp(decay*(t0.daysAfter(t)/self.daysinyear)))
                           )
                self.calcFactor = calc
            else:
                fdif /= (1-math.exp(decay*(t0.daysAfter(t1)/self.daysinyear)))
                def calc(t):
                    t = Time.Parse(t)
                    return (f0 if t <= t0 else
                            f1 if t >= t1 else
                            f0+fdif*(1-math.exp(decay*(t0.daysAfter(t)/self.daysinyear)))
                           )
                self.calcFactor = calc
        self.squareVarianceFactor = type == "velocity"

    def __str__( self ):
        return self._description

    def calcFactor( self, t ):
        '''
        Calculate the scale factor at a specific time.  This is replaced by an 
        instance method in the object initiallization
        '''
        raise RuntimeError("Time model not defined")