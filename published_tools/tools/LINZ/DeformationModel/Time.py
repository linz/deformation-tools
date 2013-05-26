
import re
from datetime import datetime, date, time
from Error import InvalidValueError

class Time( object ):

    def __init__( self, dt ):
        self._dt = dt

    def __str__( self ):
        return self.strftime()

    def __cmp__( self, dt ):
        dt = Time.Parse(dt)
        if dt == None:
            return 1
        return cmp(self._dt,dt._dt)

    def strftime( self, format='%Y-%m-%d' ):
        return self._dt.strftime(format)

    def daysAfter( self, t0 ):
        td = self._dt-t0._dt
        return td.days + float(td.seconds)/(24 * 3600)

    @staticmethod
    def Now():
        return Time(datetime.now())

    @staticmethod
    def Parse( t ):
        if type(t) == Time:
            return t
        if type(t) == datetime:
            return Time(t)
        if type(t) == date:
            return Time(datetime.combine(t,time(0,0,0)))
        if t == None or t == '' or t == '0':
            return None
        if type(t) not in (str,unicode):
            return InvalidValueError("Invalid date/time "+str(t))
        m = re.match(r'^(\d\d\d\d)\-(\d\d)\-(\d\d)$',t)
        if m:
            return Time(datetime(int(m.group(1)),int(m.group(2)),int(m.group(3)),0,0,0))
        m = re.match(r'^(\d\d\d\d)\-(\d\d)\-(\d\d)\s+(\d\d)\:(\d\d)\:(\d\d)$',t)
        if m:
            return Time(
                datetime(int(m.group(1)),int(m.group(2)),int(m.group(3)),
                int(m.group(4)),int(m.group(5)),int(m.group(6))))
        raise InvalidValueError("Invalid date/time "+str(t))
    
