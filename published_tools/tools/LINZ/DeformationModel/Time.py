
import re
from datetime import datetime, date, time, timedelta
from Error import InvalidValueError

class Time( object ):

    def __init__( self, dt ):
        if not isinstance(dt,datetime):
            dt=Time.Parse(dt)
            if dt:
                dt=dt._dt
        self._dt = dt

    def __str__( self ):
        return self.strftime()

    def __repr__( self ):
        return "Time('"+self.strftime()+"')"

    def __cmp__( self, dt ):
        dt = Time.Parse(dt)
        if dt is None:
            return 1
        return cmp(self._dt,dt._dt)

    def strftime( self, format='%Y-%m-%d' ):
        return self._dt.strftime(format)

    def daysAfter( self, t0 ):
        td = self._dt-t0._dt
        return td.days + float(td.seconds)/(24 * 3600)

    def asDateTime( self ):
        return self._dt

    def asDate( self ):
        return date(self._dt.year, self._dt.month, self._dt.day )

    def asYear( self ):
        year=self._dt.year
        y0=date(year,1,1)
        ndays=float((date(year+1,1,1)-date(year,1,1)).days)
        return year+self.daysAfter(Time(y0))/ndays

    @staticmethod
    def Now():
        return Time(datetime.now())

    @staticmethod
    def Parse( t ):
        if isinstance(t,Time):
            return t
        if isinstance(t,datetime):
            return Time(t)
        if isinstance(t,float):
            year=int(t)
            frac=t-year
            d0=datetime(year,1,1)
            d1=datetime(year+1,1,1)
            td=d1-d0
            secs=td.days*24*3600+td.seconds
            d0 = d0+timedelta(0,secs*frac)
            return Time(d0)
        if isinstance(t,date):
            return Time(datetime.combine(t,time(0,0,0)))
        if t is None or t == '' or t == '0':
            return None
        if type(t) not in (str,unicode):
            raise InvalidValueError("Invalid date/time "+str(t))
        if t.lower() == 'now':
            return Time.Now()
        m = re.match(r'^(\d\d\d\d)(\-?)(\d\d)\2(\d\d)$',t)
        if m:
            return Time(datetime(int(m.group(1)),int(m.group(3)),int(m.group(4)),0,0,0))
        m = re.match(r'^(\d\d\d\d)\-(\d\d)\-(\d\d)\s+(\d\d)\:(\d\d)\:(\d\d)$',t)
        if m:
            return Time(
                datetime(int(m.group(1)),int(m.group(2)),int(m.group(3)),
                int(m.group(4)),int(m.group(5)),int(m.group(6))))
        raise InvalidValueError("Invalid date/time "+str(t))
    
