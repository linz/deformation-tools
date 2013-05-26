
from datetime import date, datetime, timedelta
from TimeModel import TimeModel
from Time import Time

params = ('3.5',Time.Parse('2004-05-08'),'5.1',Time.Parse('2005-12-01'),0.7);
params2 = ('3.5',Time.Parse('2004-05-08'),'5.1',None,0.7);

models = (
    TimeModel('velocity',*params),
    TimeModel('step',*params),
    TimeModel('ramp',*params),
    TimeModel('decay',*params),
    TimeModel('decay',*params2),
)


f=open('factors.csv','wb');
f.write("date,velocity,step,ramp,decay,decay0\n");

td = timedelta( days=3 );
d0 = date(2004,1,1)

for i in range(300):
    f.write(d0.strftime('%Y-%m-%d'))
    for m in models:
        factor = m.calcFactor(d0)
        f.write(",%.4f" %(factor,))
    f.write("\n")
    d0 += td



