import sys
import logging
from Model import Model
from Error import Error

logging.basicConfig(level=logging.INFO)

print "\n" * 10

dm = Model('../../../model',load=True)


testpoints= (
    (167,-45),
(166.157096,-45.226101),
(166.718923,-42.786949),
(164.800489,-43.389886),
(166.294127,-45.842741),
)

testdates=(
    '2000-01-01',
    '2001-04-05',
    '2005-01-20',
    '2008-07-09',
    '2011-10-09',
    '2020-01-01',
)

for x in testpoints:
    print "==================================================="
    print "Calculating at ",x," for ",testdates[-1]
    try:
        print dm.calcDeformation(testdates[-1],x[0],x[1])
    except Error:
        print sys.exc_info()[1]

x = testpoints[0]
for d in testdates:
    print "==================================================="
    print "Calculating at ",x," for ",d
    try:
        print dm.calcDeformation(d,x[0],x[1])
    except Error:
        print sys.exc_info()[1]


