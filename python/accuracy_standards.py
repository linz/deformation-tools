
from collections import namedtuple

# Accuracy standards
#   na = network accuracy (metres)
#   lac = local accuracy (metres)
#   lap = local accuracy (proportional)

AccuracyStandard=namedtuple('AccuracyStandard','name code na lac lap')

accuracies=[
    AccuracyStandard('National reference frame','NRF',0.05,0.003,3.0E-8),
    AccuracyStandard('Deformation monitoring - national','DMN',0.05,0.003,3.0E-7),
    AccuracyStandard('Deformation monitoring - regional','DMR',0.10,0.003,1.0E-6),
    AccuracyStandard('Deformation monitoring - local','DML',0.15,0.01,1.0E-6),
    AccuracyStandard('Basic Geospatial network','BGN',0.15,0.01,5.0E-5),
]

for a in accuracies:
    globals()[a.code]=a

