
import CsvFile

fields = (
    'animal str',
    'legs int',
    'height float',
    'isfish ^(yes|no)$',
    'dob datetime'
)
f = CsvFile.CsvFile('animal','test.csv',fields)
for c in f:
    print c
