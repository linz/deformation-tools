import sys

data=[]
with open(sys.argv[1]) as of:
    header=of.readline().split()
    for i in range(5):
        data.append([float(x) for x in of.readline().split()])

de=data[2][0]-data[1][0]
dn=data[4][1]-data[3][1]

diff=[]
for ic in range(3):
    for ix,ir in enumerate(range(1,5,2)):
        dx=data[ir+1][ix]-data[ir][ix]
        diff.append((data[ir+1][ic+2]-data[ir][ic+2])/dx)

for i in range(6):
    print header[i+5],data[0][i+5],diff[i]*1000000
