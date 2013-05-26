import TIN
import sys

dir='../../../model/patch_2009_dusky/'
file1 = dir+'trig_pts_dusky.csv'
file2 = dir+'trig_trg_dusky.csv'
values=(file1,file2,164.5,170.5,-49,-41,24,38,2)
testpoints= (
    (167,-45),
(166.157096,-45.226101),
(166.718923,-42.786949),
(164.800489,-43.389886),
(166.294127,-45.842741),
)

g = TIN.TIN(*values)

if len(sys.argv) > 1:
    f = open(sys.argv[1],"w")
    f.write("id\twkt\n")
    for t in g.triangles():
        f.write(str(t['id']))
        f.write("\t")
        f.write("POLYGON((")
        f.write(",".join([str(t['points'][i][0])+' '+str(t['points'][i][1]) for i in (0,1,2,0)]))
        f.write("))\n")
    f.close()

try:
    g.load()
    for p in testpoints:
        try:
            v = g.calcDeformation(*p)
            print v
        except Exception as e:
            print e
except Exception as e:
    print e
    raise
