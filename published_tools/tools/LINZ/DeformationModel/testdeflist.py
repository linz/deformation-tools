import DeformationList

d = DeformationList.DeformationList(3,4)
d.addPoint([1,2,3])
d.addPoint([2,4,6])
d.addPoint([3,None,9])
d.addPoint([4,2,16])
d.checkValid()
print d.calcDeformation((0,1),(0.5,0.5))
print d.calcDeformation((0,1,3),(-1,1,-1))
print d.calcDeformation((3,),(2,))
print d.calcDeformation((3,2),(2,1))
