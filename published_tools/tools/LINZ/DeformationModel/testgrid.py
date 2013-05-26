import Grid

#file='../../../model/patch_2010_chch/grid_farfield.csv'
#values=(file,171.4,173.2,-44,-42.84,46,30,3)
file='../../../model/patch_2010_chch/grid_nearfield.csv'
values=(file,171.95,172.49,-43.688,-42.432,55,33,3)

testpoints= (
    (173.04,-43.92),
    (172.82,-43.88),
    (172.83,-43.88),
    (173.12,-43.90),
    (173.12,-43.91),
    (173.21,-43.91),
    (171.39,-43.91),
    (172,-44.01),
    (172,-41.99),
)

g = Grid.Grid(*values)
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
