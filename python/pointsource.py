
def points( source ):
    if source.startswith('grid:'):
        # Grid format
        try:
            parts=source.split(':')
            if len(parts) != 7 or parts[0] != 'grid':
                raise ValueError('')
            min_lon=float(parts[1])
            min_lat=float(parts[2])
            max_lon=float(parts[3])
            max_lat=float(parts[4])
            nlon=int(parts[5])
            nlat=int(parts[6])
            if max_lon<=min_lon or max_lat<=min_lat or nlon<2 or nlat < 2:
                raise ValueError('')
            dlon = (max_lon-min_lon)/(nlon-1)
            dlat = (max_lat-min_lat)/(nlat-1)
            lat = min_lat-dlat
            for ilat in range(nlat):
                lat += dlat
                lon = min_lon-dlon
                for ilon in range(nlon):
                    lon += dlon
                    yield [lon,lat]
        except:
            print("Invalid grid definition",source)
    else:
        instream = open(source,"rb")
        csvrdr = csv.reader(instream)
        headers = next(reader)
        try:
            lonncol = headers.index('lon')
            latcol = headers.index('lat')
        except:
            raise ValueError(source+' does not have lon and lat columns')
        nrec=1
        for row in csv:
            nrec += 1
            if len(row)>1 or len(row)==1 and row[0]:
                try:
                    yield [float(row[loncol]),float(row[latcol])]
                except:
                    raise "Invalid latitude and longitude at line "+str(nrec)
