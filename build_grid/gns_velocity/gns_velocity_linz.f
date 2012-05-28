       PROGRAM GNS_VELOCITY
c
c       LINZ modified version
c
c       $Id$
c
c       Chris Crook
c       23 July 2010
c
c       Modifications:
c       1) Changing output units to mm/yr
c       2) Replacing prompts with command line parameters
c       3) Not requiring point ids to be sequential
c       4) Ignoring out of range points and continuing
c 
c       Additionally, 2012
c       5) Adding option to specify input by grid:xxxx parameters
c       6) Adding option for LINZ csv output
c
       implicit real*8 (a-h,o-z)
       parameter (nlong_max=140,nlat_max=60,nout_max=99999)
c
       dimension ival(0:nlong_max,0:nlat_max),idel(0:nlong_max,0:nlat_max),
     1  rlong(0:nlong_max,0:nlat_max),dxlong(0:nlong_max,0:nlat_max),
     2  dylong(0:nlong_max,0:nlat_max),rlat(0:nlong_max,0:nlat_max),
     3  dxlat(0:nlong_max,0:nlat_max),dylat(0:nlong_max,0:nlat_max),
     7  xcap(3,nout_max),dxxcap(3,nout_max),dyxcap(3,nout_max),
     9  xlong(nout_max),ylat(nout_max),
     9  xout(nout_max),yout(nout_max),iout(nout_max)
       common/work/nlong,nlat,ival,idel,rlong,dxlong,dylong,rlat,dxlat,dylat,
     1  nval,kvar,nout,iout,xout,yout,xlong,ylat,xcap,dxxcap,dyxcap
       dimension vary(3,3,nout_max),vely(3,nout_max)
       common/filenames/ solnfile,infile,outfile
       character*256 solnfile, infile, outfile
       logical isgrid
       real*8 glonmin,glonmax,glatmin,glatmax,gdlon,gdlat
       integer gnlon,gnlat
       common /griddef/ isgrid,glonmin,glonmax,gdlon,glatmin,glatmax,gdlat,gnlon,gnlat
       integer ilat,ilon,iiout
       real*8 glat,glon

       if (iargc().lt.3) then
           print *, 'Require parameters: gns_model_file input_point_file output_file'
           print *, 'Optional additional params: euler_pole_lat euler_pole_lon  euler_rotation_rate'
           print *, 'Rotation rate in radians / 10**7 years'
           print *, 'The input file can be instead a grid definition entered as:'
           print *, '  grid:min_lon:min_lat:max_lon:max_lat:ngrid_lon:ngrid_lat'
           print *, 'The output file can be prefixed with csv: to generate a LINZ csv file'
           stop ''
       end if
       call getarg(1,solnfile)
       call getarg(2,infile)
       call getarg(3,outfile)

       isgrid = .FALSE.
       if( infile(1:5) .EQ. "grid:") then
         isgrid = .TRUE.
         do i=1,len(infile)
            if(i .lt. 5 .or. infile(i:i) .eq. ':' ) infile(i:i)=' '
         end do
         read(infile,*) glonmin,glatmin,glonmax,glatmax,gnlon,gnlat
         gdlon=(glonmax-glonmin)/gnlon
         gdlat=(glatmax-glatmin)/gnlat
         gnlon = gnlon+1
         gnlat = gnlat+1
       end if

       call CONVERT_LAT_LONG
       call MAKE_X_FOR_OUTPUT
       ny=nout
       call VELOCITY(ny,vely,vary)

       if( outfile(1:4).eq.'csv:') then
           outfile=outfile(5:len(outfile))
           open(3,file=outfile)
           write(3,41)
           if(isgrid) then
              iiout = 1
              iid = 0
              do ilat=1,gnlat
                 glat = glatmin+(ilat-1)*gdlat
                 do ilon=1,gnlon
                    glon = glonmin+(ilon-1)*gdlon
                    iid=iid+1
                    if(iout(iiout).eq.iid) then
                       write(3,42) glon,glat,vely(1,iiout)*0.6371,vely(2,iiout)*0.6371
                       if( iiout .lt. ny ) iiout=iiout+1
                    else
                       write(3,42) glon,glat,9999.0,9999.0
                    end if
                 end do
              end do

           else
              do i=1,ny
                 write(3,42) xlong(i),ylat(i),vely(1,i)*0.6371,vely(2,i)*0.6371
              end do
           end if
       else
       open(3,file=outfile)
       write(3,31)
       write(3,*)' '
       do i=1,ny
         se1=dsqrt(vary(1,1,i))
         se2=dsqrt(vary(2,2,i))
         if (se1*se2.eq.0.0d0) then
           cor12=0.0d0
         else
           cor12=vary(1,2,i)/(se1*se2)
         end if
         write(3,32)iout(i),ylat(i),xlong(i),
     1    vely(1,i)*637.1,vely(2,i)*637.1,
     2    se1*637.1,se2*637.1,cor12
       end do
       end if

       print *, 'Results written to ',TRIM(outfile)

       stop
32    format(1x,i7,2x,f10.5,2x,f10.5,2x,f10.3,2x,f10.3,2x,f10.3,
     1  2x,f10.3,2x,f8.4)
31    format(1x,'    no.',2x,'   LAT.   ',2x,'   LONG.  ',2x,
     1  'Veast     ',2x,'Vnorth    ',2x,
     2  'SEeast    ',2x,'SEnorth   ',2x,
     3  'Correlation')     
41     format(' longitude,latitude,de,dn')
42     format(' ',f9.4,',',f8.4,',',f11.6,',',f11.6)
       end
c
c
c
c
       SUBROUTINE VELOCITY(ny,vely,vary)
c
c  Input: (1) ny = number of points
c         (4) rlat(1:ny) = standard latitude in degrees of the ny points
c         (5) rlong(1:ny) = standard longitude in degrees of the ny points
c
c  Output: (6) vely(1:3,1:ny) = the 3-component velocity vectors at the ny
c             points. The first component is the velocity in the direction of
c             increasing longitude, the second is the velocity in the
c             direction of increasing latitude, and the third component is the
c             rate of rotation about the radial axis
c          (7) vary(1:3,1:3,1:ny) = the corresponding 3x3 variance-covariance
c             matrices at the ny points
c
c
c
       implicit real*8 (a-h,o-z)
       parameter (nlong_max=140,nlat_max=60,nval_max=5000,
     1  nx_max=3*nval_max,nout_max=99999)
c
       dimension vely(3,ny),vary(3,3,ny)
       dimension ival(0:nlong_max,0:nlat_max),idel(0:nlong_max,0:nlat_max),
     1  rlong(0:nlong_max,0:nlat_max),dxlong(0:nlong_max,0:nlat_max),
     2  dylong(0:nlong_max,0:nlat_max),rlat(0:nlong_max,0:nlat_max),
     3  dxlat(0:nlong_max,0:nlat_max),dylat(0:nlong_max,0:nlat_max),
     7  xcap(3,nout_max),dxxcap(3,nout_max),dyxcap(3,nout_max),
     9  xlong(nout_max),ylat(nout_max),
     9  xout(nout_max),yout(nout_max),iout(nout_max)
       common/work/nlong,nlat,ival,idel,rlong,dxlong,dylong,rlat,dxlat,dylat,
     1  nval,kvar,nout,iout,xout,yout,xlong,ylat,xcap,dxxcap,dyxcap
       dimension a1(3,0:1,0:1,0:1,0:1),
     1  a2(3,0:1,0:1,0:1,0:1),a3(3,0:1,0:1,0:1,0:1),
     2  val(3,nval_max),w0(3),
     3  a(3,nx_max),x(nx_max),
     7  sd(nx_max,nx_max),asd(3,nx_max)
       character*1 ans,yes,no
       equivalence (val,x)
       character*20 pvalstr
       data yes,no/'Y','N'/
c
       if (nval.gt.nval_max) then
         print *,'nval_max ',nval_max,' needs to be enlarged to at least ',nval
         stop 'compiled arrays too small'
       end if

       read(1,*) ((val(i,j),i=1,3),j=1,nval)
c  print *,'Number of independent rotations =',nval
c  print *,' Upper bound on RMS velocity =',varx
c  print *,'  Number of rectangles =',nrect
       if (kvar.eq.0) then
         print *,'No variances and covariances'
c10    print *,'Do you want to continue (Y/N)?'
c      read(*,501) ans
c    if (ans.eq.no) stop
c    if (ans.ne.yes) goto 10
         ans=no
       else
c20    print *,'Do you want the variance-covariance matrices (Y/N)?'
c    read(*,501) ans
c    if ((ans.ne.yes).and.(ans.ne.no)) goto 20
         ans=yes
       end if
       nx=3*nval
       if (ans.eq.yes) then
         do j=1,nx
           read(1,*) (sd(i,j),i=1,nx)
         end do
       end if
c  print *,'Enter the latitude and longitude (in degrees) of the Euler'
c  print *,' pole for the frame of reference and the rotation rate to'
c  print *,'  be removed'
       plat = 0.0
       plong = 0.0
       prot = 0.0
c  read(*,*) plat,plong,prot
       if( iargc().gt.3) then
         call getarg(4,pvalstr)
         read(pvalstr,*) plat
       end if
       if( iargc().gt.4) then
         call getarg(5,pvalstr)
         read(pvalstr,*) plong
       end if
       if( iargc().gt.5) then
         call getarg(6,pvalstr)
         read(pvalstr,*) prot
       end if
       print *, 'Reference frame rotation.'
       print *, 'Euler pole latitude ',plat,' longitude ',plong
       print *, 'Rotation rate ',prot,' radians per 10**7 years'
       conv=datan(1.0d0)/45.0d0
       plat=conv*plat
       plong=conv*plong
       w0(1)=prot*dcos(plat)*dcos(plong)
       w0(2)=prot*dcos(plat)*dsin(plong)
       w0(3)=prot*dsin(plat)
       do i=1,nval
       do k=1,3
         val(k,i)=val(k,i)-w0(k)
       end do
       end do

       ny=nout
       do iy=1,ny
         xcoord=xout(iy)
         ycoord=yout(iy)
         i0=min0(idint(xcoord),nlong-1)
         j0=min0(idint(ycoord),nlat-1)
         i1=i0+1
         j1=j0+1
         xcoord=xcoord-i0
         ycoord=ycoord-j0
         call V_COEFFS(xcoord,ycoord,a1,a2,a3,
     1    xcap(1,iy),dxxcap(1,iy),dyxcap(1,iy))
         call V_RELATE(i0,i1,j0,j1,nlong,nlat,
     1    ival,idel,a1,a2,a3,nx,a,
     2    nlong_max,nlat_max,nx_max)
         do i=1,3
           s=0.0d0
           do j=1,nx
             s=s+a(i,j)*x(j)
           end do
           vely(i,iy)=s
         end do
         if (ans.eq.yes) then
           do j=1,nx
           do i=1,3
             s=0.0d0
             do k=1,nx
               s=s+a(i,k)*sd(k,j)
             end do
             asd(i,j)=s
           end do
           end do
           do j=1,3
           do i=1,j
             s=0.0d0
             do k=1,nx
               s=s+asd(i,k)*asd(j,k)
             end do
             vary(i,j,iy)=s
             vary(j,i,iy)=s
           end do
           end do
         else
           do j=1,3
           do i=1,3
             vary(i,j,iy)=0.0d0
           end do
           end do
         end if
       end do
       return
c
501    format(a)
       end
c
c
c
c
       SUBROUTINE V_RELATE(i0,i1,j0,j1,nlong,nlat,
     1  ival,idel,a1,a2,a3,nx,a,
     2  nlong_max,nlat_max,nx_max)
c
       implicit real*8 (a-h,o-z)
       dimension ival(0:nlong_max,0:nlat_max),idel(0:nlong_max,0:nlat_max),
     2  a1(3,0:1,0:1,0:1,0:1),a2(3,0:1,0:1,0:1,0:1),
     3  a3(3,0:1,0:1,0:1,0:1),a(3,nx_max)
c
       do ix=1,nx
       do iy=1,3
         a(iy,ix)=0.0d0
       end do
       end do
c
       w=1.0d0
       if (j0.eq.0) then
         if (i0.eq.0)
     1    call V_A00(w,i0,j0,ival,idel,
     2      a1(1,0,0,0,0),a2(1,0,0,0,0),
     3      a3(1,0,0,0,0),
     4      a,nlong_max,nlat_max,nx_max)
         if (i0.gt.0)
     1    call V_A_0(w,i0,j0,ival,idel,
     2      a1(1,0,0,0,0),a2(1,0,0,0,0),
     3      a3(1,0,0,0,0),
     4      a,nlong_max,nlat_max,nx_max)
         if (i1.lt.nlong)
     1    call V_A_0(w,i1,j0,ival,idel,
     2      a1(1,0,0,1,0),a2(1,0,0,1,0),
     3      a3(1,0,0,1,0),
     4      a,nlong_max,nlat_max,nx_max)
         if (i1.eq.nlong)
     1    call V_AN0(w,i1,j0,ival,idel,
     2      a1(1,0,0,1,0),a2(1,0,0,1,0),
     3      a3(1,0,0,1,0),
     4      a,nlong_max,nlat_max,nx_max)
       end if
       if (j0.gt.0) then
         if (i0.eq.0)
     1    call V_A0_(w,i0,j0,ival,idel,
     2      a1(1,0,0,0,0),a2(1,0,0,0,0),
     3      a3(1,0,0,0,0),
     4      a,nlong_max,nlat_max,nx_max)
         if (i0.gt.0)
     1    call V_A__(w,i0,j0,ival,idel,
     2      a1(1,0,0,0,0),a2(1,0,0,0,0),
     3      a3(1,0,0,0,0),
     4      a,nlong_max,nlat_max,nx_max)
         if (i1.lt.nlong)
     1    call V_A__(w,i1,j0,ival,idel,
     2      a1(1,0,0,1,0),a2(1,0,0,1,0),
     3      a3(1,0,0,1,0),
     4      a,nlong_max,nlat_max,nx_max)
         if (i1.eq.nlong)
     1    call V_AN_(w,i1,j0,ival,idel,
     2      a1(1,0,0,1,0),a2(1,0,0,1,0),
     3      a3(1,0,0,1,0),
     4      a,nlong_max,nlat_max,nx_max)
       end if
       if (j1.lt.nlat) then
         if (i0.eq.0)
     1    call V_A0_(w,i0,j1,ival,idel,
     2      a1(1,0,0,0,1),a2(1,0,0,0,1),
     3      a3(1,0,0,0,1),
     4      a,nlong_max,nlat_max,nx_max)
         if (i0.gt.0)
     1    call V_A__(w,i0,j1,ival,idel,
     2      a1(1,0,0,0,1),a2(1,0,0,0,1),
     3      a3(1,0,0,0,1),
     4      a,nlong_max,nlat_max,nx_max)
         if (i1.lt.nlong)
     1    call V_A__(w,i1,j1,ival,idel,
     2      a1(1,0,0,1,1),a2(1,0,0,1,1),
     3      a3(1,0,0,1,1),
     4      a,nlong_max,nlat_max,nx_max)
         if (i1.eq.nlong)
     1    call V_AN_(w,i1,j1,ival,idel,
     2      a1(1,0,0,1,1),a2(1,0,0,1,1),
     3      a3(1,0,0,1,1),
     4      a,nlong_max,nlat_max,nx_max)
       end if
       if (j1.eq.nlat) then
         if (i0.eq.0)
     1    call V_A0N(w,i0,j1,ival,idel,
     2      a1(1,0,0,0,1),a2(1,0,0,0,1),
     3      a3(1,0,0,0,1),
     4      a,nlong_max,nlat_max,nx_max)
         if (i0.gt.0)
     1    call V_A_N(w,i0,j1,ival,idel,
     2      a1(1,0,0,0,1),a2(1,0,0,0,1),
     3      a3(1,0,0,0,1),
     4      a,nlong_max,nlat_max,nx_max)
         if (i1.lt.nlong)
     1    call V_A_N(w,i1,j1,ival,idel,
     2      a1(1,0,0,1,1),a2(1,0,0,1,1),
     3      a3(1,0,0,1,1),
     4      a,nlong_max,nlat_max,nx_max)
         if (i1.eq.nlong)
     1    call V_ANN(w,i1,j1,ival,idel,
     2      a1(1,0,0,1,1),a2(1,0,0,1,1),
     3      a3(1,0,0,1,1),
     4      a,nlong_max,nlat_max,nx_max)
       end if
       return
       end
c
c
c
c
       SUBROUTINE V_A00(w,i,j,ival,idel,a1,a2,a3,
     3  a,nlong_max,nlat_max,nx_max)
c
       implicit real*8 (a-h,o-z)
       dimension ival(0:nlong_max,0:nlat_max),idel(0:nlong_max,0:nlat_max),
     1  a1(3,0:1,0:1),a2(3,0:1,0:1),a3(3,0:1,0:1),
     2  a(3,nx_max),DX(0:1,-1:1),DY(0:1,-1:1)
       data DX/0.0,0.0,  1.0,-1.0,  0.0,1.0/
       data DY/0.0,0.0,  1.0,-1.0,  0.0,1.0/
c
       im=i
       ip=i+1
       jm=j
       jp=j+1
       call V_CORE(w,i,j,ival,idel,im,ip,jm,jp,DX,DY,a1,a2,a3,
     1  a,nlong_max,nlat_max,nx_max)
       return
       end
c
c
c
c
       SUBROUTINE V_A_0(w,i,j,ival,idel,a1,a2,a3,
     3  a,nlong_max,nlat_max,nx_max)
c
       implicit real*8 (a-h,o-z)
       dimension ival(0:nlong_max,0:nlat_max),idel(0:nlong_max,0:nlat_max),
     1  a1(3,0:1,0:1),a2(3,0:1,0:1),a3(3,0:1,0:1),
     2  a(3,nx_max),DX(0:1,-1:1),DY(0:1,-1:1)
       data DX/0.0,-0.5,  1.0,0.0,  0.0,0.5/
       data DY/0.0,0.0,  1.0,-1.0,  0.0,1.0/
       im=i-1
       ip=i+1
       jm=j
       jp=j+1
       if (ival(im,j).eq.0) then
         call V_A00(w,i,j,ival,idel,a1,a2,a3,
     1    a,nlong_max,nlat_max,nx_max)
         return
       end if
       if (ival(ip,j).eq.0) then
         call V_AN0(w,i,j,ival,idel,a1,a2,a3,
     1    a,nlong_max,nlat_max,nx_max)
         return
       end if
       if (ival(im,jp).eq.0) then
         call V_A00(w,i,j,ival,idel,a1,a2,a3,
     1    a,nlong_max,nlat_max,nx_max)
         return
       end if
       if (ival(ip,jp).eq.0) then
         call V_AN0(w,i,j,ival,idel,a1,a2,a3,
     1    a,nlong_max,nlat_max,nx_max)
         return
       end if
       call V_CORE(w,i,j,ival,idel,im,ip,jm,jp,DX,DY,a1,a2,a3,
     1  a,nlong_max,nlat_max,nx_max)
       return
       end
c
c
c
c
       SUBROUTINE V_AN0(w,i,j,ival,idel,a1,a2,a3,
     3  a,nlong_max,nlat_max,nx_max)
c
       implicit real*8 (a-h,o-z)
       dimension ival(0:nlong_max,0:nlat_max),idel(0:nlong_max,0:nlat_max),
     1  a1(3,0:1,0:1),a2(3,0:1,0:1),a3(3,0:1,0:1),
     2  a(3,nx_max),DX(0:1,-1:1),DY(0:1,-1:1)
       data DX/0.0,-1.0,  1.0,1.0,  0.0,0.0/
       data DY/0.0,0.0,  1.0,-1.0,  0.0,1.0/
c
       im=i-1
       ip=i
       jm=j
       jp=j+1

       call V_CORE(w,i,j,ival,idel,im,ip,jm,jp,DX,DY,a1,a2,a3,
     1  a,nlong_max,nlat_max,nx_max)
       return
       end
c
c
c
c
       SUBROUTINE V_A0_(w,i,j,ival,idel,a1,a2,a3,
     3  a,nlong_max,nlat_max,nx_max)
c
       implicit real*8 (a-h,o-z)
       dimension ival(0:nlong_max,0:nlat_max),idel(0:nlong_max,0:nlat_max),
     1  a1(3,0:1,0:1),a2(3,0:1,0:1),a3(3,0:1,0:1),
     2  a(3,nx_max),DX(0:1,-1:1),DY(0:1,-1:1)
       data DX/0.0,0.0,  1.0,-1.0,  0.0,1.0/
       data DY/0.0,-0.5,  1.0,0.0,  0.0,0.5/
c
       im=i
       ip=i+1
       jm=j-1
       jp=j+1
       if (ival(i,jm).eq.0) then
         call V_A00(w,i,j,ival,idel,a1,a2,a3,
     1    a,nlong_max,nlat_max,nx_max)
         return
       end if
       if (ival(i,jp).eq.0) then
         call V_A0N(w,i,j,ival,idel,a1,a2,a3,
     1    a,nlong_max,nlat_max,nx_max)
         return
       end if
       if (ival(ip,jm).eq.0) then
         call V_A00(w,i,j,ival,idel,a1,a2,a3,
     1    a,nlong_max,nlat_max,nx_max)
         return
       end if
       if (ival(ip,jp).eq.0) then
         call V_A0N(w,i,j,ival,idel,a1,a2,a3,
     1    a,nlong_max,nlat_max,nx_max)
         return
       end if
       call V_CORE(w,i,j,ival,idel,im,ip,jm,jp,DX,DY,a1,a2,a3,
     1  a,nlong_max,nlat_max,nx_max)
       return
       end
c
c
c
c
       SUBROUTINE V_A__(w,i,j,ival,idel,a1,a2,a3,
     3  a,nlong_max,nlat_max,nx_max)
c
       implicit real*8 (a-h,o-z)
       dimension ival(0:nlong_max,0:nlat_max),idel(0:nlong_max,0:nlat_max),
     1  a1(3,0:1,0:1),a2(3,0:1,0:1),a3(3,0:1,0:1),
     2  a(3,nx_max),DX(0:1,-1:1),DY(0:1,-1:1)
       data DX/0.0,-0.5,  1.0,0.0,  0.0,0.5/
       data DY/0.0,-0.5,  1.0,0.0,  0.0,0.5/
c
       im=i-1
       ip=i+1
       jm=j-1
       jp=j+1
       if (ival(i,jm).eq.0) then
         call V_A_0(w,i,j,ival,idel,a1,a2,a3,
     1    a,nlong_max,nlat_max,nx_max)
         return
       end if
       if (ival(im,j).eq.0) then
         call V_A0_(w,i,j,ival,idel,a1,a2,a3,
     1    a,nlong_max,nlat_max,nx_max)
         return
       end if
       if (ival(ip,j).eq.0) then
         call V_AN_(w,i,j,ival,idel,a1,a2,a3,
     1    a,nlong_max,nlat_max,nx_max)
         return
       end if
       if (ival(i,jp).eq.0) then
         call V_A_N(w,i,j,ival,idel,a1,a2,a3,
     1    a,nlong_max,nlat_max,nx_max)
         return
       end if
       if (ival(im,jm).eq.0) then
         call V_A_0(0.5d0*w,i,j,ival,idel,a1,a2,a3,
     1    a,nlong_max,nlat_max,nx_max)
         call V_A0_(0.5d0*w,i,j,ival,idel,a1,a2,a3,
     1    a,nlong_max,nlat_max,nx_max)
         return
       end if
       if (ival(ip,jm).eq.0) then
         call V_A_0(0.5d0*w,i,j,ival,idel,a1,a2,a3,
     1    a,nlong_max,nlat_max,nx_max)
         call V_AN_(0.5d0*w,i,j,ival,idel,a1,a2,a3,
     1    a,nlong_max,nlat_max,nx_max)
         return
       end if
       if (ival(im,jp).eq.0) then
         call V_A_N(0.5d0*w,i,j,ival,idel,a1,a2,a3,
     1    a,nlong_max,nlat_max,nx_max)
         call V_A0_(0.5d0*w,i,j,ival,idel,a1,a2,a3,
     1    a,nlong_max,nlat_max,nx_max)
         return
       end if
       if (ival(ip,jp).eq.0) then
         call V_A_N(0.5d0*w,i,j,ival,idel,a1,a2,a3,
     1    a,nlong_max,nlat_max,nx_max)
         call V_AN_(0.5d0*w,i,j,ival,idel,a1,a2,a3,
     1    a,nlong_max,nlat_max,nx_max)
         return
       end if
       call V_CORE(w,i,j,ival,idel,im,ip,jm,jp,DX,DY,a1,a2,a3,
     1  a,nlong_max,nlat_max,nx_max)
       return
       end
c
c
c
c
       SUBROUTINE V_AN_(w,i,j,ival,idel,a1,a2,a3,
     3  a,nlong_max,nlat_max,nx_max)
c
       implicit real*8 (a-h,o-z)
       dimension ival(0:nlong_max,0:nlat_max),idel(0:nlong_max,0:nlat_max),
     1  a1(3,0:1,0:1),a2(3,0:1,0:1),a3(3,0:1,0:1),
     2  a(3,nx_max),DX(0:1,-1:1),DY(0:1,-1:1)
       data DX/0.0,-1.0,  1.0,1.0,  0.0,0.0/
       data DY/0.0,-0.5,  1.0,0.0,  0.0,0.5/
c
       im=i-1
       ip=i
       jm=j-1
       jp=j+1
       if (ival(i,jm).eq.0) then
         call V_AN0(w,i,j,ival,idel,a1,a2,a3,
     1    a,nlong_max,nlat_max,nx_max)
         return
       end if
       if (ival(i,jp).eq.0) then
         call V_ANN(w,i,j,ival,idel,a1,a2,a3,
     1    a,nlong_max,nlat_max,nx_max)
         return
       end if
       if (ival(im,jm).eq.0) then
         call V_AN0(w,i,j,ival,idel,a1,a2,a3,
     1    a,nlong_max,nlat_max,nx_max)
         return
       end if
       if (ival(im,jp).eq.0) then
         call V_ANN(w,i,j,ival,idel,a1,a2,a3,
     1    a,nlong_max,nlat_max,nx_max)
         return
       end if
       call V_CORE(w,i,j,ival,idel,im,ip,jm,jp,DX,DY,a1,a2,a3,
     1  a,nlong_max,nlat_max,nx_max)
       return
       end
c
c
c
c
       SUBROUTINE V_A0N(w,i,j,ival,idel,a1,a2,a3,
     3  a,nlong_max,nlat_max,nx_max)
c
       implicit real*8 (a-h,o-z)
       dimension ival(0:nlong_max,0:nlat_max),idel(0:nlong_max,0:nlat_max),
     1  a1(3,0:1,0:1),a2(3,0:1,0:1),a3(3,0:1,0:1),
     2  a(3,nx_max),DX(0:1,-1:1),DY(0:1,-1:1)
       data DX/0.0,0.0,  1.0,-1.0,  0.0,1.0/
       data DY/0.0,-1.0,  1.0,1.0,  0.0,0.0/
c
       im=i
       ip=i+1
       jm=j-1
       jp=j
       call V_CORE(w,i,j,ival,idel,im,ip,jm,jp,DX,DY,a1,a2,a3,
     1  a,nlong_max,nlat_max,nx_max)
       return
       end
c
c
c
c
       SUBROUTINE V_A_N(w,i,j,ival,idel,a1,a2,a3,
     3  a,nlong_max,nlat_max,nx_max)
c
       implicit real*8 (a-h,o-z)
       dimension ival(0:nlong_max,0:nlat_max),idel(0:nlong_max,0:nlat_max),
     1  a1(3,0:1,0:1),a2(3,0:1,0:1),a3(3,0:1,0:1),
     2  a(3,nx_max),DX(0:1,-1:1),DY(0:1,-1:1)
       data DX/0.0,-0.5,  1.0,0.0,  0.0,0.5/
       data DY/0.0,-1.0,  1.0,1.0,  0.0,0.0/
c
       im=i-1
       ip=i+1
       jm=j-1
       jp=j
       if (ival(im,j).eq.0) then
         call V_A0N(w,i,j,ival,idel,a1,a2,a3,
     1    a,nlong_max,nlat_max,nx_max)
         return
       end if
       if (ival(ip,j).eq.0) then
         call V_ANN(w,i,j,ival,idel,a1,a2,a3,
     1    a,nlong_max,nlat_max,nx_max)
         return
       end if
       if (ival(im,jm).eq.0) then
         call V_A0N(w,i,j,ival,idel,a1,a2,a3,
     1    a,nlong_max,nlat_max,nx_max)
         return
       end if
       if (ival(ip,jm).eq.0) then
         call V_ANN(w,i,j,ival,idel,a1,a2,a3,
     1    a,nlong_max,nlat_max,nx_max)
         return
       end if
       call V_CORE(w,i,j,ival,idel,im,ip,jm,jp,DX,DY,a1,a2,a3,
     1  a,nlong_max,nlat_max,nx_max)
       return
       end
c
c
c
c
       SUBROUTINE V_ANN(w,i,j,ival,idel,a1,a2,a3,
     3  a,nlong_max,nlat_max,nx_max)
c
       implicit real*8 (a-h,o-z)
       dimension ival(0:nlong_max,0:nlat_max),idel(0:nlong_max,0:nlat_max),
     1  a1(3,0:1,0:1),a2(3,0:1,0:1),a3(3,0:1,0:1),
     2  a(3,nx_max),DX(0:1,-1:1),DY(0:1,-1:1)
       data DX/0.0,-1.0,  1.0,1.0,  0.0,0.0/
       data DY/0.0,-1.0,  1.0,1.0,  0.0,0.0/
c
       im=i-1
       ip=i
       jm=j-1
       jp=j
       call V_CORE(w,i,j,ival,idel,im,ip,jm,jp,DX,DY,a1,a2,a3,
     1  a,nlong_max,nlat_max,nx_max)
       return
       end
c
c
c
c
       SUBROUTINE V_CORE(w,i,j,ival,idel,im,ip,jm,jp,DX,DY,a1,a2,a3,
     1  a,nlong_max,nlat_max,nx_max)
c
       implicit real*8 (a-h,o-z)
       dimension ival(0:nlong_max,0:nlat_max),idel(0:nlong_max,0:nlat_max),
     1  DX(0:1,-1:1),DY(0:1,-1:1),a1(3,0:1,0:1),a2(3,0:1,0:1),
     2  a3(3,0:1,0:1),a(3,nx_max)
c
       jdmax=idel(i,j)/2
       idmax=mod(idel(i,j),2)
c
c  idel:  0 = no derivs
c    1 = x deriv
c    2 = y deriv
c    3 = both derivs
c
       do ja=jm,jp
       do ia=im,ip
       do jd=0,jdmax
       do id=0,idmax
         coeff=w*DX(id,ia-i)*DY(jd,ja-j)
         if (coeff.ne.0.0d0) then
         do kcomp=1,3
           ix=kcomp+3*(ival(ia,ja)-1)
             iy=1
             a(iy,ix)=a(iy,ix)
     1        +a1(kcomp,id,jd)*coeff
             iy=iy+1
             a(iy,ix)=a(iy,ix)
     1        +a2(kcomp,id,jd)*coeff
             iy=iy+1
             a(iy,ix)=a(iy,ix)
     1        +a3(kcomp,id,jd)*coeff
         end do
         end if
       end do
       end do
       end do
       end do
       return
       end
c
c
c
c
c
       SUBROUTINE V_COEFFS(x,y,a1,a2,a3,xcap,dxxcap,dyxcap)
c
       implicit real*8 (a-h,o-z)
       dimension a1(3,0:1,0:1,0:1,0:1),a2(3,0:1,0:1,0:1,0:1),
     1  a3(3,0:1,0:1,0:1,0:1),xcap(3),dxxcap(3),dyxcap(3),
     2  f_x(0:1,0:1),d_x(0:1,0:1),f_y(0:1,0:1),d_y(0:1,0:1),
     3  p1(3),p2(3)
c
       x1=x
       x2=x**2
       x3=x**3
       f_x(0,0)=1.0d0-3.0d0*x2+2.0d0*x3
       f_x(0,1)=3.0d0*x2-2.0d0*x3
       f_x(1,0)=x1-2.0d0*x2+x3
       f_x(1,1)=-x2+x3
       d_x(0,0)=-6.0d0*x1+6.0d0*x2
       d_x(0,1)=6.0d0*x1-6.0d0*x2
       d_x(1,0)=1.0d0-4.0d0*x1+3.0d0*x2
       d_x(1,1)=-2.0d0*x1+3.0d0*x2

       x1=y
       x2=y**2
       x3=y**3
       f_y(0,0)=1.0d0-3.0d0*x2+2.0d0*x3
       f_y(0,1)=3.0d0*x2-2.0d0*x3
       f_y(1,0)=x1-2.0d0*x2+x3
       f_y(1,1)=-x2+x3
       d_y(0,0)=-6.0d0*x1+6.0d0*x2
       d_y(0,1)=6.0d0*x1-6.0d0*x2
       d_y(1,0)=1.0d0-4.0d0*x1+3.0d0*x2
       d_y(1,1)=-2.0d0*x1+3.0d0*x2

       p2(3)=dsqrt(xcap(1)**2+xcap(2)**2)
       if (p2(3).ne.0.0d0) then
         p1(1)=-xcap(2)/p2(3)
         p1(2)=xcap(1)/p2(3)
       else
         p1(1)=1.0d0
         p1(2)=0.0d0
       end if
       p1(3)=0.0d0
       p2(1)=-xcap(3)*p1(2)
       p2(2)=xcap(3)*p1(1)
       p1dx=0.0d0
       p1dy=0.0d0
       p2dx=0.0d0
       p2dy=0.0d0
       do k=1,3
         p1dx=p1dx+p1(k)*dxxcap(k)
         p1dy=p1dy+p1(k)*dyxcap(k)
         p2dx=p2dx+p2(k)*dxxcap(k)
         p2dy=p2dy+p2(k)*dyxcap(k)
       end do
       area=p1dx*p2dy-p1dy*p2dx
       do jend=0,1
       do iend=0,1
       do jd=0,1
       do id=0,1
         f=f_x(id,iend)*f_y(jd,jend)
         d1=d_x(id,iend)*f_y(jd,jend)*p2dy
     1    -f_x(id,iend)*d_y(jd,jend)*p2dx
         d2=f_x(id,iend)*d_y(jd,jend)*p1dx
     2    -d_x(id,iend)*f_y(jd,jend)*p1dy
         do k=1,3
           if (p2(3).ne.0.0d0) then
             a1(k,id,jd,iend,jend)=f*p2(k)
             a2(k,id,jd,iend,jend)=-f*p1(k)
           else
             a1(k,id,jd,iend,jend)=0.0d0
             a2(k,id,jd,iend,jend)=0.0d0
           end if
           if (area.ne.0.0d0) then
             a3(k,id,jd,iend,jend)=f*xcap(k)
     1        -0.5d0*(d1*p1(k)+d2*p2(k))/area
           else
             a3(k,id,jd,iend,jend)=0.0d0
           end if
         end do
       end do
       end do
       end do
       end do
       return
       end

*******************************************************************************
       SUBROUTINE CONVERT_LAT_LONG
c
       implicit real*8 (a-h,o-z)
       parameter (nlong_max=140,nlat_max=60,ny_max=99999,
     1  nout_max=99999,ntest_max=5)
c
       dimension ival(0:nlong_max,0:nlat_max),idel(0:nlong_max,0:nlat_max),
     1  rlong(0:nlong_max,0:nlat_max),dxlong(0:nlong_max,0:nlat_max),
     2  dylong(0:nlong_max,0:nlat_max),rlat(0:nlong_max,0:nlat_max),
     3  dxlat(0:nlong_max,0:nlat_max),dylat(0:nlong_max,0:nlat_max),
     7  xcap(3),dxxcap(3),dyxcap(3),pos(ntest_max),p1(3),p2(3),
     8  xc0(ny_max),yc0(ny_max),xl0(ny_max),yl0(ny_max),
     9  dxxl0(ny_max),dyxl0(ny_max),dxyl0(ny_max),dyyl0(ny_max),
     9  xout(nout_max),yout(nout_max),iout(nout_max)
c  character*40 name
       common/work/nlong,nlat,ival,idel,rlong,dxlong,dylong,rlat,dxlat,dylat,
     1  nval,kvar,nout,iout,xout,yout
       common/filenames/ solnfile,infile,outfile
       character*256 solnfile, infile, outfile
       logical isgrid
       real*8 glonmin,glonmax,glatmin,glatmax,gdlon,gdlat
       integer gnlon,gnlat
       common /griddef/ isgrid,glonmin,glonmax,gdlon,glatmin,glatmax,gdlat,gnlon,gnlat

c
       print *,'Reading model from ',solnfile(1:LEN(TRIM(solnfile)))
       open(unit=1,file=solnfile,form='formatted',status='old')
       read(1,*) nlong,nlat,nval,kvar

       if (nlong.gt.nlong_max) then 
         print *,'nlong_max ',nlong_max,' needs to be at least ',n_long
         stop 'compiled array sizes too small'
       end if
       if (nlat.gt.nlat_max) then 
         print *,'nlat_max ',nlat_max,' needs to be at least ',n_lat
         stop 'compiled array sizes too small'
       end if

       do j=0,nlat
       do i=0,nlong
         read(1,*) ia,ja,ival(i,j),idel(i,j)
         if (ival(i,j).ne.0) then
           read(1,*) rlong(i,j),dxlong(i,j),dylong(i,j)
           read(1,*) rlat(i,j),dxlat(i,j),dylat(i,j)
         end if
       end do
       end do

c  print *,'Enter ntest the square root of the total number of test'
c  print *,' points that will be used from each rectangle in the'
c  print *,'  curvilinear grid (N.B. nxtest = ntest and nytest = ntest)'
c  read(*,*) ntest
       ntest=5
       do i=1,ntest
         pos(i)=(dfloat(i)-0.5d0)/dfloat(ntest)
       end do

       conv=datan(1.0d0)/45.0d0

       ny=0
       do j1=1,nlat
         j0=j1-1
       do i1=1,nlong
         i0=i1-1
         if (min0(ival(i0,j0),ival(i1,j0),ival(i0,j1),ival(i1,j1)).gt.0)
     1    then
           do j=1,ntest
             ycoord=j0+pos(j)
           do i=1,ntest
             xcoord=i0+pos(i)
c------------------------------------------------------------------------------
         call INIT_X(i0,i1,j0,j1,xcoord-i0,ycoord-j0,
     1    xlong0,ylat0,xcap,dxxcap,dyxcap,
     2    rlong,dxlong,dylong,rlat,dxlat,dylat,
     3    nlong_max,nlat_max)

         alat=dasin(xcap(3))
         if ((xcap(1).ne.0.0d0).or.(xcap(2).ne.0.0d0)) then
           along=datan2(xcap(2),xcap(1))
         else

           stop 'AS POLE ENCOUNTERED USE DIFFERENT ntest'

         end if
         coslong=dcos(along)
         sinlong=dsin(along)
         coslat=dcos(alat)
         sinlat=dsin(alat)
         p1(1)=-sinlong
         p1(2)=coslong
         p1(3)=0.0d0
         p2(1)=-sinlat*coslong
         p2(2)=-sinlat*sinlong
         p2(3)=coslat
         dxxlong=0.0d0
         dyxlong=0.0d0
         dxylat=0.0d0
         dyylat=0.0d0
         do k=1,3
           dxxlong=dxxlong+dxxcap(k)*p1(k)
           dyxlong=dyxlong+dyxcap(k)*p1(k)
           dxylat=dxylat+dxxcap(k)*p2(k)
           dyylat=dyylat+dyxcap(k)*p2(k)
         end do
         dxxlong=dxxlong/coslat
         dyxlong=dyxlong/coslat
c------------------------------------------------------------------------------
               ny=ny+1

               if (ny.le.ny_max) then
                 xc0(ny)=xcoord
                 yc0(ny)=ycoord
                 xl0(ny)=along
                 yl0(ny)=alat
                 dxxl0(ny)=dxxlong
                 dyxl0(ny)=dyxlong
                 dxyl0(ny)=dxylat
                 dyyl0(ny)=dyylat
               end if
           end do
           end do
         end if
       end do
       end do

               if (ny.gt.ny_max) then
               print *,'ny_max ',ny_max,' needs to be at least',ny
               stop 'compiled array sizes too small'
               end if

c  type *,'ny=',ny

c  print *,' longitude values in the following format. The first line'
c  print *,'  has the number of the points. Then each succesive line'
c  print *,'   contains the number of the point and the latitude value'
c  print *,'    and longitude value (degrees)'
c  read(*,*) name
       if( isgrid ) then
         print *,'Calculating on grid',infile
         nin = gnlon*gnlat
       else 
          print *, 'Reading computation points from ',TRIM(infile)
          open(unit=2,file=infile,form='formatted',status='old')
          read(2,*) nin
       end if

       nout=0
       nbadpt=0
       do i=1,nin
         if (isgrid) then
             ia = i
             xlong = glonmin+mod(i-1,gnlon)*gdlon
             ylat = glatmin+((i-1)/gnlon)*gdlat
         else
             read(2,*) ia,ylat,xlong
         end if

c    if (ia.ne.i) stop 'A POINT HAS THE WRONG INDEX'

         xlong=conv*xlong
         ylat=conv*ylat

c  type *,'Starting lat,long=',ylat/conv,xlong/conv

         iymin=1
         dsqmin=(xl0(1)-xlong)**2+(yl0(1)-ylat)**2
         do iy=2,ny
           dsq=(xl0(iy)-xlong)**2+(yl0(iy)-ylat)**2
           if (dsq.lt.dsqmin) then
             iymin=iy
             dsqmin=dsq
           end if
         end do
         xcoord=xc0(iymin)
         ycoord=yc0(iymin)
         along=xl0(iymin)
         alat=yl0(iymin)
         dxxlong=dxxl0(iymin)
         dyxlong=dyxl0(iymin)
         dxylat=dxyl0(iymin)
         dyylat=dyyl0(iymin)

c  type *,'x,y,lat,long=',xcoord,ycoord,alat/conv,along/conv

         dlong=along-xlong
         dlat=alat-ylat
         dsq=dsqmin
         dmax=1.0d0/dfloat(ntest)
         count=1.0d0
c******************************************************************************
10    det=dxxlong*dyylat-dyxlong*dxylat
         if (det.eq.0.0d0) goto 30
         dx=(-dlong*dyylat+dyxlong*dlat)/det
         dy=(-dxxlong*dlat+dlong*dxylat)/det
         dr=dsqrt(dx**2+dy**2)
         cx=dx/dr
         cy=dy/dr
         dmax=dmin1(dmax,dr)
         x0=xcoord
         y0=ycoord
         dlong0=dlong
         dlat0=dlat
         drxl0=-dlong/dr
         dryl0=-dlat/dr
         i0=idint(x0)
         j0=idint(y0)
         i1=i0+1
         j1=j0+1
         dsq0=dsq
c******************************************************************************
20    xcoord=x0+dmax*cx
         ycoord=y0+dmax*cy
         dr=dmax
         if (xcoord.lt.dfloat(i0)) then
           xcoord=dfloat(i0)
           dr=(xcoord-x0)/cx
           ycoord=y0+dr*cy
         end if
         if (xcoord.gt.dfloat(i1)) then
           xcoord=dfloat(i1)
           dr=(xcoord-x0)/cx
           ycoord=y0+dr*cy
         end if
         if (ycoord.lt.dfloat(j0)) then
           ycoord=dfloat(j0)
           dr=(ycoord-y0)/cy
           xcoord=x0+dr*cx
         end if
         if (ycoord.gt.dfloat(j1)) then
           ycoord=dfloat(j1)
           dr=(ycoord-y0)/cy
           xcoord=x0+dr*cx
         end if
         if (dr.eq.0.0) goto 30
c------------------------------------------------------------------------------
         call INIT_X(i0,i1,j0,j1,xcoord-i0,ycoord-j0,
     1    xlong0,ylat0,xcap,dxxcap,dyxcap,
     2    rlong,dxlong,dylong,rlat,dxlat,dylat,
     3    nlong_max,nlat_max)

         alat=dasin(xcap(3))
         if ((xcap(1).ne.0.0d0).or.(xcap(2).ne.0.0d0)) then
           along=datan2(xcap(2),xcap(1))
         else

           stop 'AS POLE ENCOUNTERED USE DIFFERENT ntest'

         end if
         coslong=dcos(along)
         sinlong=dsin(along)
         coslat=dcos(alat)
         sinlat=dsin(alat)
         p1(1)=-sinlong
         p1(2)=coslong
         p1(3)=0.0d0
         p2(1)=-sinlat*coslong
         p2(2)=-sinlat*sinlong
         p2(3)=coslat
         dxxlong=0.0d0
         dyxlong=0.0d0
         dxylat=0.0d0
         dyylat=0.0d0
         do k=1,3
           dxxlong=dxxlong+dxxcap(k)*p1(k)
           dyxlong=dyxlong+dyxcap(k)*p1(k)
           dxylat=dxylat+dxxcap(k)*p2(k)
           dyylat=dyylat+dyxcap(k)*p2(k)
         end do
         dxxlong=dxxlong/coslat
         dyxlong=dyxlong/coslat
c------------------------------------------------------------------------------
         dlong=along-xlong
         dlat=alat-ylat
         drxlong=cx*dxxlong+cy*dyxlong
         drylat=cx*dxylat+cy*dyylat
         drxlav=(dlong-dlong0)/dr
         drylav=(dlat-dlat0)/dr
         hdr=0.5d0*dr
         ddrxl0=(drxlav-drxl0)/hdr
         ddrxl1=(drxlong-drxlav)/hdr
         ddryl0=(drylav-dryl0)/hdr
         ddryl1=(drylat-drylav)/hdr
         ddrsup=dmax1(dabs(2.0d0*ddrxl0-ddrxl1),
     1    dabs(2.0d0*ddrxl1-ddrxl0))
     2    +dmax1(dabs(2.0d0*ddryl0-ddryl1),
     3    dabs(2.0d0*ddryl1-ddryl0))
         dr0=dmax/(1.0+0.5d0*ddrsup*dmax**2/dsqrt(dsq))
         xcoord=x0+dr0*cx
         ycoord=y0+dr0*cy
         dr=dr0
         if (xcoord.lt.dfloat(i0)) then
           xcoord=dfloat(i0)
           dr=(xcoord-x0)/cx
           ycoord=y0+dr*cy
         end if
         if (xcoord.gt.dfloat(i1)) then
           xcoord=dfloat(i1)
           dr=(xcoord-x0)/cx
           ycoord=y0+dr*cy
         end if
         if (ycoord.lt.dfloat(j0)) then
           ycoord=dfloat(j0)
           dr=(ycoord-y0)/cy
           xcoord=x0+dr*cx
         end if
         if (ycoord.gt.dfloat(j1)) then
           ycoord=dfloat(j1)
           dr=(ycoord-y0)/cy
           xcoord=x0+dr*cx
         end if
         dr=(dr+dr0)/2.0d0
         xcoord=x0+dr*cx
         ycoord=y0+dr*cy
         if ((xcoord.eq.x0).and.(ycoord.eq.y0)) goto 30
         if (dmin1(xcoord,ycoord).lt.0.0d0) goto 30
         i0=idint(xcoord)
         j0=idint(ycoord)
         i1=i0+1
         j1=j0+1
         if (min0(ival(i0,j0),ival(i1,j0),ival(i0,j1),ival(i1,j1)).eq.0)
     1    goto 30

c------------------------------------------------------------------------------
         call INIT_X(i0,i1,j0,j1,xcoord-i0,ycoord-j0,
     1    xlong0,ylat0,xcap,dxxcap,dyxcap,
     2    rlong,dxlong,dylong,rlat,dxlat,dylat,
     3    nlong_max,nlat_max)

         alat=dasin(xcap(3))
         if ((xcap(1).ne.0.0d0).or.(xcap(2).ne.0.0d0)) then
           along=datan2(xcap(2),xcap(1))
         else

           stop 'AS POLE ENCOUNTERED USE DIFFERENT ntest'

         end if
         coslong=dcos(along)
         sinlong=dsin(along)
         coslat=dcos(alat)
         sinlat=dsin(alat)
         p1(1)=-sinlong
         p1(2)=coslong
         p1(3)=0.0d0
         p2(1)=-sinlat*coslong
         p2(2)=-sinlat*sinlong
         p2(3)=coslat
         dxxlong=0.0d0
         dyxlong=0.0d0
         dxylat=0.0d0
         dyylat=0.0d0
         do k=1,3
           dxxlong=dxxlong+dxxcap(k)*p1(k)
           dyxlong=dyxlong+dyxcap(k)*p1(k)
           dxylat=dxylat+dxxcap(k)*p2(k)
           dyylat=dyylat+dyxcap(k)*p2(k)
         end do
         dxxlong=dxxlong/coslat
         dyxlong=dyxlong/coslat
c------------------------------------------------------------------------------
c  type *,'x,y,lat,long=',xcoord,ycoord,alat/conv,along/conv

         dlong=along-xlong
         dlat=alat-ylat
         dsq=dlong**2+dlat**2
         if (dsq.ge.dsq0) then
           if (dmax.lt.5.0d-17) goto 30
           dmax=dmax/1.5d0
           goto 20
         end if
         if (1.0d0+dsq.gt.1.0d0) then
           count=1.0d0+count*dsq/dsq0
           if (count.ge.25.0d0) goto 30
           dmax=dr
           goto 10
         end if
         nout=nout+1

         if (nout.le.nout_max) then
           xout(nout)=xcoord
           yout(nout)=ycoord
           iout(nout)=ia
         end if
         goto 40

30    nbadpt=nbadpt+1

40    continue
       end do
       if (nout.gt.nout_max) then
         print *,'nout_max ',nout_max,' needs to be enlarged to at least ',nout
         stop 'compiled arrays too small'
       end if
       if( nbadpt.gt.0 ) print *,'Cannot compute velocity for ',nbadpt,
     1  ' points out of range' 
       return
c
       end
c
c
c
c
       SUBROUTINE INIT_X(i0,i1,j0,j1,x,y,
     1  xlong,ylat,xcap,dxxcap,dyxcap,
     2  rlong,dxlong,dylong,rlat,dxlat,dylat,
     3  nlong_max,nlat_max)
c
       implicit real*8 (a-h,o-z)
       dimension xcap(3),dxxcap(3),dyxcap(3),
     1  rlong(0:nlong_max,0:nlat_max),dxlong(0:nlong_max,0:nlat_max),
     2  dylong(0:nlong_max,0:nlat_max),rlat(0:nlong_max,0:nlat_max),
     3  dxlat(0:nlong_max,0:nlat_max),dylat(0:nlong_max,0:nlat_max),
     4  f_x(0:1,0:1),d_x(0:1,0:1),f_y(0:1,0:1),d_y(0:1,0:1),
     5  p1(3),p2(3)
c
       x1=x
       x2=x**2
       x3=x**3
       f_x(0,0)=1.0d0-3.0d0*x2+2.0d0*x3
       f_x(0,1)=3.0d0*x2-2.0d0*x3
       f_x(1,0)=x1-2.0d0*x2+x3
       f_x(1,1)=-x2+x3
       d_x(0,0)=-6.0d0*x1+6.0d0*x2
       d_x(0,1)=6.0d0*x1-6.0d0*x2
       d_x(1,0)=1.0d0-4.0d0*x1+3.0d0*x2
       d_x(1,1)=-2.0d0*x1+3.0d0*x2

       x1=y
       x2=y**2
       x3=y**3
       f_y(0,0)=1.0d0-3.0d0*x2+2.0d0*x3
       f_y(0,1)=3.0d0*x2-2.0d0*x3
       f_y(1,0)=x1-2.0d0*x2+x3
       f_y(1,1)=-x2+x3
       d_y(0,0)=-6.0d0*x1+6.0d0*x2
       d_y(0,1)=6.0d0*x1-6.0d0*x2
       d_y(1,0)=1.0d0-4.0d0*x1+3.0d0*x2
       d_y(1,1)=-2.0d0*x1+3.0d0*x2

       xlong=0.0d0
       dxxlong=0.0d0
       dyxlong=0.0d0
       ylat=0.0d0
       dxylat=0.0d0
       dyylat=0.0d0
       do j=j0,j1
         ja=j-j0
       do i=i0,i1
         ia=i-i0
         xlong=xlong+f_x(0,ia)*f_y(0,ja)*rlong(i,j)
     1    +f_x(1,ia)*f_y(0,ja)*dxlong(i,j)
     2    +f_x(0,ia)*f_y(1,ja)*dylong(i,j)
         dxxlong=dxxlong+d_x(0,ia)*f_y(0,ja)*rlong(i,j)
     1    +d_x(1,ia)*f_y(0,ja)*dxlong(i,j)
     2    +d_x(0,ia)*f_y(1,ja)*dylong(i,j)
         dyxlong=dyxlong+f_x(0,ia)*d_y(0,ja)*rlong(i,j)
     1    +f_x(1,ia)*d_y(0,ja)*dxlong(i,j)
     2    +f_x(0,ia)*d_y(1,ja)*dylong(i,j)
         ylat=ylat+f_x(0,ia)*f_y(0,ja)*rlat(i,j)
     1    +f_x(1,ia)*f_y(0,ja)*dxlat(i,j)
     2    +f_x(0,ia)*f_y(1,ja)*dylat(i,j)
         dxylat=dxylat+d_x(0,ia)*f_y(0,ja)*rlat(i,j)
     1    +d_x(1,ia)*f_y(0,ja)*dxlat(i,j)
     2    +d_x(0,ia)*f_y(1,ja)*dylat(i,j)
         dyylat=dyylat+f_x(0,ia)*d_y(0,ja)*rlat(i,j)
     1    +f_x(1,ia)*d_y(0,ja)*dxlat(i,j)
     2    +f_x(0,ia)*d_y(1,ja)*dylat(i,j)
       end do
       end do
       conv=datan(1.0d0)/45.0d0
       along=conv*xlong
       alat=conv*ylat
       dxxlong=conv*dxxlong
       dyxlong=conv*dyxlong
       dxylat=conv*dxylat
       dyylat=conv*dyylat
       coslong=dcos(along)
       sinlong=dsin(along)
       coslat=dcos(alat)
       sinlat=dsin(alat)
       xcap(1)=coslat*coslong
       xcap(2)=coslat*sinlong
       xcap(3)=sinlat
       p1(1)=-sinlong
       p1(2)=coslong
       p1(3)=0.0d0
       p2(1)=-sinlat*coslong
       p2(2)=-sinlat*sinlong
       p2(3)=coslat
       do k=1,3
         dxxcap(k)=dxxlong*coslat*p1(k)+dxylat*p2(k)
         dyxcap(k)=dyxlong*coslat*p1(k)+dyylat*p2(k)
       end do
       return
       end

*******************************************************************************
       SUBROUTINE MAKE_X_FOR_OUTPUT
c
       implicit real*8 (a-h,o-z)
       parameter (nlong_max=140,nlat_max=60,nout_max=99999)
c
       dimension ival(0:nlong_max,0:nlat_max),idel(0:nlong_max,0:nlat_max),
     1  rlong(0:nlong_max,0:nlat_max),dxlong(0:nlong_max,0:nlat_max),
     2  dylong(0:nlong_max,0:nlat_max),rlat(0:nlong_max,0:nlat_max),
     3  dxlat(0:nlong_max,0:nlat_max),dylat(0:nlong_max,0:nlat_max),
     7  xcap(3,nout_max),dxxcap(3,nout_max),dyxcap(3,nout_max),
     9  xlong(nout_max),ylat(nout_max),
     9  xout(nout_max),yout(nout_max),iout(nout_max)
       logical pos1,pos2,pos3,pos4,prev_zero,dx_zero,dy_zero
       common/work/nlong,nlat,ival,idel,rlong,dxlong,dylong,rlat,dxlat,dylat,
     1  nval,kvar,nout,iout,xout,yout,xlong,ylat,xcap,dxxcap,dyxcap
c
       nx=3*nval

       conv=datan(1.0d0)/45.0d0

       ny=nout
       prev_zero=.false.
       do iy=1,ny
         xcoord=xout(iy)
         ycoord=yout(iy)

         if ((xcoord.lt.0.0d0).or.(xcoord.gt.dfloat(nlong)))
     1    stop 'xcoord OUTSIDE RANGE BEING CONSIDERED'
         if ((ycoord.lt.0.0d0).or.(ycoord.gt.dfloat(nlat)))
     1    stop 'ycoord OUTSIDE RANGE BEING CONSIDERED'

         i0=min0(idint(xcoord),nlong-1)
         j0=min0(idint(ycoord),nlat-1)
         i1=i0+1
         j1=j0+1
         im=i0-1
         jm=j0-1
         pos1=min0(ival(i0,j0),ival(i1,j0),ival(i0,j1),ival(i1,j1)).gt.0
         pos2=.true.
         pos3=.true.
         pos4=.true.
         if ((xcoord.gt.dfloat(i0)).or.(i0.eq.0)) then
           pos2=.false.
           pos4=.false.
         end if
         if ((ycoord.gt.dfloat(j0)).or.(j0.eq.0)) then
           pos3=.false.
           pos4=.false.
         end if
         if (pos2) pos2=min0(ival(im,j0),ival(i0,j0),ival(im,j1),
     1      ival(i0,j1)).gt.0
         if (pos3) pos3=min0(ival(i0,jm),ival(i1,jm),ival(i0,j0),
     1      ival(i1,j0)).gt.0
         if (pos4) pos4=min0(ival(im,jm),ival(i0,jm),ival(im,j0),
     1    ival(i0,j0)).gt.0
         if (.not.pos1) then

           if (.not.(pos2.or.pos3.or.pos4))
     1      stop 'POINT NOT IN REGION BEING CONSIDERED'

           if (pos2) then
             i1=i0
             i0=im
           else
             j1=j0
             j0=jm
             if (.not.pos3) then
               i1=i0
               i0=im
             end if
           end if
         end if
         call INIT_X(i0,i1,j0,j1,xcoord-i0,ycoord-j0,
     1    xlong0,ylat0,xcap(1,iy),dxxcap(1,iy),dyxcap(1,iy),
     2    rlong,dxlong,dylong,rlat,dxlat,dylat,
     3    nlong_max,nlat_max)

         dx_zero=dmax1(dabs(dxxcap(1,iy)),dabs(dxxcap(2,iy)),
     1    dabs(dxxcap(3,iy))).eq.0.0d0 
         dy_zero=dmax1(dabs(dyxcap(1,iy)),dabs(dyxcap(2,iy)),
     1    dabs(dyxcap(3,iy))).eq.0.0d0
         if (dx_zero.or.dy_zero) then
           if (.not.prev_zero) then
             print *,'WARNING: the following zeroes in the'
             print *,'Jacobian of the co-ordinate system'
             print *,'have been encountered. Derivatives of'
             print *,'the velocity field are not continuous'
             print *,'at these points and, as the'
             print *,'derivatives are, therefore, not'
             print *,'unique, they will not be computed by'
             print *,'the velocity and rate-of-strain'
             print *,'programs. These programs will output'
             print *,'default zero values for velocity-field'
             print *,'derivatives at these points. On the'
             print *,'other hand, as velocity and the'
             print *,'components of the underlying rotation'
             print *,'rate vector W are continuous'
             print *,'everywhere (provided the constraints'
             print *,'mentioned below on assigning values of'
             print *,'W and constraining its derivatives are'
             print *,'not violated), the programs output'
             print *,'meaningful values of velocity and W.'
             print *,'    To get meaningful values of'
             print *,'derivatives for, for instance, contour'
             print *,'plots change the x or y co-ordinate a'
             print *,'little. (In some cases changes in one'
             print *,'co-ordinate will not work, whereas'
             print *,'changes in the other co-ordinate will'
             print *,'work). Also, check that where knot'
             print *,'points in the co-ordinate system'
             print *,'coincide these points have been'
             print *,'assigned the same W vector and that'
             print *,'the appropriate derivative index has'
             print *,'been assigned - if the coincident'
             print *,'points have different x co-ordinates,'
             print *,'check that the x derivatives of W are'
             print *,'constrained to be zero, and, likewise,'
             print *,'if the coincident points have'
             print *,'different y co-ordinates, check that'
             print *,'the y derivatives of W are constrained'
             print *,'to be zero'
           end if
           if (dx_zero) print *,'x co-ordinate singularity at',
     1      ' x=',xcoord,' y=',ycoord
           if (dy_zero) print *,'y co-ordinate singularity at',
     1      ' x=',xcoord,' y=',ycoord
           prev_zero=.true.
         end if

         ylat(iy)=dasin(xcap(3,iy))/conv
         if ((xcap(1,iy).ne.0.0d0).or.(xcap(2,iy).ne.0.0d0)) then
           xlong(iy)=datan2(xcap(2,iy),xcap(1,iy))/conv
         else
           xlong(iy)=0.0d0
         end if
       end do
       return
c
       end
