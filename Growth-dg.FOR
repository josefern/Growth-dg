c***********************************************************************
c                         G R O W T H - d g            11/5/2021       *
c                                                                      *
c        Compilation with Intel Visual Fortran,                        *
c            (Quickwin Application Project)                            *
c                                                                      *
c  Aim:      3-D gravity inversion for modelling gravity changes as    *
c               aggregated growing bodies for anomalous density        *
c                                                                      *
c  Reference: A free-geometry geodynamic modelling of surface gravity  *
c                changes using Growth-dg software                      *
c             A.G.Camacho, P.Vajda, C.A.Miller, J.Fernández            *
c ---------------------------------------------------------------------*
c                                                                      *
c  INPUT                                                               *
c  -----                                                               *
c  file "Grad.txt": coordinates, altitudes, observed gravity changes   *
c         xp,yp,zp : planar coordinates (metres) and altitude (metres) *
c         dg: observed gravity changes(in microgals)                   *
c                                                                      *
c  OUTPUT                                                              *
c  ------                                                              *
c  file "mod.txt": Resulting 3-D model of the subsurface anomalous     *
c        structures as given by means of the 3-D grid of cells         *
c        Values for several characteristic parameters of the model     *
c                                                                      *
c  file "FIL.txt": regional, local, calculated and residual gravity    *
c        values resulting from inversion                               *
c                                                                      *
c ---------------------------------------------------------------------*
c  DIMENSIONS                                                          *
c  ----------                                                          *
c   ms:  max. number of gravity stations                               *
c   mc:  max. number of prismatic cells for the 3-D grid               *
c   mb:  number of points in file map.bln for map lines                *
c                                                                      *
c***********************************************************************

      use IFQWIN
      use IFPORT

c  Dimensions

      parameter (ms=5000,mc=99999,mb=2000,ml=6,mr=50,may=100)
      type (xycoord) xy
      dimension 
     & xc(mc),yc(mc),zc(mc),tx(mc),ty(mc),tz(mc),df(mc),far(mc),
     & qm(mc),isen(mc),vol(mc),lev(mc), lay(mc), iup(mc),ido(mc), 
     & rr(ms),g(ms),gc(ms),x(ms),y(ms),z(ms),re(ms),w(ms),
     & xn(ms),yn(ms),zn(ms),xk(ms),yk(ms), 
     & cov(mr),swc(mr),      xb(mb),yb(mb),ib(mb),  klev(ml),       
     & zdm(2),vdm(2),szdm(2) , ice(may)
      integer col(21),  nada
      integer*2 u1,u2,v1,v2
      real*8 gxx,gyy,gxg,gyg,gxz,gyz,gzz,gzg,ryy,ryg,
     & px,py,pz,f1g,f11,srg,su,dg,f  
      character*50 obs,mod,fil
      character*120 texto
      character*9 hoy
      character*1 ts1,ts2,ts3,ts4,ts5,ts6,bl,as,le(15)
      data bl/' '/,as/'*'/,
     & le/'a','b','c','d','e','f','g','h','i','j','k','l','m','n','o'/

      call seed(2016)
      nada=setbkcolorrgb(#f0f0f0)
      call clearscreen($gclearscreen)
      nada=initializefonts()
      nada=settextcolorrgb(#000000)

c    Default values for parameters

      fla =  30
      uw  =   0
      na  =   4
      seli=  50  
      blun= 2.2
      mlev=   4
      iud =   1
      dco =   0
      tpc=    3
      fdib= 106  
      pp  = 680
      base=-14200

c   Reading file names

      call dial1(obs,mod,fil,na,fla,tpc,dco,iud,mlev,blun)

c   Reading position of bechmarks and gravity value

      open(1,file=obs)
      ax= 9.d9
      bx=-ax
      ay= ax
      by=-ay
      xd=0.
      yd=0.
      zd=0.
      gm=0.
      i=0
      ns=0
    1 read(1,*,err=2,end=2) xp,yp,zp,dg 
      eg=1
      if(xp.eq.0.and.yp.eq.0.and.zp.eq.0) go to 2
      if(xp.lt.ax) ax=xp
      if(xp.gt.bx) bx=xp
      if(yp.lt.ay) ay=yp
      if(yp.gt.by) by=yp
      i=i+1
      if(i.ge.ms) call dim(ms)
      x(i)=xp
      y(i)=yp
      z(i)=zp  
      g(i)=dg
      w(i)=eg   
      xd=xd+xp
      yd=yd+yp
      zd=zd+zp
      gm=gm+dg
      go to 1
    2 ns=i 
      gm=gm/ns 
      ig0 =1
      igxy=0 
      if(ns.gt.50) then
          ig0=1
          igxy=1
          tpc=30
          na=50
      endif
      if(ns.le.2) then
        write(*,'(//a,a)') ' ***** Error in data file ',obs
        stop
      endif
      close(1)

c    Normalized coordinates
      
      ixm= nint((xd/ns+ax+bx)/300.)*100
      iym= nint((yd/ns+ay+by)/300.)*100.
      izm= nint(zd/ns/10.)*10.
     
      fe=(bx-ax+by-ay)/2./30000. 
      tic=3000.*fe     
      call round(tic,r)   
      tic=r
      zme=-3000*fe
      zdm(1)=0
      zdm(2)=0     
      techo=-9.d9

      do 3 i=1,ns
      x(i)=x(i)-ixm         
      y(i)=y(i)-iym
      z(i)=z(i)-izm
      zn(i)=z(i)/fe
      if(zn(i).gt.techo) techo=zn(i)  
      w(i)=1.
      xk(i)=x(i)/1000.
      yk(i)=y(i)/1000.  
      xn(i)=x(i)/fe
    3 yn(i)=y(i)/fe

c    Coordinate limits

      ax=(ax-ixm)/fe
      bx=(bx-ixm)/fe
      ay=(ay-iym)/fe
      by=(by-iym)/fe

c    Read geometrical parameters for the 3D grid

      iex=3
      cell=pp*fe
   10 call dial2(obs,mod,fil,na,ns,fla,tpc,dco,iud,mlev,cell,blun)
      pp=cell/fe*0.55                                       

c   Creating a 3D grid of cells (for normalized coordinates)

      avm=pp*pp*pp*pp*1.0d-10   
      nc=0
      kr=0
      nay=0
      write(*,'(///10x,a\)') 'Please, wait ' 
    7 ncr=nc
      kr=kr+1
      zs=techo
      ax=ax-4000
      ay=ay-4000
      bx=bx+4000
      by=by+4000
    8 dz=pp*0.3                
      call especapa(zs,avm,ms,ns,xn,yn,zn,dz) 
      if(dz.gt.1) go to 9
         zs=zs-pp*0.5
         if(zs.gt.base) go to 8
    9 if(kr.eq.2) dz=dz*1.7   
      if(kr.eq.3) dz=dz*2.0  
      write(*,'(a1\)') '.'
      call celdascapa(ax,bx,ay,by,zs,dz,avm,ms,ns,xn,yn,zn,
     & mc,nc,xc,yc,zc,tx,ty,tz,df,seli,ncr,isen)
      if(nc.lt.0) then 
        write(*,'(////80x,a,i6)') '***** Warning: Num. of cells >',mc
        write(*,'(/84x,a,i6)') '>> Try a smaller length for cells'
        pp=680
        go to 10
      endif
      zs=zs-dz                      
         if(kr.eq.1.and.nc.gt.0) then
         nay=nay+1
         ice(nay)=nc
         endif
      if(zs.gt.base) go to 8  
      if(kr.lt.iex) go to 7   
      c=0
      d=0
      do 4 i=1,nc
      xc(i)=xc(i)*fe
      yc(i)=yc(i)*fe
      zc(i)=zc(i)*fe
      tx(i)=tx(i)*fe
      ty(i)=ty(i)*fe
      tz(i)=tz(i)*fe
      vol(i)=tx(i)*ty(i)*tz(i)
      c=c+tx(i)+ty(i)+tz(i)
      d=d+abs(tx(i)-ty(i))+abs(tz(i)-ty(i))+abs(tz(i)-tx(i))
    4 continue
      c=c/3./nc
      d=d/3./nc
      siz=nint(c)
      techo=techo*fe 

      call updo(may,nay,ice,mc,nc,xc,yc,zc,tx,ty,tz,iup,ido)

c  Reading map.bln
      
      open(1,file='map.bln') 
      k=1
      l=2
      lx=192
    6 read(1,*,end=11) j
      do 5 i=k,(j+k-1)
      read(1,*) xp,yp
      ix=(xp-ixm)/fe/fdib
      iy=(yp-iym)/fe/fdib
      if(ix.gt.lx.or.ix.lt.-lx.or.iy.gt.lx.or.iy.lt.-lx) l=1
      if(i.ge.mb) call dim(mb)
      xb(i)=ix
      yb(i)=iy
      ib(i)=l
      l=0
    5 continue
      k=k+j
      l=2
      go to 6
   11 close(1)
      nb=k-1

c  Reading options and parameters for inversion

      pp=cell
      call dial3(obs,mod,fil,na,nc,ns,fla,dco,mlev,pp,
     & ig0,igxy,tpc,blun,iud)
      flg=fla*fla*fe*fe*1000   
      if(mlev.gt.ml) mlev=ml
           iud=0 
      gmi=9.d9
      gma=-9.d9
      disg=0
      do 14 i=1,ns         
      if(g(i).lt.gmi) gmi=g(i)
      if(g(i).gt.gma) gma=g(i)
      disg=disg+g(i)**2
   14 continue
      disg=sqrt(disg/ns)  
      rag=gma-gmi
      deg=gma*gmi

c  Calculating design matrix and sensitivity

      avm=siz/50.  
      do 16 j=1,nc
      sv=0.
      pen=vol(j)*6.672e-3
      do 17 i=1,ns
      xd=x(i)-xc(j)
      yd=y(i)-yc(j)
      zd=z(i)-zc(j)                            
      d=xd*xd+yd*yd+zd*zd + 10.*fe 
   17 sv= sv+zd*zd/d/d/d
      sv=sv/ns
      qm(j)=sv*pen*pen/avm/avm         
   16 isen(j)=1.d-7/sqrt(sv)/fe/fe     


c  ------Graphics and headings ------------------------------------------
      call clearscreen($gclearscreen)                                   !
      nada=setcolorrgb(#000000)                                         !
      nada=initializefonts()                                            !
      write(texto,'(a)') 'GROWTH-dg'                                    !
      nada=setfont('t''Times New Roman''h23w11')                        !
      call moveto(140,15,xy)                                            !
      call outgtext(texto)                                              !
      kdib=0
      write(texto,'(a)')                                                !
     & 'Inversion of gravity changes as 3D mass sources'                !
      nada=setfont('t''Times New Roman''h19w8')                         !
      call moveto(30,42,xy)                                             !
      call outgtext(texto)                                              !
      call date(hoy)                                                    !
      write(texto,'(a)') hoy                                            !
      call moveto(170,250,xy)                                           !
      call outgtext(texto)                                              !

      nada=setfont('t''Courier''h17w7')                                 !
      call ventana(18,124,380,142)                                      !
      nada=setcolorrgb(#a0ffff)                                         !
      nada=rectangle($gfillinterior,285,124,330,142)                    !
      nada=setcolorrgb(#000000)                                         !
      write(texto,*)                                                    !
     & ' Num   Num  C.Side  Blun  Lamb  Dep  % cell  Rand'              !
      call moveto(20,90,xy)                                             !
      call outgtext(texto)                                              !
      write(texto,*)                                                    !
     & 'data  cells   (m)    val   val  val  (kg/m3)  val'              !
      call moveto(20,103,xy)                                            !
      call outgtext(texto)                                              !
      c=tpc                                                             !
      if(dco.ne.0) c=dco                                                !
      write(texto,'(i5,i8,i6,f6.1,i6,i5,f8.1,i6)')                      !
     & ns,nc,nint(cell),blun,nint(fla),iud,c,na                         !
      call moveto(20,126,xy)                                            !
      call outgtext(texto)                                              !

      write(texto,'(a,2x,a)') ' Number of',                             !
     -'Dens.   Resid  Filled  Off.G   Reg.G '                           !
      call moveto(25,165,xy)                                            !
      call outgtext(texto)                                              !
      write(texto,'(a,a)') 'Steps Cells ',                              !
     -'(kg/m3) (uGal) (%cell)  (uG)   (uG/km)'                          !
      call moveto(20,180,xy)                                            !
      call outgtext(texto)                                              !
      call ventana(18,202,380,220)                                      !
      nada=setcolorrgb(#000000)                                         !
      write(texto,'(a)') 'Initial steps . . . . . .'                    !
      call moveto(40,204,xy)                                            !
      call outgtext(texto)                                              !
      
      ix=40                                                             !
      iy=315                                                            !
      nada=setcolorrgb(#000000)                                         !
      write(texto,'(a,i4,a)') 'Mean reliev ',izm,' m'                   !
      call moveto(ix+53,iy-2,xy)                                        !
      call outgtext(texto)                                              !
      write(texto,'(a)') 'Sea level'                                    !
      call moveto(ix+53,iy+22,xy)                                       !
      call outgtext(texto)                                              !
      write(texto,'(a)') 'Neg. dens.contr.'                             !
      call moveto(ix+53,iy+46,xy)                                       !
      call outgtext(texto)                                              !
      write(texto,'(a)') 'Pos. dens.contr.'                             !
      call moveto(ix+53,iy+69,xy)                                       !
      call outgtext(texto)                                              !
      call ventana(ix-20,iy-15,ix+42,iy+100)                            !
      nada=setcolorrgb(#d0f0d0)                                         !
      u1=ix-20                                                          !
      u2=ix+42                                                          !
      v1=iy+5                                                           !
      v2=iy+102                                                         !
      nada=rectangle($gfillinterior,u1,v1 ,u2,v2)                       !
      nada=setcolorrgb(#f57070)                                         !
      u1=ix-11                                                          !
      u2=ix-3                                                           !
      v1=iy+50                                                          !
      v2=iy+57                                                          !
      nada=rectangle($gfillinterior,u1,v1,u2,v2)                        !
      nada=setcolorrgb(#f06060)                                         !
      u1=ix+1                                                           !
      v1=iy+50                                                          !
      u2=ix+10                                                          !
      v2=iy+57                                                          !
      nada=rectangle($gfillinterior,u1,v1,u2,v2)                        !
      nada=setcolorrgb(#c05050)                                         !
      u1=ix+14                                                          !
      v1=iy+50                                                          !
      u2=ix+22                                                          !
      v2=iy+57                                                          !
      nada=rectangle($gfillinterior,u1,v1,u2,v2)                        !
      nada=setcolorrgb(#804040)                                         !
      u1=ix+26                                                          !
      v1=iy+50                                                          !
      u2=ix+34                                                          !
      v2=iy+57                                                          !
      nada=rectangle($gfillinterior,u1,v1,u2,v2)                        !
      nada=setcolorrgb(#7070f7)                                         !
      u1=ix-11                                                          !
      v1=iy+74                                                          !
      u2=ix-3                                                           !
      v2=iy+81                                                          !
      nada=rectangle($gfillinterior,u1,v1,u2,v2)                        !
      nada=setcolorrgb(#6060f0)                                         !
      u1=ix+1                                                           !
      v1=iy+74                                                          !
      u2=ix+10                                                          !
      v2=iy+81                                                          !
      nada=rectangle($gfillinterior,u1,v1,u2,v2)                        !
      nada=setcolorrgb(#5050c0)                                         !
      u1=ix+14                                                          !
      v1=iy+74                                                          !
      u2=ix+22                                                          !
      v2=iy+81                                                          !
      nada=rectangle($gfillinterior,u1,v1,u2,v2)                        !
      nada=setcolorrgb(#404080)                                         !
      u1=ix+26                                                          !
      v1=iy+74                                                          !
      u2=ix+34                                                          !
      v2=iy+81                                                          !
      nada=rectangle($gfillinterior,u1,v1,u2,v2)                        !
      nada=setcolorrgb(#f0f000)                                         !
      call moveto(ix-20,iy+30,xy)                                       !
      u1=ix+42                                                          !
      v1=iy+30                                                          !
      nada=lineto(u1,v1)                                                !

      jx1=625                                                           !
      nada=setcolorrgb(#000000)                                         !
      call setgtextrotation(900)                                        !
      write(texto,'(a,i6,a)') 'Tic=',nint(tic),' m'                     !
      call moveto(jx1-218,620,xy)                                       !
      call outgtext(texto)                                              !
      call moveto(jx1+197,620,xy)                                       !
      call outgtext(texto)                                              !
      call moveto(jx1+197,220,xy)                                       !
      call outgtext(texto)                                              !
      write(texto,'(a)') 'A.G. Camacho 2021'                            !
      call moveto(jx1+195,410,xy)                                       !
      call outgtext(texto)                                              !
      call setgtextrotation(000)                                        !

      jz=460                                                            !
      jx=jx1                                                            !
      jy=220                                                            !
      lx=190                                                            !
      ly=190                                                            !
      lz=172                                                            !
      call ventana(jx-lx,jz,jx+lx,jz+lz)                                !
      call ventana(jx-lx,jy-ly,jx+lx,jy+ly)                             !

      jx0=210                                                           !
      call ventana(jx0-lx,jz,jx0+lx,jz+lz)                              !

      col( 1)=#76a076                                                   !
      col( 2)=#80a680                                                   !
      col( 3)=#86b086                                                   !
      col( 4)=#90b690                                                   !
      col( 5)=#96c096                                                   !
      col( 6)=#a0c6a0                                                   !
      col( 7)=#a6d0a6                                                   !
      col( 8)=#b0d6b0                                                   !
      col( 9)=#b6e0b6                                                   !
      col(10)=#c0e6c0                                                   !
      col(11)=#c6f0c6                                                   !
      col(12)=#d0f6d0                                                   !
      
      jx=345                                                            !
      jy=275                                                            !
      nada=setcolorrgb(#000000)                                         !
      write(texto,'(a)') 'Sens.'                                        !
      call moveto(jx-3,jy-34,xy)                                        !
      call outgtext(texto)                                              !
      write(texto,'(a)') 'kg/uGal*e8'                                   !
      if(fe.gt.1) write(texto,'(a)') 'kg/uGal*e7'                       !
      if(fe.gt.10) write(texto,'(a)') 'kg/uGal*e6'                      !
      call moveto(jx-25,jy-20,xy)                                       !
      call outgtext(texto)                                              !
      call ventana(jx+10,jy,jx+20,jy+142)                               !
      do 21 i=1,12                                                      !
      nada=setcolorrgb(#000000)                                         !
      write(texto,'(i4)') nint(i/fe/fe/fe/2.)                           !
      if(fe.gt.1) write(texto,'(i4)') nint(i/fe/fe/fe/2.*10.)           !
      if(fe.gt.10) write(texto,'(i4)') nint(i/fe/fe/fe/2.*100.)         !
      call moveto(jx-25,jy+(i-1)*12-3,xy)                               !
      call outgtext(texto)                                              !
      nada=setcolorrgb(col(i))                                          !
      u1=jx+10                                                          !
      v1=jy+(i-1)*12                                                    !
      u2=jx+20                                                          !
      v2=jy+i*12                                                        !
   21 nada=rectangle($gfillinterior,u1,v1,u2,v2)                        !
      
      jx=jx1                                                            !
      jy=220                                                            !
      nada=setcolorrgb(#c0d0c0)                                         !
      u1=jx-lx                                                          !
      v1=jy-ly                                                          !
      u2=jx+lx                                                          !
      v2=jy+ly                                                          !
      nada=rectangle($gfillinterior,u1,v1,u2,v2)                        !

      zz=techo-8000*fe                                                  !
      do 22 i=1,39                                                      !
      ix=-lx+lx/20.*i                                                   !
      xx=ix*fe*fdib                                                     !
      do 22 l=1,39                                                      !
      iy=-ly+ly/20.*l                                                   !
      yy=iy*fe*fdib                                                     !
      call ss2(ms,ns,x,y,z,fe,xx,yy,zz,sen,k)                           !
      nada=setcolorrgb(col(k))                                          !
      k=10                                                              !
      u1=jx+ix-k                                                        !
      v1=jy-iy-k                                                        !
      u2=jx+ix+k                                                        !
      v2=jy-iy+k                                                        !
   22 nada=rectangle($gfillinterior,u1,v1,u2,v2)                        !

      nada=setcolorrgb(#c0d0c0)                                         !
      u1=jx-lx                                                          !
      v1=jz+30                                                          !
      u2=jx+lx                                                          !
      v2=jz+lz                                                          !
      nada=rectangle($gfillinterior,u1,v1,u2,v2)                        !

      yy=0                                                              !
      do 23 i=1,39                                                      !
      ix=-lx+lx/20.*i                                                   !
      xx=ix*fe*fdib                                                     !
      do 23 l=1,14                                                      !
      iz= ly/20.*l                                                      !
      zz=techo-(iz+30)*fe*fdib                                          !
      call ss2(ms,ns,x,y,z,fe,xx,yy,zz,sen,k)                           !
      nada=setcolorrgb(col(k))                                          !
      k=10                                                              !
      u1=jx+ix-k                                                        !
      v1=jz+30+iz-k                                                     !
      u2=jx+ix+k                                                        !
      v2=jz+30+iz+k                                                     !
   23 nada=rectangle($gfillinterior,u1,v1,u2,v2)                        !

      nada=setcolorrgb(#c0d0c0)                                         !
      jy=jx0                                                            !

      u1=jy-ly                                                          !
      v1=jz+30                                                          !
      u2=jy+ly                                                          !
      v2=jz+lz                                                          !
      nada=rectangle($gfillinterior,u1,v1,u2,v2)                        !

      xx=0                                                              !
      do 24 i=1,39                                                      !
      iy=-ly+ly/20.*i                                                   !
      yy=iy*fe*fdib                                                     !
      do 24 l=1,14                                                      !
      iz= ly/20.*l                                                      !
      zz=techo-(iz+30)*fe*fdib                                          !
      call ss2(ms,ns,x,y,z,fe,xx,yy,zz,sen,k)                           !
      nada=setcolorrgb(col(k))                                          !
      k=10  
      u1=jy+iy-k                                                        !
      v1=jz+30+iz-k                                                     !
      u2=jy+iy+k                                                        !
      v2=jz+30+iz+k                                                     !
   24 nada=rectangle($gfillinterior,u1,v1,u2,v2)                        !
c------------------------------------------------------------------------

      jx2=1036                                                          !
      jx3=1485                                                          !
      jx=jx2                                                            !
      jy=220                                                            !
      call ventana(jx-lx,jy-ly,jx+lx,jy+ly)                             !
      nada=setcolorrgb(#d0f0d0)                                         !
      u1=jx-lx                                                          !
      v1=jy-ly                                                          !
      u2=jx+lx                                                          !
      v2=jy+ly                                                          !
      nada=rectangle($gfillinterior,u1,v1,u2,v2)                        !

      jy= jz+ly                                                         !
      call ventana(jx2-lx,jy-ly,jx2+lx,jy+ly)                           !
      nada=setcolorrgb(#d0f0d0)                                         !
      u1=jx2-lx                                                         !
      v1=jy-ly                                                          !
      u2=jx2+lx                                                         !
      v2=jy+ly                                                          !
      nada=rectangle($gfillinterior,u1,v1,u2,v2)                        !

      jy4=560                                                           !
      call ventana(jx3-lx,jy4-100,jx3+lx,jy4+100)                       !
      u1=jx3-lx-45                                                      !
      v1=jy4-75                                                         !
      u2=u1+8                                                           !
      v2=v1-5                                                           !
      nada=setcolorrgb(#f040f0)                                         !
      nada=ellipse($gfillinterior,u1,v1,u2,v2)                          !
      write(texto,'(a)') 'obs'                                          !
      call moveto(u1+12,v1-12,xy)                                       !
      nada=setcolorrgb(#000000)                                         !
      call outgtext(texto)                                              !
      v1=v1+20                                                          !
      v2=v2+20                                                          !
      nada=setcolorrgb(#00c000)                                         !
      nada=ellipse($gfillinterior,u1,v1,u2,v2)                          !
      write(texto,'(a)') 'mod'                                          !
      call moveto(u1+12,v1-12,xy)                                       !
      nada=setcolorrgb(#000000)                                         !
      call outgtext(texto)                                              !
      write(texto,'(i6,a)') nint(rag), ' uGal'                          !
      call setgtextrotation(900)                                        !
      call moveto(u1+42,v1+100,xy)                                      !
      call outgtext(texto)                                              !
      call setgtextrotation(000)                                        !

      jx=jx3                                                            !
      jy=220                                                            !
      call ventana(jx-lx,jy-ly,jx+lx,jy+ly)                             !
      nada=setcolorrgb(#d0f0d0)                                         !
      u1=jx-lx                                                          !
      v1=jy-ly                                                          !
      u2=jx+lx                                                          !
      v2=jy+ly                                                          !
      nada=rectangle($gfillinterior,u1,v1,u2,v2)                        !
      
      v1=jz+30+izm/fe/fdib                                              !
      nada=setcolorrgb(#f0f000)                                         !
      jx=jx1                                                            !
      call moveto(jx-lx,v1,xy)                                          !
      u1=jx+lx                                                          !
      nada=lineto(u1,v1)                                                !
      jx=jx0                                                            !
      call moveto(jx-lx,v1,xy)                                          !
      u1=jx+lx                                                          !
      nada=lineto(u1,v1)                                                !

      ix=250                                                            !
      iy=355                                                            !
      nada=setcolorrgb(#000000)                                         !
      write(texto,'(i6,a)') nint(tic),' m'                              !
      call moveto(ix-2,iy-2,xy)                                         !
      call outgtext(texto)                                              !
      v1=iy+28                                                          !
      call moveto(ix+3,v1,xy)                                           !
      u1=ix+47                                                          !
      nada=lineto(u1,v1)                                                !
      u1=ix+10                                                          !
      call moveto(u1,v1-5,xy)                                           !
      nada=lineto(u1,v1)                                                !
      u1=ix+10+tic/fe/fdib                                              !
      call moveto(u1,v1-5,xy)                                           !
      nada=lineto(u1,v1)                                                !
      u1=ix-6                                                           !
      call moveto(u1,iy-17,xy)                                          !
      v1=iy+24                                                          !
      nada=lineto(u1,v1)                                                !
      v1=iy+20                                                          !
      call moveto(u1+5,v1,xy)                                           !
      nada=lineto(u1,v1)                                                !
      v1=iy+20-tic/fe/fdib                                              !
      call moveto(u1+5,v1,xy)                                           !
      nada=lineto(u1,v1)                                                !
c-------      
      j=izm/tic+2                                                       !
      do 25 i=1,20                                                      !
      jx=jx1                                                            !
      jy=220                                                            !
      ix=jx-lx                                                          !
      k=(i-1)*tic/fe/fdib                                               !
      v1=jy-ly+k                                                        !
      call moveto(ix,v1,xy)                                             !
      u1=ix+4                                                           !
      if(k.lt.2*ly) nada=lineto(u1,v1)                                  !
      ix=jx+lx                                                          !
      call moveto(ix,v1,xy)                                             !
      u1=ix-4                                                           !
      if(k.lt.2*ly) nada=lineto(u1,v1)                                  !
      ix=jx-lx+k                                                        !
      iy=jy-ly                                                          !
      call moveto(ix,iy,xy)                                             !
      u1=ix                                                             !
      v1=iy+4                                                           !
      if(k.lt.2*lx) nada=lineto(u1,v1)                                  !
      iy=jy+ly                                                          !
      call moveto(ix,iy,xy)                                             !
      v1=iy-4                                                           !
      if(k.lt.2*lx) nada=lineto(u1,v1)                                  !
     
      jx=jx3                                                            !
      jy=220                                                            !
      ix=jx-lx                                                          !
      k=(i-1)*tic/fe/fdib                                               !
      v1=jy-ly+k                                                        !
      u1=ix                                                             !
      call moveto(u1,v1,xy)                                             !
      u1=ix+4                                                           !
      if(k.lt.2*ly) nada=lineto(u1,v1)                                  !
      ix=jx+lx                                                          !
      call moveto(ix,v1,xy)                                             !
      u1=ix-4                                                           !
      if(k.lt.2*ly) nada=lineto(u1,v1)                                  !
      ix=jx-lx+k                                                        !
      iy=jy-ly                                                          !
      call moveto(ix,iy,xy)                                             !
      u1=ix                                                             !
      v1=iy+4                                                           !
      if(k.lt.2*lx) nada=lineto(u1,v1)                                  !
      iy=jy+ly                                                          !
      call moveto(ix,iy,xy)                                             !
      v1=iy-4                                                           !
      if(k.lt.2*lx) nada=lineto(u1,v1)                                  !
      
      k=((i-j)*tic+izm)/fe/fdib                                         !
      v1=jz+30+k                                                        !
      jx=jx1                                                            !
      ix=jx-lx                                                          !
      call moveto(ix,v1,xy)                                             !
      u1=ix+4                                                           !
      if(v1.lt.(jz+lz)) nada=lineto(u1,v1)                              !
      ix=jx+lx                                                          !
      call moveto(ix,v1,xy)                                             !
      u1=ix-4                                                           !
      if(v1.lt.(jz+lz)) nada=lineto(u1,v1)                              !
      k=(i-1)*tic/fe/fdib                                               !
      ix=jx-lx+k                                                        !
      call moveto(ix,jz,xy)                                             !
      u1=ix                                                             !
      v1=jz+5                                                           !
      if(k.lt.2*lx) nada=lineto(u1,v1)                                  !

      k=((i-j)*tic+izm)/fe/fdib                                         !
      v1=jz+30+k                                                        !
      jx=jx0                                                            !
      ix=jx-lx                                                          !
      call moveto(ix,v1,xy)                                             !
      u1=ix+4                                                           !
      if(v1.lt.(jz+lz)) nada=lineto(u1,v1)                              !
      ix=jx+lx                                                          !
      call moveto(ix,v1,xy)                                             !
      u1=ix-4                                                           !
      if(v1.lt.(jz+lz)) nada=lineto(u1,v1)                              !
      k=(i-1)*tic/fe/fdib                                               !
      ix=jx-lx+k                                                        !
      call moveto(ix,jz,xy)                                             !
      u1=ix                                                             !
      v1=jz+5                                                           !
      if(k.lt.2*lx) nada=lineto(u1,v1)                                  !
      
   25 continue                                                          !
      
      write(texto,'(a)') 'Plan view'                                    !
      call moveto(jx1-20,420,xy)                                        !
      call outgtext(texto)                                              !
      write(texto,'(a,55x,a)') 'Model.dg','Obser.dg'                    !
      call moveto(jx2-10,420,xy)                                        !
      call outgtext(texto)                                              !

      write(texto,'(a)') 'Residual values'                              !
      call moveto(jx2-20,855,xy)                                        !
      call outgtext(texto)                                              !
      write(texto,'(a)') 'Elevation view W - E'                         !
      call moveto(590,640,xy)                                           !
      call outgtext(texto)                                              !
      write(texto,'(a)') 'Elevation view S - N'                         !
      call moveto(160,640,xy)                                           !
      call outgtext(texto)                                              !

      col( 1)=#600000                                                   !
      col( 2)=#a00000                                                   !
      col( 3)=#d00000                                                   !
      col( 4)=#f06000                                                   !
      col( 5)=#f09000                                                   !
      col( 6)=#f0c000                                                   !
      col( 7)=#f0e000                                                   !
      col( 8)=#f0f060                                                   !
      col( 9)=#f0f0b0                                                   !
      col(10)=#f3f3d0                                                   !
      col(11)=#f6f6f6                                                   !
      col(12)=#a0f2f2                                                   !
      col(13)=#70e9f0                                                   !
      col(14)=#00e0f0                                                   !
      col(15)=#00c0f0                                                   !
      col(16)=#0090f0                                                   !
      col(17)=#0060f0                                                   !
      col(18)=#0000f0                                                   !
      col(19)=#0000d0                                                   !
      col(20)=#0000a0                                                   !
      col(21)=#000060                                                   !

      ix =jx2+236                                                       !
      ix2=ix                                                            !
      iy2=570                                                           !
      iy=50                                                             !
      pg=rag/20                                                         !
      dg=0                                                              !
      if(abs(gmi).gt.100) then                                          !
       dg=gmi+10*pg                                                     !
       write(texto,'(a)') '+'                                           !
       call moveto(ix-15,iy+275,xy)                                     !
       call outgtext(texto)                                             !
       write(texto,'(i7)') nint(dg)                                     !
       call moveto(ix-43,iy+290,xy)                                     !
       call outgtext(texto)                                             !
       write(texto,'(a)') 'uGal'                                        !
       call moveto(ix-25,iy+305,xy)                                     !
       call outgtext(texto)                                             !
      endif                                                             ! 
      do 26 i=1,21                                                      !
      k=gmi+(i-1)*pg-dg                                                 !
      nada=setcolorrgb(col(i))                                          !
      u1=ix                                                             !
      v1=iy +i*12                                                       !
      u2=ix+13                                                          !
      v2=iy +12+i*12                                                    !
      nada=ellipse($gfillinterior,u1,v1,u2,v2)                          !
      nada=setcolorrgb(#000000)                                         !
      nada=ellipse($gborder      ,u1,v1,u2,v2)                          !
      write(texto,'(i6)') k                                             !
      call moveto(ix-45,iy+i*12-2,xy)                                   !
      call outgtext(texto)                                              !
      u1=ix                                                             !
      v1=iy2+i*12                                                       !
      u2=ix +13                                                         !
      v2=iy2+12+i*12                                                    !
      nada=setcolorrgb(col(i))                                          !
      nada=ellipse($gfillinterior,u1,v1,u2,v2)                          !
      nada=setcolorrgb(#000000)                                         !
      nada=ellipse($gborder      ,u1,v1,u2,v2)                          !
   26 continue                                                          !
      write(texto,'(a)') 'uGal'                                         !
      call moveto(ix-20,iy-7,xy)                                        !
      call outgtext(texto)                                              !
      call moveto(ix-22,iy2-7,xy)                                       !
      call outgtext(texto)                                              !
      write(texto,'(a)') 'Obs.'                                         !
      call moveto(ix-20,iy-20,xy)                                       !
      call outgtext(texto)                                              !
      write(texto,'(a)') 'Res.'                                         !
      call moveto(ix-22,iy2-20,xy)                                      !
      call outgtext(texto)                                              !

      call cont(mb,nb,xb,yb,ib,jx1,220)                                 !
      call cont(mb,nb,xb,yb,ib,jx3,220)                                 !

      l=12                                                              !
      if(ns.gt.100) l=4                                                 !
      do 28 i=1,ns                                                      !
      ix=jx3 + x(i)/fe/fdib                                             !
      iy= 220 - y(i)/fe/fdib                                            !
      k=(g(i)-gmi)/pg+1                                                 !
      nada=setcolorrgb(col(k))                                          !
      u1=ix-l                                                           !
      v1=iy-l                                                           !
      u2=ix+l                                                           !
      v2=iy+l                                                           !
      nada=ellipse($gfillinterior,u1,v1,u2,v2)                          !
      nada=setcolorrgb(#104010)                                         !
      nada=ellipse($gborder,u1,v1,u2,v2)                                !
   28 continue                                                          !
!-----------------------------------------------------------------------!
      
c    More initialization
      
      call step(ms,ns,x,y,rr,dr)
      do 29 i=1,ns
      w(i)=1.
   29 gc(i)=0.
      do 18 j=1,nc
      df(j)=0.
      lev(j)=0
   18 continue
      jend=0
      nlev=0
      nu=0
      ncf=0
      sm=0.
      fit0=9.d29
      fs1=9.d9
      fsa=9.d9
      szdm(1)=0
      szdm(2)=0
      vdm(1)=0.1
      vdm(2)=0.1 
      ncc=nc/na
      sdg=0
      snu=0

c-----------------------------------------------------------------------!
c       Step by step growth of the anom anomalous structure             !
c-----------------------------------------------------------------------!

      fa1 = 0.5    
      fac1= 100.   
  
   30 em=0

c    Previous misfit value  

      fit=fit0*1.005
      if(nu.lt.200 ) fit=fit0*3.00 

c    Lambda value

      cr=1.-fa1*exp(-nu*fac1/nc)    
      fl=flg *cr 

c   Normal equations

      s11=0.
      s1x=0.
      s1y=0.
      s1g=0.
      sxx=0.
      sxy=0.
      syy=0.
      sxg=0.
      syg=0.      
      do 31 i=1,ns
      p2=w(i)
      xd=xk(i)
      yd=yk(i)
      dg=g(i)
      s11=s11+p2
      s1g=s1g+dg*p2
        s1x=s1x+xd*p2
        s1y=s1y+yd*p2
        sxx=sxx+xd*xd*p2
        sxy=sxy+xd*yd*p2
        sxg=sxg+xd*dg*p2
        syy=syy+yd*yd*p2
        syg=syg+yd*dg*p2
   31 continue

c  Exploring model space

      js=0
      r=1./mlev
      do 35 j1=1,ncc
      j=j1
      if(ncc.lt.nc) then
          call RANDOM_NUMBER(ranval)
          j=(nc-1)*ranval+1
      endif
      if(lev(j).ge.1) then 
        if(lev(j).eq.mlev) go to 35
        k=lev(j)+1
        if(klev(k).gt.klev(lev(j))*0.8) go to 35 
        d=r*(mlev-lev(j))
        d=klev(1)*d*d
        if(klev(k).gt.d) go to 35   
      endif

          dep=(zme-zc(j))*iud/1000.*1.e-3

      tt=6.672e-3*vol(j)   
      do 34 k=1,2          
      db=2*k-3                                       
        if(k.eq.1.and.db.ge.-0.001) go to 34
        if(k.eq.2.and.db.le. 0.001) go to 34

         qq=qm(j)*(1.-db*dep)

      srr=0.
      srg=0.
      s1r=0.
         px=0
         py=0
         dg=0
         srx=0.
         sry=0.
         srz=0
      c=10*fe
      do 32 i=1,ns
      xd=x(i)-xc(j)
      yd=y(i)-yc(j)
      zd=z(i)-zc(j)                 
      d=xd*xd+yd*yd+zd*zd+c 
      av=sqrt(zd*zd/d/d/d)*tt     
      rr(i)=gc(i)+av*db    
      zb=rr(i)*w(i)
      s1r=s1r+zb
      srg=srg+zb*g(i)
      if(igxy.ne.1) go to 32
         srx=srx+zb*xk(i)
         sry=sry+zb*yk(i)
   32 srr=srr+zb*rr(i)
      sv=sm+db*db*qq        
      su=srr + fl*sv

         f11=s11*su-s1r*s1r
         f1x=s1x*su-s1r*srx
         f1y=s1y*su-s1r*sry
         f1g=s1g*su-s1r*srg
         if(igxy.eq.1) then
           fxx=sxx*su-srx*srx
           fxy=sxy*su-srx*sry
           fxg=sxg*su-srx*srg
           fyy=syy*su-sry*sry
           fyg=syg*su-sry*srg
           gxx=f11*fxx-f1x*f1x
           gxy=f11*fxy-f1x*f1y
           gxg=f11*fxg-f1x*f1g
           gyy=f11*fyy-f1y*f1y
           gyg=f11*fyg-f1y*f1g
           ryy=gyy*gxx-gxy*gxy
           ryg=gyg*gxx-gxy*gxg
           py=ryg/ryy
           px=(gxg-py*gxy)/gxx
         endif 

c   Solutions: offset and scale factor

         if(ig0.eq.1) dg=(f1g-px*f1x-py*f1y)/f11
         f=(srg-dg*s1r-px*srx-py*sry)/su
         if(f.lt.0.00001) go to 34 

c   Evaluating the solution: misfit value

      sv=f*f*sv*fl
      em=0
      do 33 i=1,ns
      c=g(i)-dg-px*xk(i)-py*yk(i)-rr(i)*f     
   33 em=em+c*c*w(i)                            
      c=sv+em  

c  Keeping best values

      if(c.ge.fit) go to 34
        fit=c
        fs=f
        js=j        
        ks=k
        svs=sv
        qs=qq
         dgs=dg
         pxs=px
         pys=py
   34 continue
   35 continue      
               
c     Growth does stop
      
      iend=1 
      if(js.eq.0) go to 50
      iend=2

c     Levels of density

      lev(js)=lev(js)+1
      if(lev(js).gt.nlev) nlev=lev(js)
      dmax=nlev*fs                        
      dme=0
      k=0
      klm=0
      do 36 i=1,nlev
      if(lev(js).eq.i) then
        klev(i)=klev(i)+1    
        if(i.gt.1) klev(i-1)=klev(i-1)-1
      endif
      if(klev(i).eq.0) go to 36
        if(klev(i).gt.klm) klm=klev(i)
        k=k+klev(i)
        dme=dme+klev(i)*i
   36 continue
      dme=dme/k*fs                    
      c=fs*mlev/2.
      if(dme.lt.c) dme=c
      if(df(js).eq.0) ncf=ncf+1

c    Growth stop

      iend=2
      if(tpc.gt.0.and.nu.gt.6.and.ncf.ge.tpc/100*nc) go to 50
      if(dco.gt.0.and.nu.gt.6.and.dme.le.dco) go to 50      
      iend=0

c     Selected values    

      nu=nu+1
      if(fit.lt.fit0) fit0=fit
      dbs=2*ks-3
      df(js)=df(js)+dbs
      sm=sm+dbs*dbs*qs     
      do 37 i=1,ns
       xd=x(i)-xc(js)
       yd=y(i)-yc(js)
       zd=z(i)-zc(js)                       
       d=xd*xd+yd*yd+zd*zd+10.*fe  
       av=zd/sqrt(d*d*d)*6.672e-3*vol(js)  
   37 gc(i)=gc(i)+av*dbs  
      
       c=vol(js)                
       vdm(ks)=vdm(ks)+c 
       szdm(ks)=szdm(ks)+zc(js)*c
       zdm(ks)=szdm(ks)/vdm(ks)
      
      far(js)=fs
      if(df(js).eq.0) far(js)=0
      if(far(js).gt.99999.) far(js)=99999.

      if(nu.eq.1) fs1=fs
      fsa=fs
      zme=(zdm(1)+zdm(2))/2.

c     Modelled and residual values 
      
      ems=0.
      c=0
      rma=0
      do 41 i=1,ns
      tt=g(i)-dgs-xk(i)*pxs-yk(i)*pys
      re(i)= tt-gc(i)*fs
      rr(i)=abs(re(i))
      if(rr(i).gt.rma) rma=rr(i)   
      c=c+w(i)
   41 ems=ems+w(i)*rr(i)*rr(i) 
      if(ems.lt.0) ems=0
      ems=sqrt(ems/c)
      rma=rma/10

c    New weigthing according previous residual values 

      call rewe(ms,ns,re,blun,w)

c---------Graphics-------------------------------------------------------
      write(texto,'(i7,i6,f7.1,i7,f7.1,i7,i6,i5)')                      !
     & nu,ncf,dme*(2*ks-3),nint(ems),ncf*100./nc,nint(dgs),             !
     & nint(pxs),nint(pys)                                              !
      ly=204                                                            !
      nada=setcolorrgb(#ffffff)                                         !
      v1=ly                                                             !
      v2=ly+17                                                          !
      nada=rectangle($gfillinterior,18,v1,380,v2)                       !
      nada=setcolorrgb(#a0ffff)                                         !
      if(dco.ne.0) nada=rectangle($gfillinterior,110,v1,157,v2)         !
      if(dco.eq.0) nada=rectangle($gfillinterior,207,v1,252,v2)         !
      nada=setcolorrgb(#000000)                                         !
      call moveto(10,ly,xy)                                             !
      call outgtext(texto)                                              !
      lx=190                                                            !
      ly=190                                                            !
      ix=xc(js)/fe/fdib                                                 !
      j=lx-2                                                            !
      if(ix.lt.-j.or.ix.gt.j) go to 30                                  !
      iy=yc(js)/fe/fdib                                                 !
      j=ly-2                                                            !
      if(iy.lt.-j.or.iy.gt.j) go to 30                                  !
      jx=jx1                                                            !
      jy=220                                                            !
      ix=ix+jx                                                          !
      iy=jy-iy                                                          !
      c=0.4                                                             !
      if(mlev.gt.1) c=(lev(js)-1.)/(mlev-1.)                            !
      if(ks.eq.1) then                                                  !
        if(c.ge.0.00) nada=setcolorrgb(#f57070)                         !
        if(c.ge.0.25) nada=setcolorrgb(#f06060)                         !
        if(c.ge.0.50) nada=setcolorrgb(#c05050)                         !
        if(c.ge.0.75) nada=setcolorrgb(#804040)                         !
      endif                                                             !
      if(ks.eq.2) then                                                  !
        if(c.ge.0.00) nada=setcolorrgb(#7070f7)                         !
        if(c.ge.0.25) nada=setcolorrgb(#6060f0)                         !
        if(c.ge.0.50) nada=setcolorrgb(#5050c0)                         !
        if(c.ge.0.75) nada=setcolorrgb(#404080)                         !
      endif                                                             !
      i=nint(tx(js)/2./fe/fdib)                                         !
      u1=ix-i                                                           !
      if(u1.lt.122) u1=122                                              !
      u2=ix+i                                                           !
      if(u2.gt.1028) u2=528                                             !
      j=nint(ty(js)/2./fe/fdib)                                         !
      v1=iy-j                                                           !
      if(v1.lt.22) v1=22                                                !
      v2=iy+j                                                           !
      if(v2.gt.418) v2=418                                              !
      nada=rectangle($gfillinterior,u1,v1,u2,v2)                        !
      iz=-zc(js)/fe/fdib+jz+30                                          !
      k=nint(tz(js)/2./fe/fdib)                                         !
      v1=iz-k                                                           !
      v2=iz+k                                                           !
      if(iz.lt.650) nada=rectangle($gfillinterior,u1,v1,u2,v2)          !
      jy=jx0                                                            !
      iy=jy+yc(js)/fe/fdib                                              !
      u1=iy-i                                                           !
      u2=iy+i                                                           !
      if(iz.lt.650) nada=rectangle($gfillinterior,u1,v1,u2,v2)          !

      if(nu.lt.10) go to 30                                             !
      c=nu/20.                                                          !
      d=c-int(c)                                                        !
      if(nu.gt.10.and.d.ne.0) go to 30                                  !
c-----------------------------------------------------------------------!
 
      ix=393                                                            !
      iy=100                                                            !
      lx=20                                                             !
      ly=280                                                            !
      call ventana(ix+10,iy-30,ix+lx,iy+ly)                             !
      nada=setcolorrgb(#00c0f0)                                         !
      c=ncf/(tpc*nc/100.)                                               !
      if(dco.gt.0) c=(dco/mlev-f1)/(fs-f1)*2                            !
      j=c*(ly+32)                                                       !
      u1=ix+10                                                          !
      v1=iy+ly-j                                                        !
      u2=ix+lx                                                          !
      v2=iy+ly                                                          !
      if(v1.gt.(iy+ly)) v1=iy+ly                                        !
      nada=rectangle($gfillinterior,u1,v1,u2,v2)                        !
      
      nada=setcolorrgb(#f0f0f0)                                         !
      u1=ix2-45                                                         !
      v1=iy2+8                                                          !
      u2=ix2-3                                                          !
      v2=iy2+265                                                        !
      nada=rectangle($gfillinterior,u1,v1,u2,v2)                        !
      nada=setcolorrgb(#000000)                                         !
      do 44 i=1,21                                                      !
      k=(i-11)*rma                                                      !
      write(texto,'(i6)') k                                             !
      call moveto(ix2-45,iy2+i*12,xy)                                   !
      call outgtext(texto)                                              !
   44 continue                                                          !

      lx=190                                                            !
      ly=190                                                            !
      jy=jz+ly                                                          !
      call hor0(jx2,220,lx,lx,tic,fe,fdib,mb,nb,xb,yb,ib)               !
      call hor0(jx2,jy,lx,lx,tic,fe,fdib,mb,nb,xb,yb,ib)                !
      do 45 i=1,ns                                                      !
      ix=jx2  + x(i)/fe/fdib                                            !
      iy= jy  - y(i)/fe/fdib                                            !
      k=11+re(i)/rma                                                    !
      if(k.gt.21) k=21                                                  !
      if(k.lt.1) k=1                                                    !
      nada=setcolorrgb(col(k))                                          !
      l=w(i)*14                                                         !
      if(l.lt.2) l=2                                                    !
      if(ns.gt.100) l=4                                                 !
      u1=ix-l                                                           !
      v1=iy-l                                                           !
      u2=ix+l                                                           !
      v2=iy+l                                                           !
      nada=ellipse($gfillinterior,u1,v1,u2,v2)                          !
      nada=setcolorrgb(#104010)                                         !
      nada=ellipse($gborder,u1,v1,u2,v2)                                !
      ix=jx2  + x(i)/fe/fdib                                            !
      iy= 220 - y(i)/fe/fdib                                            !
      d=gc(i)*fs+xk(i)*pxs+yk(i)*pys+dgs                                !
      k=(d-gmi)/pg+1                                                    !
      if(k.gt.21) k=21                                                  !
      if(k.lt.1) k=1                                                    !
      nada=setcolorrgb(col(k))                                          !
      u1=ix-l                                                           !
      v1=iy-l                                                           !
      u2=ix+l                                                           !
      v2=iy+l                                                           !
      nada=ellipse($gfillinterior,u1,v1,u2,v2)                          !
      nada=setcolorrgb(#104010)                                         !
      nada=ellipse($gborder,u1,v1,u2,v2)                                !
   45 continue                                                          !

      lx=190                                                            !
      ly=100                                                            !
      nada=setcolorrgb(#ffffff)                                         !
      u1=jx3-lx+10                                                      !
      v1=jy4-ly                                                         !
      u2=jx3+lx                                                         !
      v2=jy4+ly                                                         !
      nada=rectangle($gfillinterior,u1,v1,u2,v2)                        !
      nada=setcolorrgb(#808080)                                         !
      c=1.96*ly/rag                                                     !
      if(deg.lt.0) then                                                 !
         iy=jy4+ly-(dgs-gmi)*c                                          !
         u1=jx3-lx+15                                                   !
         v1=iy                                                          !
         call moveto(jx3+lx,iy,xy)                                      !
         nada=lineto(u1,v1)                                             !
      endif                                                             !
      do 42 i=1,ns                                                      !
      tt=gc(i)*fs+xk(i)*pxs+yk(i)*pys+dgs                               !
      ix=jx3 +x(i)/fe/fdib                                              !
      qq=ly-(g(i)-gmi)*c                                                !
      if(abs(qq).le.ly) then                                            !
       iy=jy4+qq                                                        !
       u1=ix-4                                                          !
       v1=iy-2                                                          !
       u2=ix+4                                                          !
       v2=iy+2                                                          !
       nada=setcolorrgb(#f040f0)                                        !
       nada=ellipse($gfillinterior,u1,v1,u2,v2)                         !
      endif                                                             !
      qq=ly-(tt-gmi)*c                                                  !
      if(abs(qq).gt.ly) go to 42                                        !
       iy=jy4+ly-(tt-gmi)*c                                             !
       u1=ix-4                                                          !
       v1=iy-2                                                          !
       u2=ix+4                                                          !
       v2=iy+2                                                          !
       nada=setcolorrgb(#00c000)                                        !
       nada=ellipse($gfillinterior,u1,v1,u2,v2)                         !
   42 continue                                                          !

      iy=820                                                            !
      ix=260                                                            !
      lx=120                                                            !
      ly=138                                                            !
      call ventana(ix-20,iy-ly-8,ix+lx+17,iy+22)                        !
      call histog(ms,ns,w,re,kdib,ix,lx,iy,ly)                          !

      c=(nlev-1.)/8.                                                    !
      ix=54                                                             !
      lx=140/nlev                                                       !
      call ventana(ix-34,iy-146,ix+163,iy+23)                           !
      nada=setcolorrgb(#101040)                                         !
      do 46 i=1,nlev                                                    !
      k=klev(i)                                                         !
      l=k*110./klm                                                      !
      if(l.gt.110) l=110                                                !
      nada=setcolorrgb(#b080b0)                                         !
      u1=ix+2+(i-1)*lx                                                  !
      v1=iy                                                             !
      u2=ix+i*lx                                                        !
      v2=iy-l                                                           !
      nada=rectangle($gfillinterior,u1,v1,u2,v2)                        !
      write(texto,'(i4)') nint(i*fs)                                    !
      nada=setcolorrgb(#000000)                                         !
      k=ix-16+(i-0.5)*lx                                                !
      call moveto(k,iy+2,xy)                                            !
      call outgtext(texto)                                              !
   46 continue                                                          !
      nada=setcolorrgb(#000000)                                         !
      call moveto(ix,iy,xy)                                             !
      u1=ix+165                                                         !
      v1=iy                                                             !
      nada=lineto(u1,v1)                                                !
      call moveto(ix,iy+5,xy)                                           !
      u1=ix                                                             !
      v1=iy-120                                                         !
      nada=lineto(u1,v1)                                                !
      write(texto,'(i5)') klm                                           !
      call moveto(ix-37,iy-115,xy)                                      !
      call outgtext(texto)                                              !
      write(texto,'(i5)') nint(klm*2./3.)                               !
      call moveto(ix-37,iy-80,xy)                                       !
      call outgtext(texto)                                              !
      write(texto,'(i4)') nint(klm/3.)                                  !
      call moveto(ix-30,iy-45,xy)                                       !
      call outgtext(texto)                                              !
      write(texto,'(i1)') 0                                             !
      call moveto(ix-10,iy-8,xy)                                        !
      call outgtext(texto)                                              !
      write(texto,*) 'Num.Cells'                                        !
      call moveto(ix-20,iy-140,xy)                                      !
      call outgtext(texto)                                              !
      write(texto,'(a)') ' Dens.'                                       !
      call moveto(ix+114,iy-95,xy)                                      !
      call outgtext(texto)                                              !
      write(texto,'(a,i3)') 'Ave=',nint(dme)                            !
      call moveto(ix+114,iy-75,xy)                                      !
      call outgtext(texto)                                              !
      write(texto,'(a,i3)') 'Max=',nint(nlev*fs)                        !
      call moveto(ix+114,iy-55,xy)                                      !
      call outgtext(texto)                                              !
      write(texto,'(a,f9.0)') 'kg/m3'                                   !
      call moveto(ix+124,iy-35,xy)                                      !
      call outgtext(texto)                                              !
      if(kdib.eq.0) then                                                !
       write(texto,'(a,i2,a)') 'Rep.cells on',mlev,' dens.levels'       !
       call moveto(ix-25,iy+35,xy)                                      !
       call outgtext(texto)                                             !
      endif                                                             !

      call cov2(ms,ns,x,y,re,w,dr,mr,cov,swc)                           
      ix=465                                                            !
      iy=785                                                            !
      lx=350                                                            !
      call ventana(ix-45,iy-110,ix+lx,iy+55)                            !
      nada=setcolorrgb(#000000)                                         !
      call moveto(ix,iy,xy)                                             !
      u1=ix+lx-2                                                        !
      v1=iy                                                             !
      nada=lineto(u1,v1)                                                !
      call moveto(ix,iy+40,xy)                                          !
      u1=ix                                                             !
      v1=iy-80                                                          !
      nada=lineto(u1,v1)                                                !
      if(kdib.eq.0) then                                                !
      write(texto,'(a,f6.0,a)') 'Autocorrel. of Residues  (Step=',dr,')'!
      call moveto(ix+50,iy+70,xy)                                       !
      call outgtext(texto)                                              !
      endif                                                             !
      write(texto,'(a,f5.2)') 'Autocorrel(1)=',cov(1)                   !
      call moveto(ix+210,iy-105,xy)                                     !
      call outgtext(texto)                                              !
      write(texto,'(a)') 'Horiz.'                                       !
      call moveto(ix+lx-42,iy-45,xy)                                    !
      call outgtext(texto)                                              !
      write(texto,'(a)') 'Dist.'                                        !
      call moveto(ix+lx-35,iy-25,xy)                                    !
      call outgtext(texto)                                              !
      do 47 i=1,10                                                      !
      u1=ix+i*tic/fe/fdib                                               !
      call moveto(u1,iy,xy)                                             !
      v1=iy+6                                                           !
      nada=lineto(u1,v1)                                                !
      if(i.gt.7) go to 47                                               !
      v1=iy-80*(i-3)*0.25                                               !
      u1=ix-5                                                           !
      call moveto(u1,v1,xy)                                             !
      u1=ix                                                             !
      nada=lineto(u1,v1)                                                !
      write(texto,'(f5.2)') -0.75+i*0.25                                !
      call moveto(ix-42,v1-10,xy)                                       !
      call outgtext(texto)                                              !
   47 continue                                                          !
      u1=ix+0.5*dr/fe/fdib                                              !
      v1=iy-cov(1)*80                                                   !
      nada=setcolorrgb(#009000)                                         !
      call moveto(u1,v1,xy)                                             !
      do 48 i=1,mr                                                      !
      u1=ix+(i-0.5)*dr/fe/fdib                                          !
      v1=iy-cov(i)*80                                                   !
      if(abs(cov(i)).eq.0) go to 48                                     !
      k=8*sqrt(swc(i))                                                  !
      if((v1-k).lt.(iy-110).or.(v1+k).gt.(iy+55)) go to 48              !
      if(u1.gt.(ix+lx)) go to 48                                        !
      if(k.eq.0) k=1                                                    !
      nada=lineto(u1,v1)                                                !
      u1=u1-k                                                           !
      v1=v1-k                                                           !
      u2=u1+2*k                                                         !
      v2=v1+2*k                                                         !
      nada=ellipse($gfillinterior,u1,v1,u2,v2)                          !
   48 continue                                                          !
      kdib=1                                                            !
!-----------------------------------------------------------------------

      if(jend.eq.0) go to 30
      go to 60     

c-----------------------------------------------------------------------!
c           End of model growth                                         !
c-----------------------------------------------------------------------!

   50 continue

c-----------------------------------------------------------------------!
c           End of process                                              !
c-----------------------------------------------------------------------!

   60 fs=fsa
      fit0=sqrt(fit0) 
      c=0.
      em=0
      do 51 i=1,ns
      re(i)=g(i)-gc(i)*fs
      em=em+re(i)
   51 gc(i)=abs(re(i))
      em=em/ns
      syg=em/0.6745
      zb=9.d9
      nn=0
      do 52 j=1,nc
      if(zc(j).lt.zb) zb=zc(j)-tz(j)/1.5
      if(df(j).ne.0) nn=nn+1
   52 continue 
      open(1,file=mod)
      kk=0
      write(1,'(i6,5x,a10,f6.0,4x,a6,f6.1,4x,a8,f4.1,4x,a5,f5.1)') 
     & nn,'cell side=',cell,'lamda=',fla,'%-model=',tpc,'dens=',dco
      write(1,*) '  X(m)   Y(m)   Z(m asl)   sx(m)  sy(m)   sz(m)',
     & '  den(kg/m3)   sen    ord'

      srx=0.
      sry=0.
      px=0.
      py=0.
      srg=0.
          nay=0
          z0=9.d9
      do 54 j=1,nc
      if(df(j).eq.0) go to 54
      df(j)=df(j)*fs             
      dg=vol(j)*df(j)                                     ! mass
      c=abs(dg)                                           ! absol mass
      srg=srg+c                                           ! addition absol. mass
      if(df(j).lt.0) then
        srx=srx+c                                         ! addition neg. mass
        px=px+c*zc(j)                                     ! mean depth
      endif
      if(df(j).gt.0) then
        sry=sry+c                                         ! addition pos. mass
        py=py+c*zc(j)                                     ! mean depth
      endif
      call ss2(ms,ns,x,y,z,fe,xc(j),yc(j),zc(j),sen,k)
      write(1,'(i7,2i8,2x,3i7,f9.1,2x,f9.1,f7.1)') 
     -nint(xc(j)+ixm),nint(yc(j)+iym),nint(zc(j)+izm),
     -nint(tx(j)),nint(ty(j)),nint(tz(j)),df(j),
     -sen,far(j)
           z1=zc(j)
           if(z1.gt.z0) go to 54
           if((z1+1).lt.z0) then
            nay=nay+1
            z0=z1
            lev(nay)=z1
            qm(nay)=0
           endif
           qm(nay)=qm(nay)+df(j)*tx(j)*ty(j)  
   54 continue 
      zme=(px+py)/(srx+sry+1.)+izm
      if(srx.gt.0.) px=px/srx+izm                         ! weighted neg. depth
      if(sry.gt.0.) py=py/sry+izm                         ! weighted pos. depth        do 53 i=1,nay
   53 qm(i)=qm(i)/fe/fe/1.e6
       
c  Writing modelled and residual values

      open(2,file=obs)
      open(3,file=fil)
      write(3,*) '    X       Y         Off      dG    Mod-dg    Res-dg'
     & ,'   S-Res  Test'
      nms=0
      dtg=0.
      j1=9999
      j2=-j1
      disr=0
      disl=0
   55 read(2,*,end=57) xp,yp,zp,dg,eg
      if(xp.eq.0.and.yp.eq.0.and.zp.eq.0) go to 57
      xp=xp-ixm
      yp=yp-iym
      zp=zp-izm
      gloc=dg
      dtg=dtg+gloc*gloc
      nms=nms+1
      gcal=0
      do 56 j=1,nc
      xd=xp-xc(j)
      yd=yp-yc(j)
      zd=zp-zc(j)                        
      d=xd*xd+yd*yd+zd*zd+10.*fe 
      av2=zd*zd/d/d/d
      av=sqrt(av2)*6.672e-3*vol(j)
   56 gcal=gcal+av*df(j)
      gres=gloc-gcal
      disl=disl+gloc*gloc
      disr=disr+(greg-dgs)**2
      gxg=abs(gres)/syg
      ts1=bl
      ts2=bl
      ts3=bl
      ts4=bl
      ts5=bl 
      ts6=bl
      if(gxg.gt.1.0*blun) ts1=as
      if(gxg.gt.1.5*blun) ts2=as
      if(gxg.gt.2.0*blun) ts3=as
      if(gxg.gt.2.5*blun) ts4=as      
      if(gxg.gt.3.0*blun) ts5=as      
      if(gxg.gt.4.0*blun) ts6=as            
      write(3,203) xp+ixm,yp+iym,greg,gloc,gcal,gres,gxg,
     -  ts1,ts2,ts3,ts4,ts5,ts6
  203 format(f8.0,f9.0,4f9.1,f9.2,2x,6a1)
      if(gres.gt.j2) j2=gres
      if(gres.lt.j1) j1=gres
      go to 55
   57 close(3)
      close(2)
      dtg=sqrt(dtg/nms)
      disl=sqrt(disl/nms)
      disr=sqrt(disr/nms)
          
c  Writing some characteristic results in MOD.txt

      close(1)    

c----------------------- Final drawing ----------------------------------
      
      ix=jx3-100                                                        !
      call ventana(ix,870,ix+220,890)                                   !
      nada=setcolorrgb(#000000)                                         !
      write(texto,'(a)') 'Process completed'                            !
      call moveto(ix+40,845,xy)                                         !
      call outgtext(texto)                                              !
      if(iend.eq.1) write(texto,'(a)') 'No further suitable cells'      !
      if(iend.eq.2) write(texto,'(a)') 'OK  End value (den,%) reached'  !
      call moveto(ix+5,870,xy)                                          !
      call outgtext(texto)                                              !

      write(texto,'(a)') ' ===> Press (Enter)'                          !
      call moveto(ix+45,900,xy)                                         !
      call outgtext(texto)                                              !

      ix=jx3-205                                                        !
      iy=675                                                            !
      call ventana(ix+130,iy+30,ix+330,iy+70)                           !
      call ventana(ix+130,iy+110,ix+365,iy+126)                         !
      nada=setcolorrgb(#000000)                                         !
      write(texto,*) 'Tot.      Neg.     Posit.'                        !
      call moveto(ix+150,iy+10,xy)                                      !
      call outgtext(texto)                                              !
      write(texto,'(a,3e10.2)') '    Anom.mass(kg)',srg,srx,sry         !
      call moveto(ix,iy+30,xy)                                          !
      call outgtext(texto)                                              !
      write(texto,'(a,3f10.0)') '    Aver.depth(m)',zme,px,py           !
      call moveto(ix,iy+50,xy)                                          !
      call outgtext(texto)                                              !
      write(texto,'(a)') 'Ab.sea l.'                                    !
      call moveto(ix+338,iy+50,xy)                                      !
      call outgtext(texto)                                              !
      write(texto,'(a)') '  Tot.   Off.  Loc.  Resid.  Re/Lo'           !
      call moveto(ix+135,iy+90,xy)                                      !
      call outgtext(texto)                                              !
      write(texto,'(a,f8.0,2f7.0,f6.0,f6.2)') '    Scatter.(uGal)',     !
     & disg,disr,disl,ems, ems/dtg                                      !
      call moveto(ix,iy+110,xy)                                         !
      call outgtext(texto)                                              !
      if(ns.le.50) go to 70                                             !
      nada=setcolorrgb(#f0f0f0)                                         !
      u1=410                                                            !
      u2=u1+12                                                          !
      v1=jz-5                                                           !
      v2=jz+180                                                         !
      nada=rectangle($gfillinterior,u1,v1,u2,v2)                        !
                  cmax=0
      do 58 i=2,nay                                                     !
      c=qm(i-1)-qm(i)                                                   !
      if(c.gt.0) nada=setcolorrgb(#8080ff)                              !
      if(c.lt.0) nada=setcolorrgb(#ff8080)                              !
                  if(c.gt.cmax) cmax=c
      iz=-lev(i)/fe/fdib+jz+30                                          !
      k=3                                                               !
      v1=iz-k                                                           !
      v2=iz+k                                                           !
      nada=rectangle($gfillinterior,u1,v1,u2,v2)                        !
   58 continue 
c-----------------------------------------------------------------------!
      go to 70

C  Some messages for error or warning

   61 nada=setcolorrgb(#000000)                                         !
      write(texto,'(a)')  '*** Error: Unexpected execution end !!!'
      call moveto(650,740,xy)
      call outgtext(texto)
      if(nu.le.1) write(texto,'(a)')
     &'Verify data suitability (units, magnitudes) and model hypothesis'
      call moveto(660,760,xy)
      if(nu.le.1) call outgtext(texto)
      go to 99

   62 write(texto,'(3a)') '*** Error reading file ',obs,' !!'
      call moveto(650,740,xy)
      call outgtext(texto)

c ------------------ Drawing several sections --------------------------!
      iv=0                                                              !
   70 read(*,'(i1)') i                                                  !
      zi=0                                                              !
      do 71 i=1,nc                                                      !
      c=abs(df(i))                                                      !
   71 if(c.gt.zi) zi=c                                                  !
      pa=zi/10.                                                         !
      lx=190                                                            !
      ly=190                                                            !
      lz=180                                                            !
      seli=10                                                           !
      call clearscreen($gclearscreen)                                   !
      nada=setcolorrgb(#000000)                                         !
      nada=initializefonts()                                            !
      write(texto,'(a)') 'GROWTH-dg'                                    !
      nada=setfont('t''Times New Roman''h23w11')                        !
      call moveto(160,15,xy)                                            !
      call outgtext(texto)                                              !
      write(texto,'(a)')                                                !
     & 'Inversion of gravity changes as 3D mass sources'                !
      nada=setfont('t''Times New Roman''h19w9')                         !
      call moveto(50,42,xy)                                             !
      call outgtext(texto)                                              !
      nada=setfont('t''Times New Roman''h16w07')                        !
      write(texto,'(a)') hoy                                            !
      call moveto(200,70,xy)                                            !
      call outgtext(texto)                                              !
      l=0                                                               !
      do 73 i=1,nc                                                      !
      k=nint((df(i)+zi)/pa)+1                                           !
      if(rr(k).ne.-99) l=l+1                                            !
      rr(k)=-99                                                         !
   73 continue                                                          !
 
      ix=360                                                            !
      iy=825                                                            !
      nada=setfont('t''Times New Roman''h14w06')                        !
      i=0                                                               !
      j=0                                                               !
      do 72 k=1,21                                                      !
      if(rr(k).ne.-99) go to 72                                         !
      c=-zi+(k-1)*pa                                                    !
      if(c.gt.0.and.j.eq.0) i=0                                         !
      if(c.gt.0.and.j.eq.0) j=1                                         !
      i=i+1                                                             !
      nada=setcolorrgb(col(k))                                          !
      u1=ix-3                                                           !
      if(c.gt.0) u1=ix+65                                               !
      v1=iy+i*14                                                        !
      u2=u1+16                                                          !
      v2=v1+13                                                          !
      nada=setcolorrgb(#d0f0d0)                                         !
      nada=rectangle($gfillinterior,u1,v1,u2,v2)                        !
      nada=setcolorrgb(col(k))                                          !
      u1=u1+3                                                           !
      v1=v1+1                                                           !
      u2=u2+3                                                           !
      v2=v2-1                                                           !
      nada=rectangle($gfillinterior,u1,v1,u2,v2)                        !
      nada=setcolorrgb(#000000)                                         !
      write(texto,'(f6.0)') c                                           !
      call moveto(u1-37,v1,xy)                                          !
      call outgtext(texto)                                              !
   72 continue                                                          !
      nada=setfont('t''Times New Roman''h15w07')                        !
      write(texto,'(a)') 'kg/m3'                                        !
      call moveto(ix-95,iy+45,xy)                                       !
      call outgtext(texto)                                              !

      ix=145                                                            !
      iy=870                                                            !
      nada=setcolorrgb(#000000)                                         !
      write(texto,'(i6,a)') nint(tic),' m'                              !
      call moveto(ix,iy-2,xy)                                           !
      call outgtext(texto)                                              !
      v1=iy+28                                                          !
      call moveto(ix+3,v1,xy)                                           !
      u1=ix+47                                                          !
      nada=lineto(u1,v1)                                                !
      u1=ix+10                                                          !
      call moveto(u1,v1-5,xy)                                           !
      nada=lineto(u1,v1)                                                !
      u1=ix+10+tic/fe/fdib                                              !
      call moveto(u1,v1-5,xy)                                           !
      nada=lineto(u1,v1)                                                !
      u1=ix-6                                                           !
      call moveto(u1,iy-17,xy)                                          !
      v1=iy+24                                                          !
      nada=lineto(u1,v1)                                                !
      v1=iy+20                                                          !
      call moveto(u1+5,v1,xy)                                           !
      nada=lineto(u1,v1)                                                !
      v1=iy+20-tic/fe/fdib                                              !
      call moveto(u1+5,v1,xy)                                           !
      nada=lineto(u1,v1)                                                !
      nada=setfont('t''Times New Roman''h16w07')                        !
      jx=232                                                            !
      jy=330                                                            !
      call ventana(jx-lx,jy-ly,jx+lx,jy+ly)                             !
      nada=setcolorrgb(#d0f0d0)                                         !
      u1=jx-lx                                                          !
      u2=jx+lx                                                          !
      v1=jy-ly                                                          !
      v2=jy+ly                                                          !
      nada=rectangle($gfillinterior,u1,v1,u2,v2)                        !
      nada=setcolorrgb(#000000)                                         !
      write(texto,'(a)') 'Plan view and profiles'                       !
      call moveto(jx-70,jy-ly-26,xy)                                    !
      call outgtext(texto)                                              !
      call hor0(jx,jy,lx,ly,tic,fe,fdib,mb,nb,xb,yb,ib)                 !
      do 75 ix=1,nc                                                     !
      i=nc-ix+1                                                         !
      if(df(i).eq.0) go to 75                                           !
      if(isen(i).gt.4.5*seli) go to 75                                  !
      if(isen(i).gt.     0) r=2.                                        !
      if(isen(i).gt.1.5*seli) r=3.                                      !
      if(isen(i).gt.2.0*seli) r=4.                                      !
      if(isen(i).gt.2.5*seli) r=5.                                      !
      if(isen(i).gt.3.0*seli) r=6                                       !
      if(isen(i).gt.3.5*seli) r=7                                       !
      if(isen(i).gt.4.0*seli) r=8                                       !
      j=(xc(i)+tx(i)/r)/fe/fdib                                         !
      if(j.gt.lx) go to 75                                              !
      u2=jx+j                                                           !
      j=(xc(i)-tx(i)/r)/fe/fdib                                         !
      if(j.lt.-lx) go to 75                                             !
      u1=jx+j                                                           !
      j=(yc(i)+ty(i)/r)/fe/fdib                                         !
      if(j.gt.ly) go to 75                                              !
      v2=jy-j                                                           !
      j=(yc(i)-ty(i)/r)/fe/fdib                                         !
      if(j.lt.-ly) go to 75                                             !
      v1=jy-j                                                           !
      k=(df(i)+zi)/pa +1                                                !
      nada=setcolorrgb(col(k))                                          !
      nada=rectangle($gfillinterior,u1,v1,u2,v2)                        !
   75 continue                                                          !
      call cont(mb,nb,xb,yb,ib,jx,jy)                                   !
      jx1=jx                                                            !
      jy1=jy                                                            !
      nada=setfont('t''Times New Roman''h14w06')                        !
      nada=setcolorrgb(#000000)                                         !
      do 74 i=-10,10                                                    !
      k=i*tic/fe/fdib                                                   !
      if(abs(k).gt.lx) go to 74                                         !
      v1=jy+k-6                                                         !
      write(texto,'(i7)') nint(iym-i*tic)                               !
      call moveto(jx+lx+6,v1,xy)                                        !
      nada=setfont('t''Times New Roman''h14w06')                        !
      call outgtext(texto)                                              !
      u1=jx+k-6                                                         !
      write(texto,'(i6)') nint((i-1)*tic+ixm)                           !
      call moveto(u1,jy+ly+45,xy)                                       !
      call setgtextrotation(900)                                        !
      nada=setfont('t''Times New Roman''h12w05')                        !
      call outgtext(texto)                                              !
      call setgtextrotation(000)                                        !
   74 continue                                                          !
      nada=setfont('t''Times New Roman''h16w07')                        !
      jy=630                                                            !
      call ventana(jx-lx,jy,jx+lx,jy+lz)                                !
      nada=setcolorrgb(#d0f0d0)                                         !
      u1=jx-lx                                                          !
      u2=jx+lx                                                          !
      v1=jy+30                                                          !
      v2=jy+lz                                                          !
      nada=rectangle($gfillinterior,u1,v1,u2,v2)                        !
      nada=setcolorrgb(#000000)                                         !
      write(texto,'(a)') 'WE elevation view and profiles'               !
      call moveto(jx1-100,jy-25,xy)                                     !
      call outgtext(texto)                                              !
      do 76 i=1,nc                                                      !
      if(df(i).eq.0) go to 76                                           !
      j=zc(i)/fe/fdib                                                   !
      if(j.lt.-lz) go to 76                                             !
      c=ty(i)/2.                                                        !
      if(isen(i).gt.4.5*seli) go to 76                                  !
      if(isen(i).gt.     0) r=2.                                        !
      if(isen(i).gt.1.5*seli) r=3.                                      !
      if(isen(i).gt.2.0*seli) r=4.                                      !
      if(isen(i).gt.2.5*seli) r=5.                                      !
      if(isen(i).gt.3.0*seli) r=6                                       !
      if(isen(i).gt.3.5*seli) r=7                                       !
      if(isen(i).gt.4.0*seli) r=8                                       !
      j=(xc(i)+tx(i)/r)/fe/fdib                                         !
      if(j.gt.lx) go to 76                                              !
      u2=jx+j                                                           !
      j=(xc(i)-tx(i)/r)/fe/fdib                                         !
      if(j.lt.-lx) go to 76                                             !
      u1=jx+j                                                           !
      j=(zc(i)+tz(i)/r)/fe/fdib                                         !
      if(j.gt.lz) go to 76                                              !
      v2=jy-j   +30                                                     !
      j=(zc(i)-tz(i)/r)/fe/fdib                                         !
      if(j.lt.-lz) go to 76                                             !
      v1=jy-j  +30                                                      !
      k=(df(i)+zi)/pa                                                   !
      nada=setcolorrgb(col(k))                                          !
      nada=rectangle($gfillinterior,u1,v1,u2,v2)                        !
   76 continue                                                          !
      call ver0(jx1,jy,lx,ly,lz,izm,tic,fe,fdib)                        !
      jy2=jy                                                            !
      nada=setfont('t''Times New Roman''h14w06')                        !
      nada=setcolorrgb(#000000)                                         !
      j=izm/tic+2                                                       !
      do 77 i=2,20                                                      !
      k=((i-j)*tic+izm)/fe/fdib                                         !
      v1=jy+23+k                                                        !
      if(v1.gt.(jy+lz)) go to 77                                        !
      write(texto,'(i6)') -nint((i-j)*tic)                              !
      call moveto(jx+lx+6,v1,xy)                                        !
      call outgtext(texto)                                              !
   77 continue                                                          !
      nada=setfont('t''Times New Roman''h16w07')                        !
      l=3000*fe                                                         !
      lk=0                                                              !
      do 79 i=1,6                                                       !
      k=-(i+1)*(i+1)*l/17.                                              !
      jy=250                                                            ! 
      if(i.ge.4) jy=700                                                 !
      if(i.eq.1.or.i.eq.4) jx=700                                       !
      if(i.eq.2.or.i.eq.5) jx=1110                                      !
      if(i.eq.3.or.i.eq.6) jx=1520                                      !
      lk=lk+1                                                           !
      call sect(jx,jy,mc,nc,mb,nb,xc,yc,zc,tx,ty,tz,df,isen,pa,zi,      !
     &  xb,yb,ib,col,fe,iym,izm,fdib,tic,dms,seli,jx1,jy2,1,k,le,lk)    !
      call cont(mb,nb,xb,yb,ib,jx,jy)                                   !
   79 continue                                                          !

      write(texto,'(a)') ' ===> Press (Enter)'                          !
      call moveto(1200,950,xy)                                          !
      nada=setcolorrgb(#000000)                                         !
      call outgtext(texto)                                              !
      read(*,'(i1)') i                                                  !
      nada=setcolorrgb(#f0f0f0)                                         !
      nada=rectangle($gfillinterior,500,10,1860,970)                    !
      lk=0                                                              !
      do 78 i=1,15                                                      !
      if(i.ge. 1) jx=700                                                !
      if(i.ge. 6) jx=1110                                               !
      if(i.ge.11) jx=1520                                               !
      if(i.ge. 1) jy=50+(i- 1)*180                                      !
      if(i.ge. 6) jy=50+(i- 6)*180                                      !
      if(i.ge.11) jy=50+(i-11)*180                                      !
      k=(8-i)*l*0.8                                                     !
      lk=lk+1                                                           !
      call sect(jx,jy,mc,nc,mb,nb,xc,yc,zc,tx,ty,tz,df,isen,pa,zi,      !
     &  xb,yb,ib,col,fe,iym,izm,fdib,tic,dms,seli,jx1,jy1,2,k,le,lk)    !
   78 continue                                                          !
c-----------------------------------------------------------------------!
      iv=iv+1
      if(ns.lt.50.or.iv.eq.2) go to 99
      do 82 j=2,nay
      c=(qm(j-1)-qm(j))/cmax
      if(c.lt.0) go to 82 
      do 81 i=1,nc 
   81 if(lev(j).gt.zc(i)) df(i)=df(i)+fs*c  !******
   82 continue
      go to 70

   99 write(texto,'(a)') ' ===> Press (Enter)'                          !
      call moveto(1200,950,xy)                                          !
      nada=setcolorrgb(#000000)                                         !
      call outgtext(texto)                                              !
      read(*,'(i1)') i                                                  !

      stop
      end
      
c***********************************************************************
c  attrac(+uGal) of a right prism, parallel to axis, upon coord origin *
c       xc,yc: coordinates (meters) of the prism center                *
c       dx,dy: sides (meters) of the prism                             *
c       zb,zt: vertical coordinates of prism bottom and top            *
c       d: density (kg/m3)        av: vertical attraction (uGal) (+z)  *
c***********************************************************************
      subroutine VAP(xc,yc,dx,dy,zb,zt,d,av)
      dimension x(4),y(4),s(4),z(2)
         zc=(zt+zb)/2.                                  ! spherical aprox.
         a=xc*xc+yc*yc+zc*zc
              r=(dx+dy+zb-zt)/2.9
              r=r*r 
              if(a.lt.r) a=r                            
         av=dx*dy*(zb-zt)*zc/sqrt(a*a*a) *d*6.672e-3
         if(zc.ne.0.) return
      a=dx/2.
      b=dy/2.
      x(1)=xc+a
      x(2)=xc-a
      x(3)=x(1)
      x(4)=x(2)
      y(1)=yc-b
      y(2)=y(1)
      y(3)=yc+b
      y(4)=y(3)
      z(1)=zb
      z(2)=zt
      do 1 i=1,4
      s(i)=0.
      if(x(i).eq.0.or.y(i).eq.0.) go to 1
      call PIC(x(i),y(i),z,s(i))
    1 continue
      av=(s(3)+s(2)-s(1)-s(4))*d*6.672e-3
      return
c*********************************************
      end
      subroutine PIC(xa,ya,z,av)
      dimension z(2)
      x= abs(xa)
      y= abs(ya)
      xx=x*x
      yy=y*y
      av=0.
      i=1
      do 1 k=1,2
      i=-i
      zz=z(k)*z(k)
      r= sqrt(xx+yy+zz)
      rx= sqrt(xx+zz)
      ry= sqrt(yy+zz)
      if(z(k).eq.0.) go to 1
      av=av+(z(k)* atan(z(k)*r/x/y)-abs(z(k))*1.570796327)*i
    1 av=av+(x* log((y+r)/rx)+y* log((x+r)/ry))*i
      if(xa.lt.0.) av=-av
      if(ya.lt.0.) av=-av
      return
      end

c***********************************************************************
c     Autocorrelation of n<m data v(x,y) for a distance dr             c
c***********************************************************************
      subroutine cov2(m,n,x,y,v,w,dr,mr,cov,sw)
      dimension x(m),y(m),v(m),w(m),cov(mr),sw(mr)
        
      vm=0.
      suw=0.
      do 1 i=1,n
      suw=suw+w(i)
    1 vm=vm+v(i)*w(i)
      vm=vm/suw
      sd2=0
      do 2 i=1,n
      vi=v(i)-vm
    2 sd2=sd2+vi*vi*w(i)
      sd2=sd2/suw
      
      do 7 k=1,mr
      cov(k)=0.
    7 sw(k)=0. 
      dr2=dr*dr
     
      do 3 i=1,n-1
      do 4 j=i+1,n
      xx=x(i)-x(j)
      yy=y(i)-y(j)
      d=sqrt(xx*xx+yy*yy)
      do 5 k=1,mr
      if(d.ge.dr*(k-1).and.d.le.dr*k) then
        cov(k)=cov(k)+(v(j)-vm)*(v(i)-vm)*w(j)*w(i)
        sw(k)=sw(k)   +w(j)*w(i)
      endif
    5 continue  
    4 continue
    3 continue
      vm=0
      do 6 k=1,mr
      if(sw(k).gt.vm) vm=sw(k)    
    6 if(sw(k).ne.0) cov(k)=cov(k)/sw(k)/sd2
      do 8 k=1,mr
    8 sw(k)=sw(k)/vm
      
      return
      end
      
c***********************************************************************
c   median of absolute values for n<m different numbers x(i)           c
c***********************************************************************
      subroutine dmedian(m,n,x,xa)
      dimension x(m)
      xi=x(1)
      xs=xi
      do 1 i=1,n
      xx=x(i)
      if(xx.lt.xi)  xi=xx
      if(xx.gt.xs)  xs=xx
        do 6 j=i+1,n
    6   if(xx.eq.x(j)) x(j)=x(j)*1.00001
    1 continue
      j=0
    4 xa=xs
      do 2 i=1,n
      xx=x(i)
      if(xx.lt.xa.and.xx.gt.xi) xa=xx
    2 continue
      j=j+1
      if(j.eq.n) go to 9
      xi=xa
      xb=xi
      do 3 i=1,n
      xx=x(i)
      if(xx.gt.xb.and.xx.lt.xs) xb=xx
    3 continue
      j=j+1
      xa=(xa+xb)/2.
      xs=xb
      if(j.lt.n) go to 4
    9 return
      end
      
c***********************************************************************
c     Execution stops when number of data is higher than dimensions    *
c***********************************************************************
      subroutine dim(m)
      write(*,200) m
  200 format(/6x,'*** Error: data size > max. dimension',i6,' !!!'/
     -10x,'Reduce the data side or change dimension in the code'/)
      stop
      end
      
c***********************************************************************
c     Sensitivity calculus
c***********************************************************************
      subroutine at2(ms,ns,x,y,z,xx,yy,zz,dx,dy,dz,sen,av,zp)
      dimension x(ms),y(ms),z(ms)
      av=0.
      t=(dx+dy+dz)/3
      t2=t*t
      zp=0
      sp=0 
      do 1 i=1,ns
      zr=zz-z(i)
      xr=xx-x(i)
      yr=yy-y(i)
      d=xr*xr+yr*yr+zr*zr+500.
      if(d.lt.t2) d=t2
      av=av+zr*zr/d/d/d
      d=1/d/d
      sp=sp+d
    1 zp=zp+z(i)*d
      av=av/ns
      sen=3.d-7/sqrt(av)
      t2=dx*dy*dz
      av=av*t2*t2  
      zp=zp/sp
      return
      end

c***********************************************************************
c     New cells for a layer
c           zs: top of the layer    tz: width of the layer 
c***********************************************************************
      subroutine celdascapa(ax,bx,ay,by,zs,tz,avm,ms,ns,x,y,z,
     -mc,nc,xc,yc,zc,dx,dy,dz,at,seli,ncr,isen)
      dimension xc(mc),yc(mc),zc(mc),dx(mc),dy(mc),dz(mc),
     -x(ms),y(ms),z(ms),isen(mc)
      character*50 texto
      real*4 at(mc)
      zi=zs-tz                           
      zz=zs-tz/2.
      tx=tz                             
      ty1=tz/3                          
      ty2=tz*3                          
      tt=tz*tz 
      kk=nc
      
      do 7 i=1,999                       
      xx=ax+tx*(i-1)
      if(xx.gt.bx) go to 9
      ty=tz
      yr=ay
      do 6 j=1,999                        
      if(yr.ge.by) go to 7
      c=ty
      l=0
    3 yy=yr+c/2
      call at2(ms,ns,x,y,z,xx,yy,zz,tx,c,tz,sen,am,zp)
      if(sen.gt.seli*10.) go to 6             
      if((zp-100).lt.zs.or.am.le.0) go to 6   
      av=am/avm
      if(abs(av-1.).le.0.1) go to 4        
      c=c/av**0.41
      if(c.gt.ty2) go to 6           
      if(c.lt.ty1) go to 6           
      l=l+1
      if(l.le.10) go to 3
      go to 6
    4 if(c.gt.9999) go to 6
    5 ty=c                              

c   Eliminacion celdas superpuestas

      if(ncr.eq.0) go to 12
      do 11 k=1,ncr
      if(abs(xx-xc(k)).gt.(dx(k)+tx)/2.6) go to 11  
      if(abs(yy-yc(k)).gt.(dy(k)+ty)/2.6) go to 11  
      if(abs(zz-zc(k)).gt.(dz(k)+tz)/2.6) go to 11  
      go to 6
   11 continue
   12 kk=kk+1
      xc(kk)=xx
      yc(kk)=yy
      zc(kk)=zz
      dx(kk)=tx
      dy(kk)=ty
      dz(kk)=tz
      at(kk)=av
      isen(kk)=0
      if(ncr.gt.1) isen(kk)=1
      if(kk.lt.mc) go to 6 
         kk=-9999
         go to 9
    6 yr=yr+ty
    7 continue
      
    9 nc=kk
      return
      end
      
c***********************************************************************
c   Width of a layer
c***********************************************************************
      subroutine especapa(zs,avm,ms,ns,x,y,z,tz)
      dimension x(ms),y(ms),z(ms)
      l=0
      c=tz                           
    1 zz=zs-c*0.5                    
      k=0
      ama=-9.d9
      ame=0
      do 2 i=1,ns
      call at2(ms,ns,x,y,z,x(i)*1.5,y(i)*1.5,zz,c,c,c,se,aa,zp)
      if((zp-zz).le.100) go to 2
      if(aa.gt.ama) ama=aa
      ame=ame+aa
      k=k+1
    2 continue
      if(k.le.1) go to 9
      am=ame/k/avm
      c=c/am**0.13
      l=l+1
      if(abs(am-1.).gt.0.1.and.l.lt.12) go to 1
      
      if(c.gt.tz) tz=c
    9 return
      end
      
c***********************************************************************
c     Step for autocorrelation                                         c
c***********************************************************************
      subroutine step(m,n,x,y,d,dr)
      dimension x(m),y(m),d(m)
      do 1 i=1,n
      d4=9.e9
      d3=9.e9
      d2=9.e9
      de=9.e9
      do 2 j=1,n
      xx=x(i)-x(j)
      yy=y(i)-y(j)
      dd= xx*xx+yy*yy
      if(dd.gt.d4.or.dd.le.1.) go to 2
      if(dd.lt.de) de=dd
      if(dd.lt.d2.and.dd.gt.de) d2=dd
      if(dd.lt.d3.and.dd.gt.d2) d3=dd
      if(dd.lt.d4.and.dd.gt.d3) d4=dd
    2 continue
      de=sqrt(de)
      d2=sqrt(d2)
      d3=sqrt(d3)
      d4=sqrt(d4)
    1 d(i)=de*0.4+d2*0.3+d3*0.2+d4*0.1  
      call dmedian(m,n,d,dr)
      return
      end
      
c****************************************************************
c     Square window  in drawing                                 *
c****************************************************************
      subroutine ventana(i1,j1,i2,j2)
      use IFQWIN
      integer*2 x1,x2,y1,y2
      nada=setcolorrgb(#909090)
      x1=i1-2
      y1=j1-2
      x2=i1-1
      y2=j2+1
      nada=rectangle($gfillinterior,x1,y1,x2,y2)
      x2=i2+1
      y2=j1-1
      nada=rectangle($gfillinterior,x1,y1,x2,y2)
      nada=setcolorrgb(#fefefe)
      x1=i1-1
      y1=j2+2 
      x2=i2+2
      y2=j2+2
      nada=rectangle($gfillinterior,x1,y1,x2,y2)
      y1=j1
      nada=rectangle($gfillinterior,x2,y1,x2,y2)
      nada=setcolorrgb(#ffffff)
      x1=i1
      x2=i2
      y2=j2
      nada=rectangle($gfillinterior,x1,y1,x2,y2)
      return
      end
      
c**********************************************************************
c       Rounding numbers                                              *
c**********************************************************************    
      subroutine round(c,r)
      j=1
      if(c.lt.0) j=-1
      r=abs(c)
      i=log10(r)                                               
      r=j*nint(r/10**i)*10**i  
      return 
      end
      
c***********************************************************************
c     Dialog interface 1 (data files)                                  *
c***********************************************************************
      subroutine dial1(obs,mod,fil,na,fla,tpc,dco,iud,mlev,blun)
      use dflogm
      include 'Growth.fd'
      logical retlog,nada
      character*50 texto,obs,mod,fil
      external fin
      type(dialog) dlg
      retlog=DlgInit( IDD_Growth, dlg)

      retlog = DlgSet(dlg,IDC_EDIT_e1 ,"Grad.txt" )
      retlog = DlgSet(dlg,IDC_EDIT_e3 ,"Mod.txt" )
      retlog = DlgSet(dlg,IDC_EDIT_e4 ,"Fil.txt" )
      write(texto,'(f5.1)') dco 
      retlog = DlgSet(dlg,IDC_EDIT_e8,texto)
      write(texto,'(i1)') mlev  
      retlog = DlgSet(dlg,IDC_EDIT_e13,texto)
      write(texto,'(i3)') na            
      retlog = DlgSet(dlg,IDC_EDIT_e5,texto)
      write(texto,'(i3)') nint(fla)      
      retlog = DlgSet(dlg,IDC_EDIT_e6,texto)
      write(texto,'(f4.0)') tpc 
      retlog = DlgSet(dlg,IDC_EDIT_e7,texto)
      write(texto,'(f3.1)') blun 
      retlog = DlgSet(dlg,IDC_EDIT_e9,texto)

      write(texto,'(i3)') iud    
      retlog = DlgSet(dlg,IDC_EDIT_e15,texto) 
      retlog = DlgSet(dlg,IDC_CHECK_c2,.false.)
      retlog = DlgSet(dlg,IDC_CHECK_c1,.true.)      
      retlog=DlgSet(dlg,IDC_TEXT_t2 ,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_TEXT_t5 ,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_TEXT_t6 ,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_TEXT_t7 ,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_TEXT_t8 ,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_TEXT_t9 ,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_TEXT_t10,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_TEXT_t12,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_TEXT_t13,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_TEXT_t14,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_TEXT_t18,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_EDIT_e2 ,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_EDIT_e5 ,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_EDIT_e6 ,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_EDIT_e7 ,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_EDIT_e15,.false.,dlg_enable) 
      retlog=DlgSet(dlg,IDC_EDIT_e8,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_EDIT_e9,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_EDIT_e13,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_BOX_b3  ,.false.,dlg_enable) 
      retlog=DlgSet(dlg,IDC_BOX_b4  ,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_CHECK_c2,.false.,dlg_enable)      
      retlog=DlgSet(dlg,IDC_CHECK_c1,.false.,dlg_enable)
      retlog = DlgSetSub(dlg,IDCANCEL,fin)

      retint = DlgModal( dlg )

      nada=DlgGet(dlg,IDC_EDIT_e1, obs )
      nada=DlgGet(dlg,IDC_EDIT_e3, mod )
      nada=DlgGet(dlg,IDC_EDIT_e4, fil )

      return
      end      
           
c***********************************************************************
c     Dialog interface 2  (average cell size)                          *
c***********************************************************************
      subroutine dial2(obs,mod,fil,na,ns,fla,tpc,dco,iud,mlev,pp,blun)
      use dflogm
      include 'Growth.fd'
      logical retlog,nada
      character*50 texto,obs,mod,fil
      external fin
      type(dialog) dlg
      retlog=DlgInit( IDD_Growth, dlg)
                                                    
      retlog = DlgSet(dlg,IDC_EDIT_e1,obs )
      retlog = DlgSet(dlg,IDC_EDIT_e3,mod )
      retlog = DlgSet(dlg,IDC_EDIT_e4,fil )

      retlog = DlgSet(dlg,IDC_TEXT_t1,.false.,dlg_enable)
      retlog = DlgSet(dlg,IDC_TEXT_t3,.false.,dlg_enable)
      retlog = DlgSet(dlg,IDC_TEXT_t4,.false.,dlg_enable)
      retlog = DlgSet(dlg,IDC_TEXT_t9 ,.false.,dlg_enable)
      retlog = DlgSet(dlg,IDC_TEXT_t10 ,.false.,dlg_enable)
      retlog = DlgSet(dlg,IDC_BOX_b1 ,.false.,dlg_enable)
      retlog = DlgSet(dlg,IDC_EDIT_e1,.false.,dlg_enable)
      retlog = DlgSet(dlg,IDC_EDIT_e3,.false.,dlg_enable)
      retlog = DlgSet(dlg,IDC_EDIT_e4,.false.,dlg_enable)
      retlog = DlgSet(dlg,IDC_EDIT_e7,.false.,dlg_enable)
      retlog = DlgSet(dlg,IDC_EDIT_e9,.false.,dlg_enable)
      retlog = DlgSet(dlg,IDC_CHECK_c1,.false.,dlg_enable)
      retlog = DlgSet(dlg,IDC_CHECK_c2,.false.,dlg_enable)
      
      write(texto,'(i6)') ns
      retlog = DlgSet(dlg,IDC_EDIT_e30,texto)
      write(texto,'(f5.1)') dco     
      retlog = DlgSet(dlg,IDC_EDIT_e8,texto)
      write(texto,'(i1)') mlev
      retlog = DlgSet(dlg,IDC_EDIT_e13,texto)
      write(texto,'(i3)') na            
      retlog = DlgSet(dlg,IDC_EDIT_e5,texto)
      write(texto,'(i3)') nint(fla)      
      retlog = DlgSet(dlg,IDC_EDIT_e6,texto)
      write(texto,'(f4.0)') tpc     
      retlog = DlgSet(dlg,IDC_EDIT_e7,texto) 
      write(texto,'(i3)') iud     
      retlog = DlgSet(dlg,IDC_EDIT_e15,texto) 
      retlog = DlgSet(dlg,IDC_CHECK_c2,.false.)
      retlog = DlgSet(dlg,IDC_CHECK_c1,.true.)   
      write(texto,'(f3.1)') blun        
      retlog = DlgSet(dlg,IDC_EDIT_e9,texto) 
      write(texto,'(i2)') mlev  
      retlog = DlgSet(dlg,IDC_EDIT_e13,texto)

      retlog=DlgSet(dlg,IDC_TEXT_t2 ,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_TEXT_t5 ,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_TEXT_t6 ,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_TEXT_t7 ,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_TEXT_t8 ,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_TEXT_t8 ,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_TEXT_t13,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_TEXT_t14,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_EDIT_e5 ,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_EDIT_e6 ,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_EDIT_e8,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_EDIT_e15,.false.,dlg_enable)      
      retlog=DlgSet(dlg,IDC_EDIT_e12,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_EDIT_e13,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_BOX_b4  ,.false.,dlg_enable)
     
      write(texto,'(i6)') nint(pp)
      nada=DlgSet(dlg,IDC_EDIT_e2,texto)
      
      retlog = DlgSetSub(dlg,IDCANCEL,fin)

      retint = DlgModal( dlg )

      retlog = DlgGet( dlg, IDC_EDIT_e2,texto)
      read(texto,*) pp
      
      return
      end
      
      
c***********************************************************************
c     Dialog interface 3  (inversion parameters)                       *
c***********************************************************************
      subroutine dial3(obs,mod,fil,na,nc,ns,fla,dco,mlev,
     & pp,ig0,igxy,tpc,blun,iud)
      use dflogm
      include 'Growth.fd'
      logical retlog,nada
      character*50 texto,obs,mod,fil
      external fin
      type(dialog) dlg
      retlog=DlgInit( IDD_Growth, dlg)
                                                    
      retlog = DlgSet(dlg,IDC_EDIT_e1,obs )
      retlog = DlgSet(dlg,IDC_EDIT_e3,mod )
      retlog = DlgSet(dlg,IDC_EDIT_e4,fil )

      retlog = DlgSet(dlg,IDC_TEXT_t1,.false.,dlg_enable)
      retlog = DlgSet(dlg,IDC_TEXT_t3,.false.,dlg_enable)
      retlog = DlgSet(dlg,IDC_TEXT_t4,.false.,dlg_enable)
      retlog = DlgSet(dlg,IDC_TEXT_t18,.false.,dlg_enable)
      retlog = DlgSet(dlg,IDC_BOX_b1 ,.false.,dlg_enable)
      retlog = DlgSet(dlg,IDC_EDIT_e1,.false.,dlg_enable)
      retlog = DlgSet(dlg,IDC_EDIT_e2,.false.,dlg_enable)      
      retlog = DlgSet(dlg,IDC_EDIT_e3,.false.,dlg_enable)
      retlog = DlgSet(dlg,IDC_EDIT_e4,.false.,dlg_enable)
      write(texto,'(i6)') ns
      retlog = DlgSet(dlg,IDC_EDIT_e30,texto)
      write(texto,'(i6)') nc
      retlog = DlgSet(dlg,IDC_EDIT_e31,texto)
      write(texto,'(i3)') na            
      retlog = DlgSet(dlg,IDC_EDIT_e5,texto)
      write(texto,'(i3)') nint(fla)      
      retlog = DlgSet(dlg,IDC_EDIT_e6,texto)
      write(texto,'(f4.0)') tpc        
      retlog = DlgSet(dlg,IDC_EDIT_e7,texto)
      write(texto,'(f3.1)') uw     
      retlog = DlgSet(dlg,IDC_EDIT_e15,texto)
      write(texto,'(f5.1)') dco     
      retlog = DlgSet(dlg,IDC_EDIT_e8,texto)
      write(texto,'(i1)') mlev  
      retlog = DlgSet(dlg,IDC_EDIT_e13,texto)
      write(texto,'(f3.1)') blun 
      retlog = DlgSet(dlg,IDC_EDIT_e9,texto)
      if(ig0.eq.0) retlog = DlgSet(dlg,IDC_CHECK_c1,.false.)
      if(ig0.eq.1) retlog = DlgSet(dlg,IDC_CHECK_c1,.true.)
      if(igxy.eq.0) retlog = DlgSet(dlg,IDC_CHECK_c2,.false.)
      if(igxy.eq.1) retlog = DlgSet(dlg,IDC_CHECK_c2,.true.)
      write(texto,'(i6)') nint(pp)
      nada=DlgSet(dlg,IDC_EDIT_e2,texto)
      retlog = DlgSetSub(dlg,IDCANCEL,fin)
      
      retint = DlgModal( dlg )
                                                  
      retlog = DlgGet( dlg, IDC_EDIT_e5,texto)
      read (texto,*) na
      if(na.lt.1) na=1
      retlog = DlgGet( dlg, IDC_EDIT_e6, texto)
      read(texto,*) fla
      retlog = DlgGet( dlg, IDC_EDIT_e7, texto)
      read(texto,*) tpc
      retlog = DlgGet( dlg, IDC_EDIT_e8, texto)
      read(texto,*) dco
      retlog = DlgGet( dlg, IDC_EDIT_e9, texto)
      read(texto,*) blun
      retlog = DlgGet( dlg, IDC_EDIT_e13, texto)
      read(texto,*) mlev
      retlog = DlgGet( dlg, IDC_EDIT_e15, texto)
      read(texto,*) iud 
      ig0=0
      retlog = DlgGet( dlg, IDC_CHECK_c1, nada)
      if(nada) ig0=1
      igxy=0
      retlog = DlgGet( dlg, IDC_CHECK_c2, nada)
      if(nada) igxy=1

      return
      end
      
c***********************************************************************
c     Simple sensitivity calculus  Ton/uGal
c***********************************************************************
      subroutine ss2(ms,ns,x,y,z,fe,xx,yy,zz,sen,k)
      dimension x(ms),y(ms),z(ms)
      av=0
      r=100*fe*fe
      do 1 j=1,ns                          
      xd=xx-x(j)
      yd=yy-y(j)
      zd=zz-z(j)                            
      d=xd*xd+yd*yd+zd*zd + r
      c=zd*zd/d/d/d
    1 av=av+c
      av=av/ns
      av=6.672e-3 *sqrt(av)*1.e8   
      sen=1./av                    
      k=sen/fe/fe/50.

      if(k.lt.1) k=1
      if(k.gt.12) k=12

      return
      end

c***********************************************************************   
c     Histogram for residuals,  n columns
c***********************************************************************
      subroutine histog(ms,ns,w,re,kdib,ix,lx,iy,ly)
      use IFQWIN
      use IFPORT
      type (xycoord) xy
      integer*2 u1,u2,v1,v2
      character*140 texto
      dimension re(ms),w(ms),k(20)
      sw=0
      do 5 i=1,ns
    5 sw=sw+w(i)
      sw=sw/ns
      rmin=0
      rmax=0
      rmed=0
      do 1 i=1,ns
      c=re(i)  *w(i)/sw
      if(c.gt.rmax) rmax=c    
      if(c.lt.rmin) rmin=c   
      rmed=rmed+c
    1 continue
      rmed=rmed/ns   
      n=ns/4.7 
      if(n.gt.20) n=20
      p=(rmax-rmin)/n
      py=0
      xce=0
      do 2 j=1,n
      k(j)=0
      r1=rmin+(j-1)*p
      r2=r1+p
      do 3 i=1,ns
      c=re(i) *w(i)/sw
    3 if(c.gt.r1.and.c.le.r2) k(j)=k(j)+1   
      if(k(j).gt.py) py=k(j)
    2 continue

      nada=setcolorrgb(#000000)                                         !
      write(texto,'(i9,8x,i9)') nint(rmin),nint(rmax)                  !
      call moveto(ix-44,iy+5,xy)                                        !
      call outgtext(texto)                                              !
      write(texto,'(i4)') 0                                             !
      call moveto(ix-27,iy-10,xy)                                       !
      call outgtext(texto)                                              !
      write(texto,'(i4)') nint(py/2)                                    !
      j=iy-ly/1.8                                                       !
      call moveto(ix-27,j,xy)                                           !
      call outgtext(texto)                                              !
      write(texto,'(i4)') nint(py)                                      !
      call moveto(ix-27,iy-ly,xy)                                       !
      call outgtext(texto)                                              !
      if(kdib.eq.0) then                                                !
       write(texto,*) 'Resid.(uGal)'                                    !
       call moveto(ix+10,iy+35,xy)                                      !
       call outgtext(texto)                                             !
      endif                                                             !
      py=ly/py*0.9                                                      ! 
      px=lx/n                                                           !
      do 4 j=1,n                                                        !
      l=py*k(j)                                                         !
      u1=ix+5+(j-1)*p/(rmax-rmin)*lx +1                                 !
      u2=ix+5+j*p   /(rmax-rmin)*lx  -2                                 !
      v1=iy                                                             !
      v2=iy-l                                                           !
      nada=setcolorrgb(#60d060)                                         !
      nada=rectangle($gfillinterior,u1,v1,u2,v2)                        !
    4 continue                                                          !                  
      nada=setcolorrgb(#000000)                                         !
      call moveto(ix+5,iy,xy)                                           !
      u1=ix+5                                                           !
      v1=iy-ly                                                          !
      nada=lineto(u1,v1)                                                !
      i=(rmed-rmin)/(rmax-rmin)*lx                                      !
      call moveto(ix+5+i,iy,xy)                                         !
      u1=ix+5+i                                                         !
      nada=lineto(u1,v1)                                                !
      write(texto,'(i4)') nint(rmed)                                    !
      call moveto(u1-20,iy+5,xy)                                        !
      call outgtext(texto)                                              !
      call moveto(ix+5,iy,xy)                                           !
      u1=ix+5+lx                                                        !
      v1=iy                                                             !
      nada=lineto(u1,v1)                                                !
   99 return
      end

c***********************************************************************   
c     Draw a section of the 3D anomalous density model
c***********************************************************************
      subroutine sect(jx,jy,mc,nc,mb,nb,xc,yc,zc,tx,ty,tz,df,isen,pa,
     &  zi,xb,yb,ib,col,fe,iym,izm,fdib,tic,dms,seli,jx1,jy1,itipo,
     & ival,le,lk)
      use IFQWIN
      use IFPORT
      type (xycoord) xy
      dimension 
     & xc(mc),yc(mc),zc(mc),tx(mc),ty(mc),tz(mc),df(mc),isen(mc),
     & xb(mb),yb(mb),ib(mb)
      integer col(21),nada
      integer*2 u1,u2,v1,v2
      character*120 texto
      character*1 le(15)

      lx=190                                                            !
      ly=190                                                            !
      lz=140                                                            !
      if(itipo.eq.1) then                                               !
      call ventana(jx-lx,jy-ly,jx+lx,jy+ly)                             !
      nada=setcolorrgb(#d0f0d0)                                         !
      u1=jx-lx                                                          !
      u2=jx+lx                                                          !
      v1=jy-ly                                                          !
      v2=jy+ly                                                          !
      nada=rectangle($gfillinterior,u1,v1,u2,v2)                        !
      call hor0(jx,jy,lx,ly,tic,fe,fdib,mb,nb,xb,yb,ib)                 !
      do 3 i=1,nc                                                       !
      c=tz(i)/2.                                                        !
      if(isen(i).gt.4.5*seli) go to 3                                   !
      if(isen(i).gt.     0) r=2.                                        !                                          
      if(isen(i).gt.1.5*seli) r=3.                                      !
      if(isen(i).gt.2.0*seli) r=4.                                      !
      if(isen(i).gt.2.5*seli) r=5.                                      !
      if(isen(i).gt.3.0*seli) r=6                                       !
      if(isen(i).gt.3.5*seli) r=7                                       !
      if(isen(i).gt.4.0*seli) r=8                                       !
      if((zc(i)+c).lt.ival) go to 3                                     !
      if((zc(i)-c).ge.ival) go to 3                                     !
      j=(xc(i)+tx(i)/r)/fe/fdib                                         !
      if(j.gt.lx) go to 3                                               !
      u2=jx+j                                                           !
      j=(xc(i)-tx(i)/r)/fe/fdib                                         !
      if(j.lt.-lx) go to 3                                              !
      u1=jx+j                                                           !
      j=(yc(i)+ty(i)/r)/fe/fdib                                         !
      if(j.gt.ly) go to 3                                               !
      v2=jy-j                                                           !
      j=(yc(i)-ty(i)/r)/fe/fdib                                         !
      if(j.lt.-ly) go to 3                                              !
      v1=jy-j                                                           !
      k=(df(i)+zi)/pa+1                                                 !
      nada=setcolorrgb(col(k))                                          !
      nada=rectangle($gfillinterior,u1,v1,u2,v2)                        !
    3 continue                                                          !
      write(texto,'(a1,8x,a,i9)') le(lk),'Depth(m) bsl =',ival+izm      ! 
      call moveto(jx-90,jy-ly-20,xy)                                    !
      nada=setcolorrgb(#000000)                                         !
      call outgtext(texto)                                              !
      k=ival/fe/fdib                                                    !
      write(texto,'(a1)') le(lk)                                        ! 
      call moveto(jx1-lx-15,jy1+22-k,xy)                                ! 
      call outgtext(texto)                                              !
      nada=setcolorrgb(#c000c0)                                         !
      ix=jx-lx                                                          !
      call moveto(jx1-lx,jy1+30-k,xy)                                   !
      u1=jx1+lx                                                         !
      v1=jy1+30-k                                                       !
      if(k.lt.2*lx) nada=lineto(u1,v1)                                  !
      endif                                                             !

      if(itipo.eq.2) then                                               !
      call ventana(jx-lx,jy,jx+lx,jy+lz)                                !
      nada=setcolorrgb(#d0f0d0)                                         !
      u1=jx-lx                                                          !
      u2=jx+lx                                                          !        
      v1=jy+30                                                          !
      v2=jy+lz                                                          !        
      nada=rectangle($gfillinterior,u1,v1,u2,v2)                        !
      do 6 i=1,nc                                                       !
      j=zc(i)/fe/fdib                                                   ! 
      if(j.lt.-lz+30) go to 6                                           !
      c=ty(i)/2.                                                        !
      if(isen(i).gt.4.5*seli) go to 6                                   !
      if(isen(i).gt.     0) r=2.                                        !                                          
      if(isen(i).gt.1.5*seli) r=3.                                      !
      if(isen(i).gt.2.0*seli) r=4.                                      !
      if(isen(i).gt.2.6*seli) r=5.                                      !
      if(isen(i).gt.3.0*seli) r=6                                       !
      if(isen(i).gt.3.5*seli) r=7                                       !
      if(isen(i).gt.4.0*seli) r=8                                       !
      if((yc(i)+c).lt.ival) go to 6                                     !
      if((yc(i)-c).ge.ival) go to 6                                     !
      j=(xc(i)+tx(i)/r)/fe/fdib                                         !
      if(j.gt.lx) go to 6                                               !
      u2=jx+j                                                           !
      j=(xc(i)-tx(i)/r)/fe/fdib                                         !
      if(j.lt.-lx) go to 6                                              !
      u1=jx+j                                                           !
      j=(zc(i)+tz(i)/r)/fe/fdib                                         !
      if(j.gt.lz) go to 6                                               !
      v2=jy-j   +30                                                     !
      j=(zc(i)-tz(i)/r)/fe/fdib                                         !
      if(j.lt.-lz) go to 6                                              !
      v1=jy-j  +30                                                      !
      k=(df(i)+zi)/pa+1                                                 !
      nada=setcolorrgb(col(k))                                          !
      nada=rectangle($gfillinterior,u1,v1,u2,v2)                        !
    6 continue                                                          !
      call ver0(jx,jy,lx,ly,lz,izm,tic,fe,fdib)                         !
      write(texto,'(a1,8x,a,i9)') le(lk),'Y UTM (m) =',ival+iym         !
      call moveto(jx-90,jy-20,xy)                                       !
      nada=setcolorrgb(#000000)                                         !
      call outgtext(texto)                                              !
      k=ival/fe/fdib                                                    !
      write(texto,'(a1)') le(lk)                                        ! 
      call moveto(jx1-lx-15,jy1-k-10,xy)                                ! 
      call outgtext(texto)                                              !
      nada=setcolorrgb(#c000c0)                                         ! 
      ix=jx-lx                                                          !
      call moveto(jx1-lx,jy1-k,xy)                                      !
      u1=jx1+lx                                                         ! 
      v1=jy1-k                                                          !
      if(k.lt.2*lx) nada=lineto(u1,v1)                                  !
      endif                                                             !

   99 return
      end

c************************************************************************
c     Draw a background for a plan view of the model
c************************************************************************
      subroutine hor0(jx,jy,lx,ly,tic,fe,fdib,mb,nb,xb,yb,ib)
      use IFQWIN
      use IFPORT
      type (xycoord) xy
      dimension xb(mb),yb(mb),ib(mb)
      integer*2 u1,u2,v1,v2
      nada=setcolorrgb(#d0f0d0)                                         !
      u1=jx-lx                                                          !
      u2=jx+lx                                                          !
      v1=jy-ly                                                          !
      v2=jy+ly                                                          !
      nada=rectangle($gfillinterior,u1,v1 ,u2,v2)                       !
      nada=setcolorrgb(#000000)                                         !
      do 1 i=-10,10                                                     !
      ix=jx-lx                                                          !
      k=i*tic/fe/fdib                                                   !
      v1=jy+k                                                           !
      call moveto(ix,v1,xy)                                             !
      u1=ix+4                                                           !
      if(abs(k).lt.ly) nada=lineto(u1,v1)                               !
      ix=jx+lx                                                          !
      call moveto(ix,v1,xy)                                             !
      u1=ix-4                                                           !
      if(abs(k).lt.ly) nada=lineto(u1,v1)                               !
      ix=jx+k                                                           !
      iy=jy-ly                                                          !
      call moveto(ix,iy,xy)                                             !
      u1=ix                                                             !
      v1=iy+4                                                           !
      if(abs(k).lt.lx) nada=lineto(u1,v1)                               !
      iy=jy+ly                                                          !
      call moveto(ix,iy,xy)                                             !
      v1=iy-4                                                           !
      if(abs(k).lt.lx) nada=lineto(u1,v1)                               !
    1 continue
      call cont(mb,nb,xb,yb,ib,jx,jy)                                   !
      return
      end

c***********************************************************************
c     Draw a background for a elevation view of the model
c***********************************************************************
      subroutine ver0(jx,jy,lx,ly,lz,izm,tic,fe,fdib)
      use IFQWIN
      use IFPORT
      type (xycoord) xy
      integer*2 u1,u2,v1,v2

      v1=jy+30+izm/fe/fdib                                         
      nada=setcolorrgb(#f0c000)                                        
      call moveto(jx-lx,v1,xy)  
      u1=jx+lx                                       
      nada=lineto(u1,v1) 

      nada=setcolorrgb(#000000)  
      j=izm/tic+2                                                      
      do 1 i=-10,10                                               
      k=((i-j)*tic+izm)/fe/fdib
      v1=jy+30+k            
      if(v1.lt.(jy+lz).and.v1.gt.jy+10) then
        ix=jx-lx            
        call moveto(ix,v1,xy) 
        u1=ix+4               
        nada=lineto(u1,v1)    
        ix=jx+lx              
        call moveto(ix,v1,xy) 
        u1=ix-4               
        nada=lineto(u1,v1)                              
      endif

      k=i*tic/fe/fdib
      ix=jx+k         
      call moveto(ix,jy,xy) 
      u1=ix                 
      v1=jy+5               
      if(abs(k).lt.lx) nada=lineto(u1,v1)   
    1 continue

      return
      end

c*********************************************************************
c       End
c*********************************************************************
      subroutine fin
      use dialogm
      include 'Growth.fd'
      stop
      end 

c*********************************************************************
c   Iterative reweighting of the data according to previous residuals
c*********************************************************************
      subroutine rewe(m,n,re,blun,w)
      dimension re(m),w(m)  
      rme=0
      do 1 i=1,n
      ra=abs(re(i))
      rme=rme+ra
    1 w(i)=ra
      rme=rme/n/0.6745 *blun
      c=4  !5
      sw=0
      do 2 i=1,n
      t=abs(re(i))/rme-1
      w(i)=1./(1+exp(c*t))
      if(w(i).lt.0.01) w(i)=0.01  
    2 sw=sw+w(i)
      sw=sw/n
      do 3 i=1,n
    3 w(i)=w(i)/sw
      return
      end

c************************************************************************
c     Contour line xb,yb
c************************************************************************
      subroutine cont(mb,nb,xb,yb,ib,jx,jy)
      use IFQWIN
      use IFPORT
      type (xycoord) xy
      dimension xb(mb),yb(mb),ib(mb)
      integer*2 u1,v1
      nada=setcolorrgb(#000000)                                         !
      do 2 i=2,nb                                                       !
      if(ib(i-1).eq.1) go to 2                                          !
      if(ib(i).eq.1.or.ib(i).eq.2) go to 2                              !
      ix=xb(i-1)+jx                                                     !
      iy=jy-yb(i-1)                                                     !
      call moveto(ix,iy,xy)                                             !
      u1=xb(i)+jx                                                       !
      v1=jy-yb(i)                                                       !
      nada=lineto(u1,v1)                                                !
    2 continue
      return
      end

c*********************************************************************
c   Up and down cells 
c*********************************************************************
      subroutine updo(may,nay,ice,mc,nc,xc,yc,zc,tx,ty,tz,iup,ido)
      dimension xc(mc),yc(mc),zc(mc),tx(mc),ty(mc),tz(mc),
     - iup(mc),ido(mc),ice(may)

      do 1 k=1,nay
      n1=1
      if(k.gt.1) n1=ice(k-1)+1 
      do 2 i=n1,ice(k)
      iup(i)=0
      ido(i)=0
      tm=tx(i)*ty(i)
      if(k.eq.1) go to 5 
      n=1
      if(k.gt.2) n=ice(k-2)+1
      dm=tm
      do 3 j=n,ice(k-1)
      dx=xc(j)-xc(i)
      dy=yc(j)-yc(i)
      d=dx*dx+dy*dy
      if(d.ge.dm) go to 3
      dm=d
      iup(i)=j 
    3 continue
    5 continue

      if(k.eq.nay) go to 2 
      dm=tm
      do 4 j=ice(k)+1,ice(k+1)
      dx=xc(j)-xc(i)
      dy=yc(j)-yc(i)
      d=dx*dx+dy*dy
      if(d.ge.dm) go to 4
      dm=d
      ido(i)=j 
    4 continue
    2 continue
    1 continue
      return
      end
c************************************************************************