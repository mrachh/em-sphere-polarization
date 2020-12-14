      implicit real *8 (a-h,o-z)
      complex *16, allocatable :: cmevalsin(:,:,:),cnevalsin(:,:,:)
      complex *16, allocatable :: cmovalsin(:,:,:),cnovalsin(:,:,:)
      complex *16, allocatable :: cmevalsout(:,:,:),cnevalsout(:,:,:)
      complex *16, allocatable :: cmovalsout(:,:,:),cnovalsout(:,:,:)
      complex *16, allocatable :: cmevalsinc(:,:,:),cnevalsinc(:,:,:)
      complex *16, allocatable :: cmovalsinc(:,:,:),cnovalsinc(:,:,:)

      complex *16, allocatable :: cevalsin(:,:),chvalsin(:,:)
      complex *16, allocatable :: cevalsout(:,:),chvalsout(:,:)
      complex *16, allocatable :: cevalsinc(:,:),chvalsinc(:,:)
      complex *16, allocatable :: acoefs(:),bcoefs(:),ccoefs(:),
     1   dcoefs(:)
      complex *16 zk,zkin,zkout,ima
      real *8, allocatable :: qin(:,:),win(:)
      real *8, allocatable :: qout(:,:),wout(:)
      complex *16 etest(3),zfac
      data ima/(0.0d0,1.0d0)/


      call prini(6,13)
      
      done = 1.0d0
      pi = atan(done)*4.0d0
      rlam = 0.5d0
      zk = rlam*2*pi/2



      rfac = 1.33
      zkout = zk
      zkin = zkout*rfac

      kthet = 30
      kphi = 20
      kr = 20
      nin = kthet*kphi*kr
      allocate(qin(3,nin),win(nin))

      nout = kthet*kphi*kr
      allocate(qout(3,nout),wout(nout))

      r1 = 0.0d0
      r2 = 1.0d0
      print *, kr,kthet,kphi,nin,nout
      print *, r1,r2

      call get_qnodes(kthet,kphi,kr,nin,r1,r2,qin,win)
      call prin2('qin=*',qin,24)
      call prin2('win=*',win,24)
c
c  test win
c
      ra = 0
      do i=1,nin
        ra = ra + win(i)
      enddo
      call prin2('volume of sphere=*',ra,1)
      rex = 4.0d0/3.0d0*pi*(r2**3-r1**3)
      erra = abs(rex-ra)/abs(rex)
      call prin2('error in volume of sphere=*',erra,1)
      print *, ra,rex,ra/rex
      

      r1 = 1
      r2 = 2.0d0
      call get_qnodes(kthet,kphi,kr,nout,r1,r2,qout,wout)
c
c  test win
c
      ra = 0
      do i=1,nin
        ra = ra + wout(i)
      enddo
      call prin2('volume of sphere=*',ra,1)
      erra = abs(4.0d0/3.0d0*pi*(r2**3-r1**3)-ra)/abs(ra)
      call prin2('error in volume of spherical shell=*',erra,1)


      nmax = 50
      allocate(cmevalsinc(3,nmax,nout),cnevalsinc(3,nmax,nout))
      allocate(cmovalsinc(3,nmax,nout),cnovalsinc(3,nmax,nout))
      ifjh = 1

      call prin2('zkout=*',zkout,2)

      call prin2('qout=*',qout,24)
      nout = 1
      qout(1:3,1) = 0
      qout(1,1) = 1
      wout(1) = 1
      call get_cmnvals(nout,qout,nmax,zkout,ifjh,cmevalsinc,cnevalsinc,
     1    cmovalsinc,cnovalsinc)

c
c  test electric field
c
      allocate(cevalsinc(3,nout),chvalsinc(3,nout))
      erra = 0
      ra = 0
      do i=1,nout
        z = qout(3,i)
        etest(1:3) = 0
        etest(1) = exp(ima*zkout*z)
        cevalsinc(1:3,i) = 0
        zfac = ima
        do n=1,nmax
          rfac = (2*n+1.0d0)/(n+0.0d0)/(n+1.0d0)
          cevalsinc(1:3,i) = cevalsinc(1:3,i) +
     1       zfac*rfac*(cmovalsout(1:3,n,i) - 
     1         ima*cnevalsout(1:3,n,i))
          zfac = zfac*ima
        enddo

        if(i.lt.5)
     1   print *, etest(1),cevalsinc(1,i),etest(1)/cevalsinc(1,i)
        if(i.lt.5)
     1   print *, etest(2),cevalsinc(2,i),etest(2)/cevalsinc(2,i)
        
        rmag = abs(etest(1))**2
        ra = ra + rmag*wout(i)

        rmag = abs(etest(1)-cevalsinc(1,i))**2 + abs(cevalsinc(2,i))**2+
     1   abs(cevalsinc(3,i))**2
         
        erra = erra + rmag*wout(i) 
      enddo

      erra = sqrt(erra/ra)
      call prin2('error in incident electric field=*',erra,1)
      
      


      stop
      end
c
c
c
c
c
c
      subroutine get_cmnvals(nq,q,nmax,zk,ifjh,cme,cne,cmo,cno)
      implicit real *8 (a-h,o-z)
      real *8 q(3,nq)
      complex *16 zk,cme(3,nmax,nq),cne(3,nmax,nq),z,ima
      complex *16 cmo(3,nmax,nq),cno(3,nmax,nq)
      real *8, allocatable :: wlege(:)
      real *8, allocatable :: ynm(:,:),ynmd(:,:)
      complex *16, allocatable :: fjs(:),fjder(:),fjdivr(:)
      complex *16 zpsi(3),zynm(3),zphi(3),zn(3)
      complex *16 zr,zt,zp

      data ima/(0.0d0,1.0d0)/


      nlege = nmax + 10
      lw7 = (nlege+1)**2*4
      allocate(wlege(lw7))
      call ylgndrfwini(nlege,wlege,lw7,lused7)


      allocate(fjs(0:nmax+10),fjder(0:nmax+10),fjdivr(0:nmax+10))
      fjs = 0
      fjder = 0
      fjdivr = 0
      rscale = 1.0d0
      ifder = 1
      allocate(ynm(0:nmax,0:nmax),ynmd(0:nmax,0:nmax))

      do i=1,nq
        r = 0
        thet = 0 
        phi = 0
        call cart2polar(q(1,i),r,thet,phi)
        ctheta = cos(thet)
        rx = sin(thet)*cos(phi)
        ry = sin(thet)*sin(phi)
        rz = cos(thet)
        thetx = cos(thet)*cos(phi)
        thety = cos(thet)*sin(phi)
        thetz = -sin(thet)
        phix = -sin(phi)
        phiy = cos(phi)
        phiz = 0
        z = zk*r

        if(ifjh.eq.1) then
          call besseljs3d(nmax+1,z,rscale,fjs,ifder,fjder)
          fjder(0) = fjder(0)*zk
          do n=1,nmax
            fjder(n) = fjder(n)*zk
            fjdivr(n) = fjs(n+1)*rscale + fjs(n-1)/rscale
            fjdivr(n) = fjdivr(n)*zk/(2*n+1.0d0)
          enddo
        else
          call h3dall(nmax+1,z,rscale,fjs,ifder,fjder)
          do n=0,nmax
            fjder(n) = fjder(n)*zk
            fjdivr(n) = fjs(n)/r
          enddo
        endif

        call prin2('fjs=*',fjs,2*nmax)
        call prin2('fjder=*',fjder,24)
        call prin2('fjdivr=*',fjdivr,24)
        call prin2('thet=*',thet,1)
        call prin2('phi=*',phi,1)
        call prin2('ima=*',ima,2)
        ynm = 0
        ynmd = 0
        call ylgndr2sfw(nmax,ctheta,ynm,ynmd,wlege,nlege)
        zn(1:3) = q(1:3,i)/r

        print *, rx,ry,rz
        print *, thetx,thety,thetz
        print *, phix,phiy,phiz

        do n=1,nmax

c
c  compute me
c
          zr = 0
          zt = -sin(phi)*ynm(n,1)*fjs(n)
          zp = cos(phi)*ynmd(n,1)*fjs(n)
          cme(1,n,i) = zr*rx + zt*thetx + zp*phix
          cme(2,n,i) = zr*ry + zt*thety + zp*phiy
          cme(3,n,i) = zr*rz + zt*thetz + zp*phiz


          zr = 0
          zt = cos(phi)*ynm(n,1)*fjs(n)
          zp = sin(phi)*ynmd(n,1)*fjs(n)
          cmo(1,n,i) = zr*rx + zt*thetx + zp*phix
          cmo(2,n,i) = zr*ry + zt*thety + zp*phiy
          cmo(3,n,i) = zr*rz + zt*thetz + zp*phiz


          zr = fjdivr(n)*cos(phi)*(n+0.0d0)*(n+1.0d0)*ynm(n,1)*sin(thet)
          zt = -cos(phi)*ynmd(n,1)*(fjdivr(n)/zk + fjder(n)/zk)
          zp = -sin(phi)*ynm(n,1)*(fjdivr(n)/zk + fjder(n)/zk)
          cne(1,n,i) = zr*rx + zt*thetx + zp*phix
          cne(2,n,i) = zr*ry + zt*thety + zp*phiy
          cne(3,n,i) = zr*rz + zt*thetz + zp*phiz

          zr = fjdivr(n)*sin(phi)*(n+0.0d0)*(n+1.0d0)*ynm(n,1)*sin(thet)
          zt = -sin(phi)*ynmd(n,1)*(fjdivr(n)/zk + fjder(n)/zk)
          zp = cos(phi)*ynm(n,1)*(fjdivr(n)/zk + fjder(n)/zk)
          cno(1,n,i) = zr*rx + zt*thetx + zp*phix
          cno(2,n,i) = zr*ry + zt*thety + zp*phiy
          cno(3,n,i) = zr*rz + zt*thetz + zp*phiz
        enddo
      enddo



      return
      end
c
c
c
c
c
c
      subroutine get_qnodes(kt,kp,kr,n,r1,r2,q,w)
      implicit real *8 (a-h,o-z)
      real *8 q(3,n),w(n)
      real *8, allocatable :: xt(:),wt(:),ut(:,:),vt(:,:)
      real *8, allocatable :: xr(:),wr(:),ur(:,:),vr(:,:)

      done = 1.0d0
      pi = atan(done)*4

      allocate(xt(kt),wt(kt),ut(kt,kt),vt(kt,kt))
      allocate(xr(kr),wr(kr),ur(kr,kr),vr(kr,kr))

      itype = 1
      call legeexps(itype,kt,xt,ut,vt,wt)
      call legeexps(itype,kr,xr,ut,vt,wr)


      i = 0
      do i1=1,kr
        r = r1 + (r2-r1)*(xr(i1)+1)/2
        do i2 = 1,kt
          thet = pi*(xt(i2)+1)/2
          do i3 = 1,kp
            i = i+1
            phi = 2*pi*(i3-1.0d0)/(kp+0.0d0)
            q(1,i) = r*sin(thet)*cos(phi)
            q(2,i) = r*sin(thet)*sin(phi)
            q(3,i) = r*cos(thet)
            w(i) = wt(i2)*pi/2*wr(i1)*(r2-r1)/2*2*pi/(kp+0.0d0)*
     1         r**2*sin(thet)

          enddo
        enddo
      enddo



      

      return
      end
