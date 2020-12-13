      implicit real *8 (a-h,o-z)
      complex *16, allocatable :: cmvalsin(:,:,:),cnvalsin(:,:,:)
      complex *16, allocatable :: cmvalsout(:,:,:),cnvalsout(:,:,:)
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
      allocate(cmvalsout(3,nmax,nout),cnvalsout(3,nmax,nout))
      ifjh = 1

      call prin2('zkout=*',zkout,2)

      call prin2('qout=*',qout,24)
      nout = 5
      qout(1:3,1) = 0
      qout(1,1) = 1
      wout(1) = 1
      call get_cmnvals(nout,qout,nmax,zkout,ifjh,cmvalsout,cnvalsout)

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
     1       zfac*rfac*(imag(cmvalsout(1:3,n,i)) - 
     1         ima*real(cnvalsout(1:3,n,i)))
          print *, n, zfac,rfac
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
      subroutine get_cmnvals(nq,q,nmax,zk,ifjh,cm,cn)
      implicit real *8 (a-h,o-z)
      real *8 q(3,nq)
      complex *16 zk,cm(3,nmax,nq),cn(3,nmax,nq),z,ima
      real *8, allocatable :: wlege(:)
      real *8, allocatable :: ynm(:,:),ynmd(:,:)
      complex *16, allocatable :: fjs(:),fjder(:),fjdivr(:)
      complex *16 rpsi(3),ry(3),rphi(3),rn(3)
      complex *16 rpsithet,rpsiphi

      data ima/(0.0d0,1.0d0)/


      nlege = nmax + 10
      lw7 = (nlege+1)**2*4
      allocate(wlege(lw7))
      call ylgndrfwini(nlege,wlege,lw7,lused7)


      allocate(fjs(0:nmax+10),fjder(0:nmax+10),fjdivr(0:nmax+10))
      rscale = 1.0d0
      ifder = 1
      allocate(ynm(0:nmax,0:nmax),ynmd(0:nmax,0:nmax))

      do i=1,nq
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

        call prin2('fjs=*',fjs,24)
        call prin2('fjder=*',fjder,24)
        call prin2('fjdivr=*',fjdivr,24)
        call prin2('thet=*',thet,1)
        call prin2('phi=*',phi,1)

        call ylgndr2sfw(nmax,ctheta,ynm,ynmd,wlege,nlege)
        rn(1:3) = q(1:3,i)/r
        print *, thetx,thety,thetz
        print *, phix,phiy,phiz

        do n=1,nmax
          rpsithet = ynmd(n,1)*exp(ima*phi)
          rpsiphi = ima*exp(ima*phi)*ynm(n,1)
          ry(1:3) = rn(1:3)*ynm(n,1)*sin(thet)*exp(ima*phi)
          rpsi(1) = rpsithet*thetx + rpsiphi*phix
          rpsi(2) = rpsithet*thety + rpsiphi*phiy
          rpsi(3) = rpsithet*thetz + rpsiphi*phiz
          call zcross_prod3d(rn,rpsi,rphi)

          cm(1:3,n,i) = -rphi(1:3)*fjs(n)*rscale**n
          cn(1:3,n,i) = ry(1:3)*n*(n+1.0d0)*fjdivr(n) + 
     1        rpsi(1:3)*(fjder(n) + fjdivr(n)) 
          cn(1:3,n,i) = cn(1:3,n,i)*rscale**n
        enddo
        call prin2('cn=*',cn,24)
        call prin2('cm=*',cm,24)
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
