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
      complex *16 etest(3),zfac,htest(3)
      real *8 v1(3),v2(3),v3(3)
      real *8 rnx,rdx,rix
      data ima/(0.0d0,1.0d0)/


      call prini(6,13)
      
      done = 1.0d0
      pi = atan(done)*4.0d0
      rlam = 2.0d0
      zk = rlam*2*pi/2



      reta = 1.33
      reta = 1.1d0
      reta = 2
      zkout = zk
      zkin = zkout*reta

      kthet = 60
      kphi = 60
      kr = 30
      nin = kthet*kphi*kr
      allocate(qin(3,nin),win(nin))

      nout = kthet*kphi*kr
      allocate(qout(3,nout),wout(nout))

      r1 = 0.0d0
      r2 = 1.0d0

      call get_qnodes(kthet,kphi,kr,nin,r1,r2,qin,win)
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


      nmax = 100
      allocate(cmevalsinc(3,nmax,nout),cnevalsinc(3,nmax,nout))
      allocate(cmovalsinc(3,nmax,nout),cnovalsinc(3,nmax,nout))
      ifjh = 1

      call get_cmnvals(nout,qout,nmax,zkout,ifjh,cmevalsinc,cnevalsinc,
     1    cmovalsinc,cnovalsinc)

c
c  test electric field
c
      allocate(cevalsinc(3,nout),chvalsinc(3,nout))
      erra = 0
      ra = 0
      errah = 0
      do i=1,nout
        z = qout(3,i)
        etest(1:3) = 0
        htest(1:3) = 0
        etest(1) = exp(ima*zkout*z)
        htest(2) = exp(ima*zkout*z)
        cevalsinc(1:3,i) = 0
        chvalsinc(1:3,i) = 0
        zfac = ima
        do n=1,nmax
          rfac = (2.0d0*n+1.0d0)/(n+0.0d0)/(n+1.0d0)
          cevalsinc(1:3,i) = cevalsinc(1:3,i) +
     1       zfac*rfac*(cmovalsinc(1:3,n,i) - 
     1         ima*cnevalsinc(1:3,n,i))
          chvalsinc(1:3,i) = chvalsinc(1:3,i) + 
     1      zfac*rfac*(cmevalsinc(1:3,n,i) + ima*
     2        cnovalsinc(1:3,n,i))
          zfac = zfac*ima
        enddo

        rmag = abs(etest(1))**2
        ra = ra + rmag*wout(i)

        rmag = abs(etest(1)-cevalsinc(1,i))**2 + abs(cevalsinc(2,i))**2+
     1   abs(cevalsinc(3,i))**2
         
        erra = erra + rmag*wout(i)
        rmag = abs(chvalsinc(1,i))**2 + abs(htest(2)-chvalsinc(2,i))**2+
     1   abs(chvalsinc(3,i))**2
        errah = errah + rmag*wout(i)
         
      enddo

      erra = sqrt(erra/ra)
      call prin2('error in incident electric field=*',erra,1)
      
      errah = sqrt(errah/ra)
      call prin2('error in incident magnetic field=*',erra,1)


      allocate(ccoefs(nmax),dcoefs(nmax))

      call get_trans_coeff(nmax,zkout,zkin,reta,ccoefs,dcoefs)
      call prin2('ccoefs=*',ccoefs,2*nmax)
      call prin2('dcoefs=*',dcoefs,2*nmax)
c
c
c
c
      allocate(cmevalsin(3,nmax,nin),cnevalsin(3,nmax,nin))
      allocate(cmovalsin(3,nmax,nin),cnovalsin(3,nmax,nin))
      allocate(cevalsin(3,nin),chvalsin(3,nin))
      ifjh = 1

      nn = 100
      do ii=1,nn
        reta = 1.0d0 + (ii-1)/(nn-1.0d0)
        zkout = zk
        zkin = zkout*reta
        call prin2('zkin=*',zkin,2)
        call prin2('zkout=*',zkout,2)
        ccoefs = 0
        dcoefs = 0
        call get_trans_coeff(nmax,zkout,zkin,reta,ccoefs,dcoefs)
        call prin2('ccoefs=*',ccoefs,2*nmax)
        call prin2('dcoefs=*',dcoefs,2*nmax)



        cmevalsin = 0
        cnevalsin = 0
        cmovalsin = 0
        cnovalsin = 0

        call prinf('ifjh=*',ifjh,1)

        call get_cmnvals(nin,qin,nmax,zkin,ifjh,cmevalsin,cnevalsin,
     1      cmovalsin,cnovalsin)

        

        rint = 0
        rrint = 0
        riint = 0
        reint = 0
        rsint = 0
        cevalsin = 0
        chvalsin = 0
        rvol = 0
        rmax = 0
        imax = 1
        do i=1,nin
          zfac = ima
          do n=1,nmax
            rfac = (2.0d0*n+1.0d0)/(n+0.0d0)/(n+1.0d0)
            cevalsin(1:3,i) = cevalsin(1:3,i) +
     1         zfac*rfac*(ccoefs(n)*cmovalsin(1:3,n,i) - 
     1         ima*dcoefs(n)*cnevalsin(1:3,n,i))
            chvalsin(1:3,i) = chvalsin(1:3,i) - 
     1        zfac*rfac*reta*(dcoefs(n)*cmevalsin(1:3,n,i) + ima*
     2        ccoefs(n)*cnovalsin(1:3,n,i))
            zfac = zfac*ima
          enddo
          v1(1:3) = real(cevalsin(1:3,i))
          v2(1:3) = imag(cevalsin(1:3,i))
          call cross_prod3d(v1,v2,v3)
          rmag3 = v3(1)**2 + v3(2)**2 + v3(3)**2 
          rint = rint + rmag3*win(i)
          rmag1 = v1(1)**2 + v1(2)**2 + v1(3)**2
          rrint = rrint + rmag1*win(i)
          rmag2 = v2(1)**2 + v2(2)**2 + v2(3)**2
          riint = riint + rmag2*win(i)


          r1 = sqrt(rmag3)/sqrt(rmag1*rmag2)
          if(rmag1*rmag2.le.1.0d-32) r1 = 0

          rsint = rsint + r1*win(i)


          rnx = rmag1 + rmag2
          rdx = rmag1 - rmag2

          rix = v1(1)*v2(1) + v1(2)*v2(2) + v1(3)*v2(3)

          rdis = sqrt(rdx**2 + 4*rix**2)


          rex = (rnx - rdis)/(rnx + rdis)

          reint = reint + rex*win(i)

          rvol = rvol + win(i)
          

        enddo

        print *, imax,rmax
        rint = sqrt(rint)
        rrint = sqrt(rrint)
        riint = sqrt(riint)
        call prin2('rvol=*',rvol,1)

        rint = rint/rrint/riint
        rsint = sqrt(rsint)/rvol
        reint = sqrt(reint)/rvol
        write(40,*) reta,rint,rsint,reint,rrint,riint
        call prin2('l2 norm of integral of cross product=*',rint,1)
        call prin2('rel l2 norm of sin=*',rsint,1)
        call prin2('rel l2 norm of ecc=*',reint,1)
      enddo
      


      stop
      end
c
c
c
c
c
c
      subroutine get_trans_coeff(nmax,zk0,zk1,reta,ccoefs,dcoefs)
      implicit real *8 (a-h,o-z)
      complex *16 zk0,zk1,ccoefs(nmax),dcoefs(nmax)
      complex *16, allocatable :: fjs0(:),fhs0(:),fjder0(:),fhder0(:)
      complex *16, allocatable :: fjs1(:),fhs1(:),fjder1(:),fhder1(:)
      complex *16 z1,z2,z3,z4

      allocate(fjs0(0:nmax+10),fhs0(0:nmax+10))
      allocate(fjder0(0:nmax+10),fhder0(0:nmax+10))


      allocate(fjs1(0:nmax+10),fhs1(0:nmax+10))
      allocate(fjder1(0:nmax+10),fhder1(0:nmax+10))

      ifder = 1
      rscale = 1.0d0
      call besseljs3d(nmax+1,zk0,rscale,fjs0,ifder,fjder0)
      call h3dall(nmax+1,zk0,rscale,fhs0,ifder,fhder0)

      call prin2('fjs0=*',fjs0,2*nmax)
      call prin2('fhs0=*',fhs0,2*nmax)


      call besseljs3d(nmax+1,zk1,rscale,fjs1,ifder,fjder1)
      call h3dall(nmax+1,zk1,rscale,fhs1,ifder,fhder1)
      call prin2('fjs1=*',fjs1,2*nmax)
      call prin2('fhs1=*',fhs1,2*nmax)
      print *, zk0
      print *, zk1

      do n=1,nmax
        z1 = (zk0*fhder0(n) + fhs0(n))*fjs0(n)
        z2 = (zk0*fjder0(n) + fjs0(n))*fhs0(n)
        z3 = (zk0*fhder0(n) + fhs0(n))*fjs1(n)
        z4 = (zk1*fjder1(n) + fjs1(n))*fhs0(n)

        print *, n, abs(z3-z4),abs(reta**2*z3-z4)


        ccoefs(n) = (z1-z2)/(z3-z4)
        dcoefs(n) = (reta*z1 - reta*z2)/(reta**2*z3 - z4)
      enddo


      return
      end



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

        ynm = 0
        ynmd = 0
        call ylgndr2sfw(nmax,ctheta,ynm,ynmd,wlege,nlege)

        do n=1,nmax
          ynm(n,1) = -ynm(n,1)*sqrt((n+0.0d0)*(n+1.0d0)/(2*n+1.0d0))
          ynmd(n,1) = -ynmd(n,1)*sqrt((n+0.0d0)*(n+1.0d0)/(2*n+1.0d0))

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
          zr = zr/zk
          zt = -cos(phi)*ynmd(n,1)*(fjdivr(n)/zk + fjder(n)/zk)
          zp = -sin(phi)*ynm(n,1)*(fjdivr(n)/zk + fjder(n)/zk)
          cne(1,n,i) = zr*rx + zt*thetx + zp*phix
          cne(2,n,i) = zr*ry + zt*thety + zp*phiy
          cne(3,n,i) = zr*rz + zt*thetz + zp*phiz

          zr = fjdivr(n)*sin(phi)*(n+0.0d0)*(n+1.0d0)*ynm(n,1)*sin(thet)
          zr = zr/zk
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
