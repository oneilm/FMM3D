      subroutine hndiv(dbsize,eps,ns,nt,ifcharge,ifdipole,ifpgh,
     1   ifpghtarg,ndiv,idivflag)
c
c
c       this subroutine estimates ndiv and idivflag 
c       based on geometry parameters
c     
c       ndiv depends only on the precision. It is not increased
c       with the size of the domain in wavelengths (dbsize, which
c       is unused): with the diagonal form of the multipole to
c       local translation at high frequency, a deeper tree is
c       cheaper than the extra direct interactions of large leaf
c       boxes.
c
c
      implicit none
      real *8 eps,dbsize
      integer *8 ns,nt,ifcharge,ifdipole,ifpgh,ifpghtarg,ndiv
      integer *8 idivflag

      idivflag = 0

       if(eps.ge.0.5d-0) then
         ndiv = 40
       else if(eps.ge.0.5d-1) then
         ndiv = 40
       else if(eps.ge.0.5d-2) then
         ndiv = 40
       else if(eps.ge.0.5d-3) then
         ndiv = 120
       else if(eps.ge.0.5d-6) then
         ndiv = 300
       else if(eps.ge.0.5d-9) then
         ndiv = 600
       else if(eps.ge.0.5d-12) then
         ndiv = 900
       else if(eps.ge.0.5d-15) then
         ndiv = 1050
       else
         ndiv = ns+nt
       endif

       ndiv = min(ndiv,10000)

      return
      end
