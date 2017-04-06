c gauss_subs.f Subroutine for the gauss_create package.
c April 17, 2002 Version.
c
c invert.f
c concat.f
!=========================================================================
      SUBROUTINE INVERT(A,NN,N,M,C)
!
! INPUT A(N,N) THIS IS DESTROYED AND AINV GOES INTO A.
! M(N,N), C(N,N) ARE WORK SPACE MATRICES.
! N=DIM OF A,M,C IN CALLING ROUTINE.
! NN=DIM OF MATRIX TO INVERT.
!          GENERAL MATRIX INVERSION SUBROUTINE
!
 
        implicit none
 
        real*8 a,c,d,temp
 
        integer nn,n,m,itemp,i,j,k,kd,l,ld,nh,nr
 
        DIMENSION A(*),M(*),C(*)
!
      if(nn .le. 0)then
        write(*,*)'nn= ',nn,' in invert!'
        stop
      else if(nn .eq.1)then
        A(1)=1.0d0/A(1)
        return
      end if
      do I=1,NN
        M(I)=-I
      end do
!
      DO 140 I = 1,NN
!
!     LOCATE LARGEST ELEMENT
!
        D=0.0d0
        DO 112 L=1,NN
          IF (M(L)) 100,100,112
  100     J=L
          DO 110 K=1,NN
            IF (M(K)) 103,103,108
  103       IF (ABS(D)-ABS(A(J))) 105,105,108
  105       LD=L
            KD=K
            D=A(J)
  108       J=J+N
  110     CONTINUE
  112   CONTINUE
!
!       INTERCHANGE ROWS
!
        ITEMP=-M(LD)
        M(LD)=M(KD)
        M(KD)=ITEMP
        L=LD
        K=KD
        DO 114 J=1,NN
          C(J)=A(L)
          A(L)=A(K)
          A(K)=C(J)
          L=L+N
          K=K+N
  114   continue
!
!       DIVIDE COLUNM BY LARGEST ELEMENT
!
        NR=(KD-1)*N+1
        NH=NR+N-1
        DO 115 K=NR,NH
          A(K)=A(K)/D
  115   continue
!
!       REDUCE REMAINING ROWS AND COLUMNS
!
        L=1
        DO 135 J=1,NN
          IF (J-KD) 130,125,130
  125     L=L+N
          GO TO 135
  130     DO 134 K=NR,NH
            A(L)=A(L)-C(J)*A(K)
            L=L+1
  134     continue
  135   CONTINUE
!
!       REDUCE ROW
!
        C(KD)=-1.0D+00
        J=KD
        DO 140 K=1,NN
          A(J)=-C(K)/D
          J=J+N
  140 continue
!
!     INTERCHANGE COLUMNS
!
      DO 200 I=1,NN
        L=0
  150   L=L+1
        IF(M(L)-I) 150,160,150
  160   K=(L-1)*N+1
        J=(I-1)*N+1
        M(L)=M(I)
        M(I)=I
        DO 200 L=1,NN
          TEMP=A(K)
          A(K)=A(J)
          A(J)=TEMP
          J=J+1
          K=K+1
  200 continue
!
      RETURN
      END
c ===================================================
c
      subroutine concat(a,b,c)
c concatenate words a and b to form word c,
c but first strips off trailing blanks.
c
      implicit none
c
c...passed and global variables
      character*(*) a,b,c
c
c...local variables
      integer lna,lnb
c
c . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
c
      lna = index(a," ")-1
      if (lna.eq.-1) lna = len(a)
      lnb = index(b," ")-1
      if (lnb.eq.-1) lnb = len(b)
c
      if (lna.eq.0) then
        c = b
      elseif (lnb.eq.0) then
        c = a
      else
        c = a(1:lna) // b(1:lnb)
      endif
      return
        end
c ======================================================
