! copyright info:
!
!                             @Copyright 1999
!                          Fireball2000 Committee
!
! ASU - Otto F. Sankey
!       Kevin Schmidt
!       Jian Jun Dong
!       John Tomfohr
!       Gary B. Adams
!
! Motorola - Alex A. Demkov
!
! University of Regensburg - Juergen Fritsch
!
! University of Utah - James P. Lewis
!                      Kurt R. Glaesemann
!
! Universidad de Madrid - Jose Ortega
!
!                  made possible by support from Motorola
!                         fireball@fermi.la.asu.edu

! RESTRICTED RIGHTS LEGEND
!
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in
! subdivsion { (b) (3) (ii) } of the Rights in Technical Data and
! Computer Software clause at 52.227-7013.

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
