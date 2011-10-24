!  FORTRAN90  driver adapted from:
!
!     ***********************
!     SIMPLE DRIVER FOR LBFGS
!     ***********************
!
!     Example of driver for LBFGS routine, using a
!     simple test problem. The solution point is at 
!     X=(1,...,1) and the optimal function value of 0.
!
!                          JORGE NOCEDAL
!                        *** July 1990 ***
!

program sdrive 
  use rsf

  implicit none 
  integer      :: NDIM,MSAVE,NWORK 
  double precision, allocatable,dimension (:) :: X,G,DIAG,W
  double precision :: F,EPS,XTOL,GTOL,T1,T2,STPMIN,STPMAX
  integer      :: IPRINT(2),IFLAG,ICALL,N,M,MP,LP,J,ITER
  logical      :: DIAGCO
  EXTERNAL LB2
  COMMON /LB3/MP,LP,GTOL,STPMIN,STPMAX



  NDIM=2000 ;MSAVE=7;NWORK=NDIM*(2*MSAVE +1)+2*MSAVE
  
!Dictionary of variables:
! NDIM    :   number of variables.
! M       :   number of corrections for BFGS update Nocedal recomends 3<= M <=7
! X       :   state variable (passed as input to LBFGS 
!                             and its returned updated from it)
! F       :   value of the objective function
! G       :   gradient passed to the LBFGS ROUTINE
! DIAGCO  :   logic variable (.true. means user is going to pass the diagonal of
!             the HESSIAN, .false. means it won't pass it --the normal thing to do)
! DIAG    :   DIAGONAL of the HESSIAN 0 if 





  allocate(X(NDIM), G(NDIM), DIAG(NDIM), W(NWORK))
  

!      The driver for LBFGS must always declare LB2 as EXTERNAL
!
!
  N=100
  M=5
  IPRINT(1)= 1
  IPRINT(2)= 0
!
!     We do not wish to provide the diagonal matrices Hk0, and 
!     therefore set DIAGCO to FALSE.
!
  DIAGCO= .FALSE.
  EPS= 1.0D-5
  XTOL= 1.0D-16
  ICALL=0
  IFLAG=0

  do j=1,N,2
     X(J)=-1.2D0
     X(J+1)=1.D0
  end do
  
  do while (ICALL.lt.2000 ) 	
  	 		 
    F= 0.D0
    do j=1,N,2
      T1= 1.D0-X(J)
      T2= 1.D1*(X(J+1)-X(J)**2)
      G(J+1)= 2.D1*T2
      G(J)= -2.D0*(X(J)*G(J+1)+T1)
      F= F+T1**2+T2**2
	end do
    CALL LBFGS(N,M,X,F,G,DIAGCO,DIAG,IPRINT,EPS,XTOL,W,IFLAG,ITER)
    ICALL=ICALL + 1
    !WRITE(0,*)'ICALL:',ICALL,'ITER:',ITER


	if(IFLAG.eq.0 ) exit
  end do
  

end program sdrive 	 