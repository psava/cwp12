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
!
! Esteban's comment:
! Nocedal's license is very flexible, it only ask to reference
! his '82 paper about the method.
!
program sdrive 
  use rsf

  implicit none 
  integer      :: NDIM,MSAVE,NWORK,cit 
  double precision, allocatable,dimension (:) :: X,G,DIAG,W
  double precision :: F,EPS,XTOL,GTOL,T1,T2,STPMIN,STPMAX
  integer      :: IPRINT(2),IFLAG,ICALL,N,M,MP,LP,J,ITER
  logical      :: DIAGCO
  real ::  Fo, grad,xo,yo,sx,sy
  EXTERNAL LB2
  COMMON /LB3/MP,LP,GTOL,STPMIN,STPMAX


  NDIM=2 ;MSAVE=5;NWORK=NDIM*(2*MSAVE +1)+2*MSAVE
  
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
  M=5
  IPRINT(1)= -1
  IPRINT(2)= -1
!
!     We do not wish to provide the diagonal matrices Hk0, and 
!     therefore set DIAGCO to FALSE.
!
  DIAGCO= .FALSE.
  EPS= 1.0D-7
  XTOL= 1.0D-16
  ICALL=0
  IFLAG=0

  xo=100.0; yo=100.0;
  sx=180; sy=180;

  X(1)=-4000.0; X(2)=+1293100.0 ;
 

  
  cit=1; ITER=0
  WRITE(0, "('iter: ',I3,' x= ',F15.7,' y=',F15.7)")ITER,X(1),X(2)

  do while (ICALL.lt.2000)
    F= 0.D0
    F= Fo(real(X(1)),real(X(2)))
    call gradient(G,real(X(1)),real(X(2)))

    cit=ITER
    CALL LBFGS(NDIM,M,X,F,G,DIAGCO,DIAG,IPRINT,EPS,XTOL,W,IFLAG,ITER)
    ICALL=ICALL + 1

    if(IFLAG.le.0 ) exit
    if(cit.ne.ITER) WRITE(0, "('iter: ',I3,' x= ',F15.7,' y=',F15.7,I3)")ITER,X(1),X(2),ICALL
    if(cit.ne.ITER) WRITE(0,*) W
  end do

end program sdrive 

function Fo(x,y)
    real x,y;
    Fo = x*x +y*y 
    return 
end function 

subroutine gradient(grad,x,y)
    real x,xo,sx,y,yo,sy;
    double precision::  grad(2) 
    real Fo;

    grad(1)= 2*x
    grad(2)= 2*y
end subroutine 
