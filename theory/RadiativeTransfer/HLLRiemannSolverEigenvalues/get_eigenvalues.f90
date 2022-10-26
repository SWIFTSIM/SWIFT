!-------------------------------------------------------
! Compute the Eigenvalues for the M1 radiative
! transfer HLL Riemann solver.
!
! Original file sent by  Matthias González, was
! used in HERACLES code 
! (ui.adsabs.harvard.edu/abs/2007A%26A...464..429G)
! under the CeCILL licence, which is compatible with
! the GNU GPL licence.
!-------------------------------------------------------
program compute_RT_riemann_solver_eigenvalues
  implicit none

  integer, parameter  :: dp = SELECTED_REAL_KIND(15)
  integer, parameter  :: ndim = 3
  real(dp), parameter :: pi = 4.D0*DATAN(1.d0)
  real(dp), parameter :: c = 1.0

  ! How many points to go through?
  ! NOTE: we always also include zero, so the resulting
  ! arrays will have size n_points + 1
  integer, parameter :: n_points = 99

  ! Make this array allocatable because it might be too
  ! large for the stack.
  real(dp), dimension(:,:,:), allocatable :: eigenvals ! eigenvalues.
  allocate(eigenvals(0:n_points, 0:n_points, 1:4))

  call compute_eigenvals()
  call write_eigenvals()


contains

!  Subroutine COMPUTE_eigenvals
!
!> Computes the tabulated eigenvalues for the M1 radiative
!! transfer solver.
!<
subroutine compute_eigenvals
  
  implicit none
  integer  :: i,j
  real(dp) :: normef ! |F|
  real(dp) :: theta

  ! Arguments needed for external eigenvalue calculation routine
  real(dp), dimension(4)      :: WR     
  real(dp), dimension(4)      :: WI
  integer,parameter           :: LDVL=1
  integer,parameter           :: LDVR=1
  real(dp), dimension(LDVL,4) :: VL
  real(dp), dimension(LDVR,4) :: VR
  integer                     :: info
  integer                     :: lda=4
  integer,parameter           :: lwork=20
  real(dp),dimension(lwork)   :: work

  real(dp)                    :: E        ! Radiation energy
  real(dp),dimension(1:3    ) :: F        ! Radiation flux
  real(dp),dimension(3,3    ) :: Dedd     ! Eddington tensor
  real(dp),dimension(3,3    ) :: Dedd_dE  ! dDedd/dE * E
  real(dp),dimension(3,3,1:3) :: Dedd_dF  ! dDedd/dF * E
  real(dp), dimension(4,4)    :: mat      ! Jacobi matrix

  E=1. ! Photon Energy
  F=0. ! Photon Flux

  do i=0,n_points
     ! loop over ||f||
     normef=float(i)/n_points

     do j=0,n_points
        ! loop over 0 < theta < pi: angle w.r.t. interface
        theta=j*pi/n_points
        
        F(1)=c*E*cos(theta)*normef
        F(2)=c*E*sin(theta)*normef
        
        ! Compute Eddington tensor and its derivatives
        ! (derivatives are multiplied by Ei)
        call cal_Dedd(E,F,Dedd,Dedd_dE,Dedd_dF)

        ! Get Jacobi Matrix (for (E, Fx) only, but using
        ! gradients in up to all 3 dimensions)
        mat=0.

        mat(1,1) = 0.
        mat(1,2) = 1.
        mat(2,1) = c**2*(Dedd(1,1)+Dedd_dE(1,1  ))
        mat(2,2) = c**2*           Dedd_dF(1,1,1)

        if(ndim.gt.1) then
           mat(1,3) = 0.
           mat(2,3) = c**2*           Dedd_dF(1,1,2)

           mat(3,1) = c**2*(Dedd(1,2)+Dedd_dE(1,2  ))
           mat(3,2) = c**2*           Dedd_dF(1,2,1)
           mat(3,3) = c**2*           Dedd_dF(1,2,2)
        endif
        if(ndim.gt.2) then
           mat(1,4) = 0.
           mat(2,4) = c**2*           Dedd_dF(1,1,3)
           mat(3,4) = c**2*           Dedd_dF(1,2,3)
         
           mat(4,1) = c**2*(Dedd(1,3)+Dedd_dE(1,3  ))
           mat(4,2) = c**2*           Dedd_dF(1,3,1)
           mat(4,3) = c**2*           Dedd_dF(1,3,2)
           mat(4,4) = c**2*           Dedd_dF(1,3,3)
        endif

        ! normalized the matrix to get eigenvalues in units of c
        mat(2:4,1)=mat(2:4,1)/c**2
        mat(2:4,2:4)=mat(2:4,2:4)/c

        ! DGEEV computes the eigenvalues and, optionally, the left 
        ! and/or right eigenvectors for GE matrices. It's a function
        ! from LAPACK
        ! Warning : dsyev works only for symmetric matrix -> use dgeev instead
        call dgeev('N',  & ! Don't compute left eigenvectors
                   'N',  & ! Don't compute right eigenvectors
                   4,    & ! Order of matrix
                   mat,  & ! Matrix for which to calculate the Eigenvalues
                   lda,  & ! leading dimension of array mat. (which index to start with)
                   WR,   & ! Real parts of eigenvalues
                   WI,   & ! Imaginary parts of eigenvalues
                   VL,   & ! array for left eigenvectors (unused)
                   LDVL, & ! leading dimension for VL array
                   VR,   & ! right eigenvectors (unused)
                   LDVR, & ! leading dimension for VR array
                   WORK, & ! workspace for this routine to work with. Should be empty array
                   LWORK,& ! Dimension of array WORK. LWORK >= max(1, 3*N)
                   info  & ! error code. Success if info = 0
                   )

         ! sort the real part of the eigenvalues in ascending order
         !
         ! warning : the eigenvalues may be not real but complex...
         !           in particular it happens when eps->0 and f->1
         ! note for Mladen: I took out the eps part
         call bubble_sort(4,WR)
         eigenvals(i,j,1:4)=WR(1:4)
     enddo
  enddo

  return

end subroutine compute_eigenvals

!###########################################################
!###########################################################
!###########################################################
!###########################################################

subroutine bubble_sort(n,arr)

  implicit none

  integer                :: n,i,j
  real(dp)               :: a
  real(dp), dimension(n) :: arr

  do j = 1,n
     do i = 2,n
        if(arr(i) < arr(i-1))then
           a        = arr(i  )
           arr(i  ) = arr(i-1)
           arr(i-1) = a
        endif
     enddo
  enddo

  return

end subroutine bubble_sort

!###########################################################
!###########################################################
!###########################################################
!###########################################################

!  Subroutine CAL_DEDD
!
!> Computes the Eddington tensor Dedd and its derivatives
!! with respect to Er and Fr, Dedd_dE and Dedd_dF,
!! respectively.
!! We'll need them to compute the Jacobi-Matrix dF/dU,
!! where 
!!    U = (E_photon, F_photon); 
!!    F = (F_photon, c^2 P)
!!    P = D E_photon.
!! The derivatives of P will be done by part, e.g.
!!    dP/dE_photon = D + dD/dE * E
!! which is why we directly multiply the derivatives by E.
!<
subroutine cal_Dedd(E,F,Dedd,Dedd_dE,Dedd_dF)

  implicit none

  real(dp)                      , intent(in ) :: E  ! Radiation energy
  real(dp),dimension(1:ndim    ), intent(in ) :: F  ! Radiation flux
  real(dp),dimension(3,3       ), intent(out) :: Dedd ! Eddington tensor
  real(dp),dimension(3,3       ), intent(out) :: Dedd_dE ! dD / dE * Ei
  real(dp),dimension(3,3,1:ndim), intent(out) :: Dedd_dF ! dD / dF * Ei

  integer :: i
  real(dp):: normeF   ! || F ||
  real(dp):: chi      ! chi parameter for Eddington tensor
  real(dp):: chiprime ! dchi/df
  real(dp):: g        ! diagonal part of Edd. tensor
  real(dp):: h        ! anisotropic part of Edd. tensor
  real(dp):: gprime   ! dg / df
  real(dp):: hprime   ! dh / df
  real(dp):: nx       ! x component of unit vector
  real(dp):: ny       ! y component of unit vector
  real(dp):: nz       ! z component of unit vector
  real(dp):: fx       ! reduced flux x component
  real(dp):: fy       ! reduced flux y component
  real(dp):: fz       ! reduced flux z component
  real(dp):: cE       ! c * E

  Dedd = 0.d0 ; Dedd_dE = 0.d0 ; Dedd_dF = 0.d0
  fx = 0.d0 ; fy = 0.d0 ; fz = 0.d0
  cE = c * E

                fx = F(1)/cE
  if(ndim.gt.1) fy = F(2)/cE
  if(ndim.gt.2) fz = F(3)/cE

  normef = 0.
  do i=1,ndim
     normef = normef + (F(i)/cE)**2
  enddo
  normef = sqrt(normef)

  if(normef.gt.1) normef = 1.d0 ! For chi, f <= 1 is required

  chi      = (3.d0+4.d0*normef**2) / (5.d0+2.d0*sqrt(4.d0-3.d0*normef**2))
  chiprime = 2.d0*normef / sqrt(4.d0 - 3.d0 * normef**2) ! dchi/df

  g      = 0.5d0 * (1.d0 - chi)      ! diagonal part of tensor
  h      = 0.5d0 * (3.d0 * chi-1.d0) ! anisotropic part of tensor
  gprime = -chiprime * 0.5d0         ! dg/df
  hprime = 3.d0 * chiprime * 0.5d0   ! dh/df

  nx = 0.d0 ; ny = 0.d0 ; nz = 0.d0  ! unit vector elements

  ! n is a unit vector
  normef = 0.d0
  do i=1,ndim
    normef = normef + (F(i)/cE)**2
  enddo
  normef = sqrt(normef)

  if(normef.gt.1.e-8) then
    nx = fx/normef
    ny = fy/normef
    nz = fz/normef
  endif

  ! Eddington tensor
  Dedd(1,1) = g + h*nx*nx ; Dedd(1,2) =     h*nx*ny ; Dedd(1,3) =     h*nx*nz
  Dedd(2,1) =     h*ny*nx ; Dedd(2,2) = g + h*ny*ny ; Dedd(2,3) =     h*ny*nz
  Dedd(3,1) =     h*nz*nx ; Dedd(3,2) =     h*nz*ny ; Dedd(3,3) = g + h*nz*nz

  ! (Eddington tensor derivative with respect to Ei) * Ei
  ! remember dchi/dE = dchi/df df/dE = - dchi/df |F|/(cE^2) = -dchi/df * f / E
  Dedd_dE(1,1) = - normef * (gprime + hprime*nx*nx)
  Dedd_dE(1,2) = - normef * (         hprime*nx*ny)
  Dedd_dE(1,3) = - normef * (         hprime*nx*nz)
  Dedd_dE(2,1) = - normef * (         hprime*ny*nx)
  Dedd_dE(2,2) = - normef * (gprime + hprime*ny*ny)
  Dedd_dE(2,3) = - normef * (         hprime*ny*nz)
  Dedd_dE(3,1) = - normef * (         hprime*nz*nx)
  Dedd_dE(3,2) = - normef * (         hprime*nz*ny)
  Dedd_dE(3,3) = - normef * (gprime + hprime*nz*nz)

  if(normef.gt.1.e-8) then
     ! (Eddington tensor derivative with respect to Fx) * Ei
     ! remember 
     !    dchi/dFx = dchi/df df/dFx
     !    df/dFx = nx / (c E)
     !    dfx/dFx = 1 / (c E)
     !    dg/dFx = dg/df df/dFx
     !    d(h nx^2)/dFx = dh/df df/dFx nx^2 + h d(nx^2) / dFx
     !    d(nx^2)/dFx = 2 nx dnx/dfx dfx/dFx
     !    dnx / dfx = 1/f (1 - nx^2)
     !    dny / dFx = - nx ny / (f c E)
     ! so that
     ! d/dFx (g + h nx^2) = 1/(c E) [ dg/df nx + dh/df nx^3 + 2 h nx (1 - nx^2) / f ] 
     ! d/dFx (h nx ny) = 1/(c E) [ dh/df nx^2 ny + h ny (1 - 2 nx^2) / f ] 
     ! d/dFx (g + h ny ny) = 1/(c E) [ dg/df nx + dh/df nx ny^2 - 2 h nx ny^2 / f ] 
     Dedd_dF(1,1,1) = ( gprime*nx + hprime*nx*nx*nx + h*(-2.d0*nx*nx*nx /normef+2.d0*nx/normef) )/c
     Dedd_dF(1,2,1) = (             hprime*nx*nx*ny + h*(-2.d0*nx*nx*ny /normef+     ny/normef) )/c
     Dedd_dF(1,3,1) = (             hprime*nx*nx*nz + h*(-2.d0*nx*nx*nz /normef+     nz/normef) )/c 

     Dedd_dF(2,1,1) = (             hprime*nx*ny*nx + h*(-2.d0*nx*ny*nx /normef+     ny/normef) )/c
     Dedd_dF(2,2,1) = ( gprime*nx + hprime*nx*ny*ny + h*(-2.d0*nx*ny*ny /normef               ) )/c
     Dedd_dF(2,3,1) = (             hprime*nx*ny*nz + h*(-2.d0*nx*ny*nz /normef               ) )/c

     Dedd_dF(3,1,1) = (             hprime*nx*nz*nx + h*(-2.d0*nx*nz*nx /normef+     nz/normef) )/c
     Dedd_dF(3,2,1) = (             hprime*nx*nz*ny + h*(-2.d0*nx*nz*ny /normef               ) )/c
     Dedd_dF(3,3,1) = ( gprime*nx + hprime*nx*nz*nz + h*(-2.d0*nx*nz*nz /normef               ) )/c

     if(ndim.gt.1) then
        ! (Eddington tensor derivative with respect to Fy) * Ei
        Dedd_dF(1,1,2) = ( gprime*ny + hprime*ny*nx*nx + h*(-2.d0*ny*nx*nx /normef               ) )/c
        Dedd_dF(1,2,2) = (             hprime*ny*nx*ny + h*(-2.d0*ny*nx*ny /normef+     nx/normef) )/c
        Dedd_dF(1,3,2) = (             hprime*ny*nx*nz + h*(-2.d0*ny*nx*nz /normef               ) )/c 

        Dedd_dF(2,1,2) = (             hprime*ny*ny*nx + h*(-2.d0*ny*ny*nx /normef+     nx/normef) )/c
        Dedd_dF(2,2,2) = ( gprime*ny + hprime*ny*ny*ny + h*(-2.d0*ny*ny*ny /normef+2.d0*ny/normef) )/c
        Dedd_dF(2,3,2) = (             hprime*ny*ny*nz + h*(-2.d0*ny*ny*nz /normef+     nz/normef) )/c

        Dedd_dF(3,1,2) = (             hprime*ny*nz*nx + h*(-2.d0*ny*nz*nx /normef               ) )/c
        Dedd_dF(3,2,2) = (             hprime*ny*nz*ny + h*(-2.d0*ny*nz*ny /normef+     nz/normef) )/c
        Dedd_dF(3,3,2) = ( gprime*ny + hprime*ny*nz*nz + h*(-2.d0*ny*nz*nz /normef               ) )/c

        if(ndim.gt.2) then
           ! (Eddington tensor derivative with respect to Fz) * Ei
           Dedd_dF(1,1,3) = ( gprime*nz + hprime*nz*nx*nx + h*(-2.d0*nz*nx*nx /normef               ) )/c
           Dedd_dF(1,2,3) = (             hprime*nz*nx*ny + h*(-2.d0*nz*nx*ny /normef               ) )/c
           Dedd_dF(1,3,3) = (             hprime*nz*nx*nz + h*(-2.d0*nz*nx*nz /normef+     nx/normef) )/c 

           Dedd_dF(2,1,3) = (             hprime*nz*ny*nx + h*(-2.d0*nz*ny*nx /normef               ) )/c
           Dedd_dF(2,2,3) = ( gprime*nz + hprime*nz*ny*ny + h*(-2.d0*nz*ny*ny /normef               ) )/c
           Dedd_dF(2,3,3) = (             hprime*nz*ny*nz + h*(-2.d0*nz*ny*nz /normef+     ny/normef) )/c

           Dedd_dF(3,1,3) = (             hprime*nz*nz*nx + h*(-2.d0*nz*nz*nx /normef+     nx/normef) )/c
           Dedd_dF(3,2,3) = (             hprime*nz*nz*ny + h*(-2.d0*nz*nz*ny /normef+     ny/normef) )/c
           Dedd_dF(3,3,3) = ( gprime*nz + hprime*nz*nz*nz + h*(-2.d0*nz*nz*nz /normef+2.d0*nz/normef) )/c

        endif
     endif
  endif

  return

end subroutine cal_Dedd

!###########################################################
!###########################################################
!###########################################################
!###########################################################

! Write the results.
subroutine write_eigenvals
  
  implicit none
  integer :: i, j

  ! Write for C input
  !--------------------
  ! Note that the eigenvalues are sorted by increasing value,
  ! and we only need the min and the max. So we only write
  ! the indices 1 and 4.
  ! You can use this output to directly copy-paste eigenvalues
  ! into /src/rt/GEAR/rt_riemann_HLL_eigenvalues.h
  !
  ! write(*, "(A2)") "{ "
  ! do i = 0, n_points
  !   write(*, "(A4)") "  { "
  !   do j = 0, n_points
  !     write(*, "(A5, E12.6, A3, E12.6, A3)") "    {", eigenvals(i,j,1), " , ", eigenvals(i,j,4), " },"
  !   enddo
  !   write(*, "(A4)") "  },"
  ! enddo
  ! write(*, "(A2)") "};"

  ! Write for yourself
  !--------------------
  ! do i = 0, n_points
  !   do j = 0, n_points
  !     write(*, "(I3,x,I3,x,4(E12.6,x))") i, j, eigenvals(i, j, :)
  !   enddo
  ! enddo

  ! Write for python3 plotting
  !--------------------------
  open(unit=1, file="eigenvals.txt")
  write(1, "(A2,I5)") "# ", n_points + 1
  do i = 0, n_points
    do j = 0, n_points
      write(1, "(4(E12.6,x))") eigenvals(i, j, :)
    enddo
  enddo

end subroutine write_eigenvals

end program compute_RT_riemann_solver_eigenvalues
