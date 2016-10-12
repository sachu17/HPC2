    module precision
!  PRECISION - select the basic precision for the computations.

!  DESCRIPTION
!   This module defines public names for floating-point KINDs.
!   It assumes IEEE arithmetic.
!
!  PUBLIC CONSTANTS DEFINED
!   SINGLE, DOUBLE - KIND names for floating-point types
!
!  REVISION HISTORY
!    09/17/13 - First implementation.
!
!  PROGRAMMER
!   Your name, sshivak8@asu.edu.
!              Sachin Shivakumar
!
    integer, parameter:: SINGLE = kind(1.0)    ! IEEE single precision
    integer, parameter:: DOUBLE = kind(1.0d0)  ! IEEE double precision
    end module precision
