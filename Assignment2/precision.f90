    module precision
!  PRECISION - select the basic precision for the computations.
!
!  SYNOPSIS
!   use precision
!
!  DESCRIPTION
!   This module defines public names for floating-point KINDs.
!   It assumes IEEE arithmetic.
!
!  PUBLIC CONSTANTS DEFINED
!   SINGLE, DOUBLE - KIND names for floating-point types
!
!  REQUIRED DEPENDENCIES
!   None.
!
!  REVISION HISTORY
!    09/17/13 - First implementation.
!
!  PROGRAMMER
!   Your name, email@asu.edu.
!
!  PROGRAMMING NOTE
!   This is not documentation for the final code, but is given here as
!   explanation to get you started with Fortran.
!   Unless indicated otherwise, module variables are PUBLIC.
!   You may include the keyword PUBLIC to emphasize the public nature of any
!   declaration, at your option. In other words, the first statement below
!   is equivalant to    integer, parameter, private:: SINGLE = kind(1.0)
!   Fortran is not case sensitive, so SINGLE, Single, and single are
!   synonymous.  I recommnd, however, that you have some consistent and
!   easy method to distinguish named constants from other variables.
!
!  If you don't want to assume IEEE arithmetic, you can also say
!    integer, parameter, private:: SINGLE = selected_real_kind(6)
!    integer, parameter, private:: DOUBLE = selected_real_kind(15)
!
    integer, parameter:: SINGLE = kind(1.0)    ! IEEE single precision
    integer, parameter:: DOUBLE = kind(1.0d0)  ! IEEE double precision
    end module precision
