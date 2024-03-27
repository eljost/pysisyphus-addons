! [1] Rational Chebyshev Approximations for the Inverse of the Error Function
!     Blair, Edwards, Johnson, 1976

module mod_pa_erfinv

  use mod_pa_constants, only: dp

  implicit none

  private

  public :: erfinv

  ! Table 17, P00 to P06 for |x| <= 0.75
  real(dp), parameter :: P_0_075(7) = (/ &
       0.160304955844066229311d2, &
      -0.90784959262960326650d2, &
       0.18644914861620987391d3, &
      -0.16900142734642382420d3, &
       0.6545466284794487048d2, &
      -0.864213011587247794d1, &
       0.176058782139590d0 &
  /)

  ! Table 17, Q00 to Q06 for |x| <= 0.75
  real(dp), parameter :: Q_0_075(7) = (/ &
       0.147806470715138316110d2, &
      -0.91374167024260313936d2, &
       0.21015790486205317714d3, &
      -0.22210254121855132366d3, &
       0.10760453916055123830d3, &
      -0.206010730328265443d2, &
       0.1000000000000000d1 &
  /)

  ! Table 37, P00 to P07 for  0.75 <= |x| <= 0.9375
  real(dp), parameter :: P_075_09375(8) = (/ &
      -0.152389263440726128d-1, &
       0.3444556924136125216d0, &
      -0.29344398672542478687d1, &
       0.11763505705217627302d2, &
      -0.22655292823101104193d2, &
       0.19121334396580330163d2, &
      -0.5478927619598318769d1, &
       0.237516689024448000d0 &
  /)

  ! Table 37, Q00 to Q07 for  0.75 <= |x| <= 0.9375
  real(dp), parameter :: Q_075_09375(8) = (/ &
      -0.108465169602059954d-1, &
       0.2610628885843078511d0, &
      -0.24068318104393757995d1, &
       0.10695129973387014469d2, &
      -0.23716715521596581025d2, &
       0.24640158943917284883d2, &
      -0.10014376349783070835d2, &
       0.1000000000000000000d1 &
  /)

  ! Table 57, P00 to P08 for  0.9375 <= x < 1.0
  real(dp), parameter :: P_09375_1(9) = (/ &
       0.3100808562552958465d-4, &
       0.40974876030119399241d-2, &
       0.1214902662897276161664d0, &
       0.11091676946390282026357d1, &
       0.32283798556639239484387d1, &
       0.2881691815651598826211d1, &
       0.2047972087262996049729d1, &
       0.854592208197214832806d0, &
       0.3551095884622383139d-2 &
  /)

  ! Table 57, Q00 to Q08 for  0.9375 <= x < 1.0
  real(dp), parameter :: Q_09375_1(9) = (/ &
       0.3100809298564522487d-4, &
       0.40975286786639145041d-2, &
       0.1215907800748757143068d0, &
       0.11186271676316964210428d1, &
       0.34323639843052901805054d1, &
       0.4140284677116202150554d1, &
       0.4119797271272204121285d1, &
       0.2162961962641434560920d1, &
       0.10000000000000000000000d1 &
  /)

contains

  real(dp) function erfinv_inner(x, base0, P, Q)
    real(dp), intent(in) :: x, base0
    real(dp), intent(in) :: P(:), Q(:)
    real(dp) :: base
    real(dp) :: num, denom
    integer :: sizepq, i

    sizepq = size(P)
    base = base0
    num = 0
    denom = 0
    base = 1
    do i = 1, sizepq
      num = num + P(i) * base
      denom = denom + Q(i) * base
      base = base * base0
    end do
    erfinv_inner = x * num / denom
  end function erfinv_inner

  real(dp) function erfinv(x)
    ! Inverse error function for real arguments -1 < x < 1.
    real(dp), intent(in) :: x
    real(dp) :: abs_x, base

    ! TODO: return nan for invalid arguments
    abs_x = abs(x)
    if (abs_x <= 0.75d0) then
      base = x**2-0.5625d0
      erfinv = erfinv_inner(x, base, P_0_075, Q_0_075)
    elseif (abs_x > 0.75d0 .and. abs_x <= 0.9375d0) then
      base = x**2-0.87890625
      erfinv = erfinv_inner(x, base, P_075_09375, Q_075_09375)
    ! 0.9375 <= abs(x) < 1
    else
      base = 1.0d0 / sqrt(-log(1.0d0 - abs_x))
      erfinv = erfinv_inner(sign(1.0d0, x) / base, base, P_09375_1, Q_09375_1)
    end if

  end function erfinv
end module mod_pa_erfinv
