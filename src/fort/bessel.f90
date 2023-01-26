subroutine cik01 ( z, cbi0, cdi0, cbi1, cdi1, cbk0, cdk0, cbk1, cdk1 )

!*****************************************************************************80
!
!! CIK01: modified Bessel I0(z), I1(z), K0(z) and K1(z) for complex argument.
!
!  Discussion:
!
!    This procedure computes the modified Bessel functions I0(z), I1(z),
!    K0(z), K1(z), and their derivatives for a complex argument.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However,
!    they give permission to incorporate this routine into a user program
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    31 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, complex ( kind = 8 ) Z, the argument.
!
!    Output, complex ( kind = 8 ) CBI0, CDI0, CBI1, CDI1, CBK0, CDK0, CBK1,
!    CDK1, the values of I0(z), I0'(z), I1(z), I1'(z), K0(z), K0'(z), K1(z),
!    and K1'(z).
!
  implicit none

  real ( kind = 8 ), save, dimension ( 12 ) :: a = (/ &
    0.125D+00,           7.03125D-02,&
    7.32421875D-02,      1.1215209960938D-01,&
    2.2710800170898D-01, 5.7250142097473D-01,&
    1.7277275025845D+00, 6.0740420012735D+00,&
    2.4380529699556D+01, 1.1001714026925D+02,&
    5.5133589612202D+02, 3.0380905109224D+03 /)
  real ( kind = 8 ) a0
  real ( kind = 8 ), save, dimension ( 10 ) :: a1 = (/ &
    0.125D+00,            0.2109375D+00, &
    1.0986328125D+00,     1.1775970458984D+01, &
    2.1461706161499D+002, 5.9511522710323D+03, &
    2.3347645606175D+05,  1.2312234987631D+07, &
    8.401390346421D+08,   7.2031420482627D+10 /)
  real ( kind = 8 ), save, dimension ( 12 ) :: b = (/ &
   -0.375D+00,           -1.171875D-01, &
   -1.025390625D-01,     -1.4419555664063D-01, &
   -2.7757644653320D-01, -6.7659258842468D-01, &
   -1.9935317337513D+00, -6.8839142681099D+00, &
   -2.7248827311269D+01, -1.2159789187654D+02, &
   -6.0384407670507D+02, -3.3022722944809D+03 /)
  complex ( kind = 8 ) ca
  complex ( kind = 8 ) cb
  complex ( kind = 8 ) cbi0
  complex ( kind = 8 ) cbi1
  complex ( kind = 8 ) cbk0
  complex ( kind = 8 ) cbk1
  complex ( kind = 8 ) cdi0
  complex ( kind = 8 ) cdi1
  complex ( kind = 8 ) cdk0
  complex ( kind = 8 ) cdk1
  complex ( kind = 8 ) ci
  complex ( kind = 8 ) cr
  complex ( kind = 8 ) cs
  complex ( kind = 8 ) ct
  complex ( kind = 8 ) cw
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k0
  real ( kind = 8 ) pi
  real ( kind = 8 ) w0
  complex ( kind = 8 ) z
  complex ( kind = 8 ) z1
  complex ( kind = 8 ) z2
  complex ( kind = 8 ) zr
  complex ( kind = 8 ) zr2

  pi = 3.141592653589793D+00
  ci = cmplx ( 0.0D+00, 1.0D+00, kind = 8 )
  a0 = abs ( z )
  z2 = z * z
  z1 = z

  if ( a0 == 0.0D+00 ) then
    cbi0 = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
    cbi1 = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
    cdi0 = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
    cdi1 = cmplx ( 0.5D+00, 0.0D+00, kind = 8 )
    cbk0 = cmplx ( 1.0D+30, 0.0D+00, kind = 8 )
    cbk1 = cmplx ( 1.0D+30, 0.0D+00, kind = 8 )
    cdk0 = - cmplx ( 1.0D+30, 0.0D+00, kind = 8 )
    cdk1 = - cmplx ( 1.0D+30, 0.0D+00, kind = 8 )
    return
  end if

  if ( real ( z, kind = 8 ) < 0.0D+00 ) then
    z1 = -z
  end if

  if ( a0 <= 18.0D+00 ) then

    cbi0 = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
    cr = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
    do k = 1, 50
      cr = 0.25D+00 * cr * z2 / ( k * k )
      cbi0 = cbi0 + cr
      if ( abs ( cr / cbi0 ) < 1.0D-15 ) then
        exit
      end if
    end do

    cbi1 = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
    cr = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
    do k = 1, 50
      cr = 0.25D+00 * cr * z2 / ( k * ( k + 1 ) )
      cbi1 = cbi1 + cr
      if ( abs ( cr / cbi1 ) < 1.0D-15 ) then
        exit
      end if
    end do

    cbi1 = 0.5D+00 * z1 * cbi1

  else

    if ( a0 < 35.0D+00 ) then
      k0 = 12
    else if ( a0 < 50.0D+00 ) then
      k0 = 9
    else
      k0 = 7
    end if

    ca = exp ( z1 ) / sqrt ( 2.0D+00 * pi * z1 )
    cbi0 = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
    zr = 1.0D+00 / z1
    do k = 1, k0
      cbi0 = cbi0 + a(k) * zr ** k
    end do
    cbi0 = ca * cbi0
    cbi1 = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
    do k = 1, k0
      cbi1 = cbi1 + b(k) * zr ** k
    end do
    cbi1 = ca * cbi1

  end if

  if ( a0 <= 9.0D+00 ) then

    cs = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
    ct = - log ( 0.5D+00 * z1 ) - 0.5772156649015329D+00
    w0 = 0.0D+00
    cr = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
    do k = 1, 50
      w0 = w0 + 1.0D+00 / k
      cr = 0.25D+00 * cr / ( k * k ) * z2
      cs = cs + cr * ( w0 + ct )
      if ( abs ( ( cs - cw ) / cs ) < 1.0D-15 ) then
        exit
      end if
      cw = cs
    end do

    cbk0 = ct + cs

  else

    cb = 0.5D+00 / z1
    zr2 = 1.0D+00 / z2
    cbk0 = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
    do k = 1, 10
      cbk0 = cbk0 + a1(k) * zr2 ** k
    end do
    cbk0 = cb * cbk0 / cbi0

  end if

  cbk1 = ( 1.0D+00 / z1 - cbi1 * cbk0 ) / cbi0

  if ( real ( z, kind = 8 ) < 0.0D+00 ) then

    if ( imag ( z ) < 0.0D+00 ) then
      cbk0 = cbk0 + ci * pi * cbi0
      cbk1 = - cbk1 + ci * pi * cbi1
    else
      cbk0 = cbk0 - ci * pi * cbi0
      cbk1 = - cbk1 - ci * pi * cbi1
    end if

    cbi1 = - cbi1

  end if

  cdi0 = cbi1
  cdi1 = cbi0 - 1.0D+00 / z * cbi1
  cdk0 = - cbk1
  cdk1 = - cbk0 - 1.0D+00 / z * cbk1

  return
end
subroutine cjynb ( n, z, nm, cbj, cdj, cby, cdy )

!*****************************************************************************80
!
!! CJYNB: Bessel functions, derivatives, Jn(z) and Yn(z) of complex argument.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However,
!    they give permission to incorporate this routine into a user program
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    03 August 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of Jn(z) and Yn(z).
!
!    Input, complex ( kind = 8 ) Z, the argument of Jn(z) and Yn(z).
!
!    Output, integer ( kind = 4 ) NM, the highest order computed.
!
!    Output, complex ( kind = 8 ) CBJ(0:N), CDJ(0:N), CBY(0:N), CDY(0:N),
!    the values of Jn(z), Jn'(z), Yn(z), Yn'(z).
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ), save, dimension ( 4 ) :: a = (/ &
    -0.7031250000000000D-01, 0.1121520996093750D+00, &
    -0.5725014209747314D+00, 0.6074042001273483D+01 /)
  real ( kind = 8 ) a0
  real ( kind = 8 ), save, dimension ( 4 ) :: a1 = (/ &
    0.1171875000000000D+00,-0.1441955566406250D+00, &
    0.6765925884246826D+00,-0.6883914268109947D+01 /)
  real ( kind = 8 ), save, dimension ( 4 ) :: b = (/  &
    0.7324218750000000D-01,-0.2271080017089844D+00, &
    0.1727727502584457D+01,-0.2438052969955606D+02 /)
  real ( kind = 8 ), save, dimension ( 4 ) :: b1 = (/ &
   -0.1025390625000000D+00,0.2775764465332031D+00, &
   -0.1993531733751297D+01,0.2724882731126854D+02 /)
  complex ( kind = 8 ) cbj(0:n)
  complex ( kind = 8 ) cbj0
  complex ( kind = 8 ) cbj1
  complex ( kind = 8 ) cbjk
  complex ( kind = 8 ) cbs
  complex ( kind = 8 ) cby(0:n)
  complex ( kind = 8 ) cby0
  complex ( kind = 8 ) cby1
  complex ( kind = 8 ) cdj(0:n)
  complex ( kind = 8 ) cdy(0:n)
  complex ( kind = 8 ) ce
  complex ( kind = 8 ) cf
  complex ( kind = 8 ) cf1
  complex ( kind = 8 ) cf2
  complex ( kind = 8 ) cp0
  complex ( kind = 8 ) cp1
  complex ( kind = 8 ) cq0
  complex ( kind = 8 ) cq1
  complex ( kind = 8 ) cs0
  complex ( kind = 8 ) csu
  complex ( kind = 8 ) csv
  complex ( kind = 8 ) ct1
  complex ( kind = 8 ) ct2
  complex ( kind = 8 ) cu
  complex ( kind = 8 ) cyy
  real ( kind = 8 ) el
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) msta1
  integer ( kind = 4 ) msta2
  integer ( kind = 4 ) nm
  real ( kind = 8 ) pi
  real ( kind = 8 ) r2p
  real ( kind = 8 ) y0
  complex ( kind = 8 ) z

  el = 0.5772156649015329D+00
  pi = 3.141592653589793D+00
  r2p = 0.63661977236758D+00
  y0 = abs ( imag ( z ) )
  a0 = abs ( z )
  nm = n

  if ( a0 < 1.0D-100 ) then
    do k = 0, n
      cbj(k) = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
      cdj(k) = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
      cby(k) = - cmplx ( 1.0D+30, 0.0D+00, kind = 8 )
      cdy(k) = cmplx ( 1.0D+30, 0.0D+00, kind = 8 )
    end do
    cbj(0) = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
    cdj(1) = cmplx ( 0.5D+00, 0.0D+00, kind = 8 )
    return
  end if

  if ( a0 <= 300.0D+00 .or. 80 < n ) then

    if ( n == 0 ) then
      nm = 1
    end if
    m = msta1 ( a0, 200 )
    if ( m < nm ) then
      nm = m
    else
      m = msta2 ( a0, nm, 15 )
    end if

    cbs = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
    csu = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
    csv = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
    cf2 = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
    cf1 = cmplx ( 1.0D-30, 0.0D+00, kind = 8 )

    do k = m, 0, -1
      cf = 2.0D+00 * ( k + 1.0D+00 ) / z * cf1 - cf2
      if ( k <= nm ) then
        cbj(k) = cf
      end if
      if ( k == 2 * int ( k / 2 ) .and. k .ne. 0 ) then
        if ( y0 <= 1.0D+00 ) then
          cbs = cbs + 2.0D+00 * cf
        else
          cbs = cbs + ( -1.0D+00 ) ** ( k / 2 ) * 2.0D+00 * cf
        end if
        csu = csu + ( -1.0D+00 ) ** ( k / 2 ) * cf / k
      else if ( 1 < k ) then
        csv = csv + ( -1.0D+00 ) ** ( k / 2 ) * k / ( k * k - 1.0D+00 ) * cf
      end if
      cf2 = cf1
      cf1 = cf
    end do

    if ( y0 <= 1.0D+00 ) then
      cs0 = cbs + cf
    else
      cs0 = ( cbs + cf ) / cos ( z )
    end if

    do k = 0, nm
      cbj(k) = cbj(k) / cs0
    end do

    ce = log ( z / 2.0D+00 ) + el
    cby(0) = r2p * ( ce * cbj(0) - 4.0D+00 * csu / cs0 )
    cby(1) = r2p * ( - cbj(0) / z + ( ce - 1.0D+00 ) * cbj(1) &
      - 4.0D+00 * csv / cs0 )

  else

    ct1 = z - 0.25D+00 * pi
    cp0 = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
    do k = 1, 4
      cp0 = cp0 + a(k) * z ** ( - 2 * k )
    end do
    cq0 = -0.125D+00 / z
    do k = 1, 4
      cq0 = cq0 + b(k) * z ** ( - 2 * k - 1 )
    end do
    cu = sqrt ( r2p / z )
    cbj0 = cu * ( cp0 * cos ( ct1 ) - cq0 * sin ( ct1 ) )
    cby0 = cu * ( cp0 * sin ( ct1 ) + cq0 * cos ( ct1 ) )
    cbj(0) = cbj0
    cby(0) = cby0
    ct2 = z - 0.75D+00 * pi
    cp1 = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
    do k = 1, 4
      cp1 = cp1 + a1(k) * z ** ( - 2 * k )
    end do
    cq1 = 0.375D+00 / z
    do k = 1, 4
      cq1 = cq1 + b1(k) * z ** ( - 2 * k - 1 )
    end do
    cbj1 = cu * ( cp1 * cos ( ct2 ) - cq1 * sin ( ct2 ) )
    cby1 = cu * ( cp1 * sin ( ct2 ) + cq1 * cos ( ct2 ) )
    cbj(1) = cbj1
    cby(1) = cby1
    do k = 2, nm
      cbjk = 2.0D+00 * ( k - 1.0D+00 ) / z * cbj1 - cbj0
      cbj(k) = cbjk
      cbj0 = cbj1
      cbj1 = cbjk
    end do
  end if

  cdj(0) = -cbj(1)
  do k = 1, nm
    cdj(k) = cbj(k-1) - k / z * cbj(k)
  end do

  if ( 1.0D+00 < abs ( cbj(0) ) ) then
    cby(1) = ( cbj(1) * cby(0) - 2.0D+00 / ( pi * z ) ) / cbj(0)
  end if

  do k = 2, nm
    if ( abs ( cbj(k-2) ) <= abs ( cbj(k-1) ) ) then
      cyy = ( cbj(k) * cby(k-1) - 2.0D+00 / ( pi * z ) ) / cbj(k-1)
    else
      cyy = ( cbj(k) * cby(k-2) - 4.0D+00 * ( k - 1.0D+00 ) &
        / ( pi * z * z ) ) / cbj(k-2)
    end if
    cby(k) = cyy
  end do

  cdy(0) = -cby(1)
  do k = 1, nm
    cdy(k) = cby(k-1) - k / z * cby(k)
  end do

  return
end
function msta1 ( x, mp )

!*****************************************************************************80
!
!! MSTA1 determines a backward recurrence starting point for Jn(x).
!
!  Discussion:
!
!    This procedure determines the starting point for backward
!    recurrence such that the magnitude of
!    Jn(x) at that point is about 10^(-MP).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However,
!    they give permission to incorporate this routine into a user program
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    08 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Input, integer ( kind = 4 ) MP, the negative logarithm of the
!    desired magnitude.
!
!    Output, integer ( kind = 4 ) MSTA1, the starting point.
!
  implicit none

  real ( kind = 8 ) a0
  real ( kind = 8 ) envj
  real ( kind = 8 ) f
  real ( kind = 8 ) f0
  real ( kind = 8 ) f1
  integer ( kind = 4 ) it
  integer ( kind = 4 ) mp
  integer ( kind = 4 ) msta1
  integer ( kind = 4 ) n0
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) nn
  real ( kind = 8 ) x

  a0 = abs ( x )
  n0 = int ( 1.1D+00 * a0 ) + 1
  f0 = envj ( n0, a0 ) - mp
  n1 = n0 + 5
  f1 = envj ( n1, a0 ) - mp
  do it = 1, 20
    nn = n1 - int ( real ( n1 - n0, kind = 8 ) / ( 1.0D+00 - f0 / f1 ) )
    f = envj ( nn, a0 ) - mp
    if ( abs ( nn - n1 ) < 1 ) then
      exit
    end if
    n0 = n1
    f0 = f1
    n1 = nn
    f1 = f
  end do

  msta1 = nn

  return
end
function msta2 ( x, n, mp )

!*****************************************************************************80
!
!! MSTA2 determines a backward recurrence starting point for Jn(x).
!
!  Discussion:
!
!    This procedure determines the starting point for a backward
!    recurrence such that all Jn(x) has MP significant digits.
!
!    Jianming Jin supplied a modification to this code on 12 January 2016.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However,
!    they give permission to incorporate this routine into a user program
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    14 January 2016
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of Jn(x).
!
!    Input, integer ( kind = 4 ) N, the order of Jn(x).
!
!    Input, integer ( kind = 4 ) MP, the number of significant digits.
!
!    Output, integer ( kind = 4 ) MSTA2, the starting point.
!
  implicit none

  real ( kind = 8 ) a0
  real ( kind = 8 ) ejn
  real ( kind = 8 ) envj
  real ( kind = 8 ) f
  real ( kind = 8 ) f0
  real ( kind = 8 ) f1
  real ( kind = 8 ) hmp
  integer ( kind = 4 ) it
  integer ( kind = 4 ) mp
  integer ( kind = 4 ) msta2
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n0
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) nn
  real ( kind = 8 ) obj
  real ( kind = 8 ) x

  a0 = abs ( x )
  hmp = 0.5D+00 * mp
  ejn = envj ( n, a0 )

  if ( ejn <= hmp ) then
    obj = mp
!
!  Original code:
!
!   n0 = int ( 1.1D+00 * a0 )
!
!  Updated code:
!
    n0 = int ( 1.1D+00 * a0 ) + 1
  else
    obj = hmp + ejn
    n0 = n
  end if

  f0 = envj ( n0, a0 ) - obj
  n1 = n0 + 5
  f1 = envj ( n1, a0 ) - obj

  do it = 1, 20
    nn = n1 - int ( real ( n1 - n0, kind = 8 ) / ( 1.0D+00 - f0 / f1 ) )
    f = envj ( nn, a0 ) - obj
    if ( abs ( nn - n1 ) < 1 ) then
      exit
    end if
    n0 = n1
    f0 = f1
    n1 = nn
    f1 = f
  end do

  msta2 = nn + 10

  return
end
function envj ( n, x )

!*****************************************************************************80
!
!! ENVJ is a utility function used by MSTA1 and MSTA2.
!
!  Discussion:
!
!    ENVJ estimates -log(Jn(x)) from the estimate
!    Jn(x) approx 1/sqrt(2*pi*n) * ( e*x/(2*n))^n
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However,
!    they give permission to incorporate this routine into a user program
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    14 January 2016
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!    Modifications suggested by Vincent Lafage, 11 January 2016.
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the Bessel function.
!
!    Input, real ( kind = 8 ) X, the absolute value of the argument.
!
!    Output, real ( kind = 8 ) ENVJ, the value.
!
  implicit none

  real ( kind = 8 ) envj
  real ( kind = 8 ) logten
  integer ( kind = 4 ) n
  real ( kind = 8 ) n_r8
  real ( kind = 8 ) r8_gamma_log
  real ( kind = 8 ) x
!
!  Original code
!
  if ( .true. ) then

    envj = 0.5D+00 * log10 ( 6.28D+00 * n ) &
      - n * log10 ( 1.36D+00 * x / n )
!
!  Modification suggested by Vincent Lafage.
!
  else

    n_r8 = real ( n, kind = 8 )
    logten = log ( 10.0D+00 )
    envj = r8_gamma_log ( n_r8 + 1.0D+00 ) / logten - n_r8 * log10 ( x )

  end if

  return
end
