real(kind = 8) function pes(sml_r, big_r, theta)
    implicit none

!   NOTE: on entry, sml_r and big_r have the same units used in the
!	wavefunctions, and theta has units of rad. Units of energy used
!	herein are also, by default, the units for all results.

    real(kind = 8), intent(in) :: sml_r, big_r, theta

    real(kind = 8) :: r12, r13, r23, x(3), y(3), z(3), der(3), inf

!
!   Set the He position (in the yz-plane):
!

    x(1) = 0.0d0
    y(1) = big_r*dcos(theta)
    z(1) = big_r*dsin(theta)

!
!   Set the H1 position (along the y-axis):
!

    x(2) = 0.0d0
    y(2) = sml_r/2.0d0
    z(2) = 0.0d0

!
!   Set the H2 position (along the y-axis):
!

    x(3) =  0.0d0
    y(3) = -y(2)
    z(3) =  0.0d0

!
!   Resolve the internuclear distances:
!

    r12 = atomic_distance(x(2), y(2), z(2), x(1), y(1), z(1))
    r13 = atomic_distance(x(3), y(3), z(3), x(1), y(1), z(1))
    r23 = atomic_distance(x(3), y(3), z(3), x(2), y(2), z(2))

!
!   NOTE: in this example both sml_r and big_r are in Angstrom
!   whereas the original PES expects them in atomic units.
!

    r12 = r12/0.52917721092d0 ! Angstrom to atomic units
    r13 = r13/0.52917721092d0
    r23 = r23/0.52917721092d0

!
!   Invoke the original PES:
!

    call fit3d(r12, r13, r23, pes, der)                ! 1 = He, 2 = 3 = H
    call fit3d(1000.0d0, 1000.0d0, 1000.0d0, inf, der) ! asymptotic energy

!
!   NOTE: for convenience we want to print all energies in wavenumber,
!   thus the PES value is returned in wavenumber.
!

    pes = (pes - inf)*219474.63137054d0 ! Atomic units to wavenumber
    return

    contains
    real(kind = 8) function atomic_distance(x1, y1, z1, x2, y2, z2)
        implicit none

        real(kind = 8), intent(in) :: x1, y1, z1, x2, y2, z2

        atomic_distance = sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)
    end function atomic_distance

end function pes

!
!     Original PES
!     Chemical Physics Letters 469 (2009) 26-30; doi: 10.1016/j.cplett.2008.12.035
!

      subroutine fit3d(r12,r13,r23,e,der)
      implicit real * 8 (a-h,o-z)
      dimension der(3)
      call diat12(r12,e12,d12)
      call diat12(r13,e13,d13)
      call diat23(r23,e23,d23)
      call triabb(r12,r13,r23,e123,der)
      e=e12+e13+e23+e123
      der(1)=d12+der(1)
      der(2)=d13+der(2)
      der(3)=d23+der(3)
      return
      end

!************************************************************************
      subroutine diat12(r,ener,der)
!************************************************************************
!*     This subroutine computes the energies of a diatomic potential
!*     fitted to    35 points
!*     rms =      0.00366602 kcal/mol
!*     emax =      0.01163009 kcal/mol
!************************************************************************
      implicit real*8 (a-h,o-z)
      dimension cf(  8)
      data cf(  1)/0.319875347303E+01/
      data cf(  2)/-.146567580366E-01/
      data cf(  3)/0.167399527625E+00/
      data cf(  4)/-.414208967596E+01/
      data cf(  5)/0.231491832608E+02/
      data cf(  6)/-.879814595809E+02/
      data cf(  7)/0.163457283605E+03/
      data cf(  8)/-.121233099595E+03/
      e0= 0.0000000E+00
      der=0.d0
      vex1= 0.8712143E+00
      vex2= 0.3534441E+01
      aux = 1.d0/r
      bux = dexp(-vex2*r)*aux
      cux = dexp(-vex1*r)
      ener=e0+cf(1)*bux
      dux=1.d0
      eux=r*cux
      do 1 i=2,  8
         der=der+(i-1)*cf(i)*dux
         dux=dux*eux
         ener=ener+cf(i)*dux
    1 continue
      der=der*(1.d0-vex1*r)*cux
      der=der-cf(1)*(vex2+aux)*bux
      return
      end

!************************************************************************
      subroutine diat23(r,ener,der)
!************************************************************************
!*     This subroutine computes the energies of a diatomic potential
!*     fitted to    37 points
!*     rms =      0.01436855 kcal/mol
!*     emax =      0.04507941 kcal/mol
!************************************************************************
      implicit real*8 (a-h,o-z)
      dimension cf(  8)
      data cf(  1)/0.102209827544E+01/
      data cf(  2)/-.721593775696E+00/
      data cf(  3)/0.184602439147E+01/
      data cf(  4)/-.223991681442E+01/
      data cf(  5)/-.193779665566E+02/
      data cf(  6)/0.111426109890E+03/
      data cf(  7)/-.234052990254E+03/
      data cf(  8)/0.180644082101E+03/
      e0= 0.0000000E+00
      der=0.d0
      vex1= 0.9880616E+00
      vex2= 0.1650299E+01
      aux = 1.d0/r
      bux = dexp(-vex2*r)*aux
      cux = dexp(-vex1*r)
      ener=e0+cf(1)*bux
      dux=1.d0
      eux=r*cux
      do 1 i=2,  8
         der=der+(i-1)*cf(i)*dux
         dux=dux*eux
         ener=ener+cf(i)*dux
    1 continue
      der=der*(1.d0-vex1*r)*cux
      der=der-cf(1)*(vex2+aux)*bux
      return
      end

!*************************************************************
      subroutine triabb(r12,r13,r23,ener,der)
!*************************************************************
!*     This subroutine computes the energies of a 3D PES
!*     for the ABB system class fitted to 1511 points
!*     rms =      0.14140340 kcal/mol
!*     emax =      1.01579755 kcal/mol
!*************************************************************
      implicit real*8(a-h,o-z)
      dimension i1(   78),i2(   78),i3(   78),i4(   78),cf(   78)
      dimension f12(0: 7),f13(0: 7),f23(0: 7)
      dimension der(3)
      data cf(  1)/0.1529045260044177E+01/
      data i1(  1)/ 1/,i2(  1)/ 1/,i3(  1)/ 0/,i4(  1)/ 2/
      data cf(  2)/0.2808607245907632E+01/
      data i1(  2)/ 1/,i2(  2)/ 0/,i3(  2)/ 1/,i4(  2)/ 1/
      data cf(  3)/0.1829411392554238E+02/
      data i1(  3)/ 1/,i2(  3)/ 1/,i3(  3)/ 1/,i4(  3)/ 1/
      data cf(  4)/0.1231966683878021E+02/
      data i1(  4)/ 2/,i2(  4)/ 1/,i3(  4)/ 0/,i4(  4)/ 2/
      data cf(  5)/-.4084703847036538E+02/
      data i1(  5)/ 2/,i2(  5)/ 0/,i3(  5)/ 1/,i4(  5)/ 2/
      data cf(  6)/-.4218399151150105E+02/
      data i1(  6)/ 0/,i2(  6)/ 2/,i3(  6)/ 1/,i4(  6)/ 2/
      data cf(  7)/0.6130817071836099E+03/
      data i1(  7)/ 2/,i2(  7)/ 1/,i3(  7)/ 1/,i4(  7)/ 2/
      data cf(  8)/-.8163596095489539E+03/
      data i1(  8)/ 1/,i2(  8)/ 2/,i3(  8)/ 1/,i4(  8)/ 1/
      data cf(  9)/-.5071101136073071E+02/
      data i1(  9)/ 2/,i2(  9)/ 2/,i3(  9)/ 0/,i4(  9)/ 2/
      data cf( 10)/0.4470022778867855E+03/
      data i1( 10)/ 2/,i2( 10)/ 0/,i3( 10)/ 2/,i4( 10)/ 1/
      data cf( 11)/-.2538360732067636E+02/
      data i1( 11)/ 3/,i2( 11)/ 1/,i3( 11)/ 0/,i4( 11)/ 2/
      data cf( 12)/0.2990848037524144E+03/
      data i1( 12)/ 3/,i2( 12)/ 0/,i3( 12)/ 1/,i4( 12)/ 2/
      data cf( 13)/0.4362731106893947E+03/
      data i1( 13)/ 0/,i2( 13)/ 3/,i3( 13)/ 1/,i4( 13)/ 2/
      data cf( 14)/-.1608985928952323E+04/
      data i1( 14)/ 2/,i2( 14)/ 2/,i3( 14)/ 1/,i4( 14)/ 2/
      data cf( 15)/-.4303319237725492E+04/
      data i1( 15)/ 2/,i2( 15)/ 1/,i3( 15)/ 2/,i4( 15)/ 1/
      data cf( 16)/-.2221718033436820E+04/
      data i1( 16)/ 3/,i2( 16)/ 1/,i3( 16)/ 1/,i4( 16)/ 2/
      data cf( 17)/0.5850081000463163E+04/
      data i1( 17)/ 1/,i2( 17)/ 3/,i3( 17)/ 1/,i4( 17)/ 1/
      data cf( 18)/0.1357484392614960E+04/
      data i1( 18)/ 3/,i2( 18)/ 2/,i3( 18)/ 0/,i4( 18)/ 2/
      data cf( 19)/-.2869125867638858E+04/
      data i1( 19)/ 3/,i2( 19)/ 0/,i3( 19)/ 2/,i4( 19)/ 2/
      data cf( 20)/-.1222599522273035E+04/
      data i1( 20)/ 0/,i2( 20)/ 3/,i3( 20)/ 2/,i4( 20)/ 2/
      data cf( 21)/-.4381487277752184E+03/
      data i1( 21)/ 4/,i2( 21)/ 1/,i3( 21)/ 0/,i4( 21)/ 2/
      data cf( 22)/-.1928118806358198E+04/
      data i1( 22)/ 4/,i2( 22)/ 0/,i3( 22)/ 1/,i4( 22)/ 2/
      data cf( 23)/-.2011488886074635E+04/
      data i1( 23)/ 0/,i2( 23)/ 4/,i3( 23)/ 1/,i4( 23)/ 2/
      data cf( 24)/0.1178682038193966E+05/
      data i1( 24)/ 2/,i2( 24)/ 2/,i3( 24)/ 2/,i4( 24)/ 1/
      data cf( 25)/0.1836389080999048E+05/
      data i1( 25)/ 3/,i2( 25)/ 2/,i3( 25)/ 1/,i4( 25)/ 2/
      data cf( 26)/0.1387412230221543E+05/
      data i1( 26)/ 3/,i2( 26)/ 1/,i3( 26)/ 2/,i4( 26)/ 2/
      data cf( 27)/-.9248027381820273E+04/
      data i1( 27)/ 1/,i2( 27)/ 3/,i3( 27)/ 2/,i4( 27)/ 2/
      data cf( 28)/-.5431269549656603E+04/
      data i1( 28)/ 3/,i2( 28)/ 3/,i3( 28)/ 0/,i4( 28)/ 2/
      data cf( 29)/0.1556197182655094E+05/
      data i1( 29)/ 3/,i2( 29)/ 0/,i3( 29)/ 3/,i4( 29)/ 1/
      data cf( 30)/-.2343093127542272E+04/
      data i1( 30)/ 4/,i2( 30)/ 1/,i3( 30)/ 1/,i4( 30)/ 2/
      data cf( 31)/-.1548344249834794E+05/
      data i1( 31)/ 1/,i2( 31)/ 4/,i3( 31)/ 1/,i4( 31)/ 1/
      data cf( 32)/-.1403512113776377E+04/
      data i1( 32)/ 4/,i2( 32)/ 2/,i3( 32)/ 0/,i4( 32)/ 2/
      data cf( 33)/0.9256314860653862E+04/
      data i1( 33)/ 4/,i2( 33)/ 0/,i3( 33)/ 2/,i4( 33)/ 2/
      data cf( 34)/0.8774018590169402E+04/
      data i1( 34)/ 0/,i2( 34)/ 4/,i3( 34)/ 2/,i4( 34)/ 2/
      data cf( 35)/0.1853340251583129E+04/
      data i1( 35)/ 5/,i2( 35)/ 1/,i3( 35)/ 0/,i4( 35)/ 2/
      data cf( 36)/0.8069823138998814E+04/
      data i1( 36)/ 5/,i2( 36)/ 0/,i3( 36)/ 1/,i4( 36)/ 2/
      data cf( 37)/0.4335123357729985E+04/
      data i1( 37)/ 0/,i2( 37)/ 5/,i3( 37)/ 1/,i4( 37)/ 2/
      data cf( 38)/-.1027837672581574E+05/
      data i1( 38)/ 3/,i2( 38)/ 2/,i3( 38)/ 2/,i4( 38)/ 2/
      data cf( 39)/-.5119978369381281E+05/
      data i1( 39)/ 2/,i2( 39)/ 3/,i3( 39)/ 2/,i4( 39)/ 1/
      data cf( 40)/-.3249654346467018E+05/
      data i1( 40)/ 3/,i2( 40)/ 3/,i3( 40)/ 1/,i4( 40)/ 2/
      data cf( 41)/-.3763209689619457E+05/
      data i1( 41)/ 3/,i2( 41)/ 1/,i3( 41)/ 3/,i4( 41)/ 1/
      data cf( 42)/-.1753931304040602E+05/
      data i1( 42)/ 4/,i2( 42)/ 2/,i3( 42)/ 1/,i4( 42)/ 2/
      data cf( 43)/-.1078215199912922E+05/
      data i1( 43)/ 4/,i2( 43)/ 1/,i3( 43)/ 2/,i4( 43)/ 2/
      data cf( 44)/0.4857428954696070E+05/
      data i1( 44)/ 1/,i2( 44)/ 4/,i3( 44)/ 2/,i4( 44)/ 2/
      data cf( 45)/0.1268841983642096E+05/
      data i1( 45)/ 4/,i2( 45)/ 3/,i3( 45)/ 0/,i4( 45)/ 2/
      data cf( 46)/-.3439906266667868E+05/
      data i1( 46)/ 4/,i2( 46)/ 0/,i3( 46)/ 3/,i4( 46)/ 2/
      data cf( 47)/0.1335508674905098E+04/
      data i1( 47)/ 0/,i2( 47)/ 4/,i3( 47)/ 3/,i4( 47)/ 2/
      data cf( 48)/0.8444128628183980E+04/
      data i1( 48)/ 5/,i2( 48)/ 1/,i3( 48)/ 1/,i4( 48)/ 2/
      data cf( 49)/0.1244676349023022E+05/
      data i1( 49)/ 1/,i2( 49)/ 5/,i3( 49)/ 1/,i4( 49)/ 1/
      data cf( 50)/-.3400118143957623E+04/
      data i1( 50)/ 5/,i2( 50)/ 2/,i3( 50)/ 0/,i4( 50)/ 2/
      data cf( 51)/-.1829633313535577E+05/
      data i1( 51)/ 5/,i2( 51)/ 0/,i3( 51)/ 2/,i4( 51)/ 2/
      data cf( 52)/-.1946752507903347E+05/
      data i1( 52)/ 0/,i2( 52)/ 5/,i3( 52)/ 2/,i4( 52)/ 2/
      data cf( 53)/-.2522586605989701E+04/
      data i1( 53)/ 6/,i2( 53)/ 1/,i3( 53)/ 0/,i4( 53)/ 2/
      data cf( 54)/-.1578145157501612E+05/
      data i1( 54)/ 6/,i2( 54)/ 0/,i3( 54)/ 1/,i4( 54)/ 2/
      data cf( 55)/-.3870024177831694E+04/
      data i1( 55)/ 0/,i2( 55)/ 6/,i3( 55)/ 1/,i4( 55)/ 2/
      data cf( 56)/-.1039922919110952E+06/
      data i1( 56)/ 3/,i2( 56)/ 3/,i3( 56)/ 2/,i4( 56)/ 2/
      data cf( 57)/0.8994670020057680E+05/
      data i1( 57)/ 3/,i2( 57)/ 2/,i3( 57)/ 3/,i4( 57)/ 1/
      data cf( 58)/0.1207873693639480E+06/
      data i1( 58)/ 4/,i2( 58)/ 2/,i3( 58)/ 2/,i4( 58)/ 2/
      data cf( 59)/0.1442186518579657E+06/
      data i1( 59)/ 2/,i2( 59)/ 4/,i3( 59)/ 2/,i4( 59)/ 1/
      data cf( 60)/0.3309059146593446E+05/
      data i1( 60)/ 4/,i2( 60)/ 3/,i3( 60)/ 1/,i4( 60)/ 2/
      data cf( 61)/-.1966208653695793E+05/
      data i1( 61)/ 4/,i2( 61)/ 1/,i3( 61)/ 3/,i4( 61)/ 2/
      data cf( 62)/0.1287939980085112E+05/
      data i1( 62)/ 1/,i2( 62)/ 4/,i3( 62)/ 3/,i4( 62)/ 2/
      data cf( 63)/-.1728183661036821E+05/
      data i1( 63)/ 4/,i2( 63)/ 4/,i3( 63)/ 0/,i4( 63)/ 2/
      data cf( 64)/0.4756468063498892E+05/
      data i1( 64)/ 4/,i2( 64)/ 0/,i3( 64)/ 4/,i4( 64)/ 1/
      data cf( 65)/-.1541909232570999E+05/
      data i1( 65)/ 5/,i2( 65)/ 2/,i3( 65)/ 1/,i4( 65)/ 2/
      data cf( 66)/-.5390709299722570E+05/
      data i1( 66)/ 5/,i2( 66)/ 1/,i3( 66)/ 2/,i4( 66)/ 2/
      data cf( 67)/-.6523309348251732E+05/
      data i1( 67)/ 1/,i2( 67)/ 5/,i3( 67)/ 2/,i4( 67)/ 2/
      data cf( 68)/0.1437739128177098E+04/
      data i1( 68)/ 5/,i2( 68)/ 3/,i3( 68)/ 0/,i4( 68)/ 2/
      data cf( 69)/0.3889049133794519E+05/
      data i1( 69)/ 5/,i2( 69)/ 0/,i3( 69)/ 3/,i4( 69)/ 2/
      data cf( 70)/0.1071996242727722E+05/
      data i1( 70)/ 0/,i2( 70)/ 5/,i3( 70)/ 3/,i4( 70)/ 2/
      data cf( 71)/0.7266271705665054E+04/
      data i1( 71)/ 6/,i2( 71)/ 1/,i3( 71)/ 1/,i4( 71)/ 2/
      data cf( 72)/0.4667504911442661E+04/
      data i1( 72)/ 1/,i2( 72)/ 6/,i3( 72)/ 1/,i4( 72)/ 1/
      data cf( 73)/0.1801284874493422E+04/
      data i1( 73)/ 6/,i2( 73)/ 2/,i3( 73)/ 0/,i4( 73)/ 2/
      data cf( 74)/0.2042197785608096E+05/
      data i1( 74)/ 6/,i2( 74)/ 0/,i3( 74)/ 2/,i4( 74)/ 2/
      data cf( 75)/0.1377349972736063E+05/
      data i1( 75)/ 0/,i2( 75)/ 6/,i3( 75)/ 2/,i4( 75)/ 2/
      data cf( 76)/0.1255357689466514E+04/
      data i1( 76)/ 7/,i2( 76)/ 1/,i3( 76)/ 0/,i4( 76)/ 2/
      data cf( 77)/0.1025949841141918E+05/
      data i1( 77)/ 7/,i2( 77)/ 0/,i3( 77)/ 1/,i4( 77)/ 2/
      data cf( 78)/0.7095229454036040E+03/
      data i1( 78)/ 0/,i2( 78)/ 7/,i3( 78)/ 1/,i4( 78)/ 2/
      vex1=0.1037770930000000E+01
      vex2=0.8809411423909257E+00
      f12(0)=1.d0
      f13(0)=1.d0
      f23(0)=1.d0
      bux12=r12*dexp(-vex1*r12)
      bux13=r13*dexp(-vex1*r13)
      bux23=r23*dexp(-vex2*r23)
      do 1 i=1, 7
         f12(i)=f12(i-1)*bux12
         f13(i)=f13(i-1)*bux13
         f23(i)=f23(i-1)*bux23
1     continue
      ener = 0.d0
      der12 = 0.d0
      der13 = 0.d0
      der23 = 0.d0
      do 2 l=1, 78
         if (i4(l).eq.1) then
            aux=f12(i1(l))*f13(i3(l))*f23(i2(l))
            dux12=i1(l)*f12(i1(l)-1)*f13(i3(l))*f23(i2(l))
            dux13=i3(l)*f12(i1(l))*f13(i3(l)-1)*f23(i2(l))
            dux23=i2(l)*f12(i1(l))*f13(i3(l))*f23(i2(l)-1)
         else
            aux1=f12(i1(l))*f13(i3(l))
            aux2=f12(i3(l))*f13(i1(l))
            aux=(aux1+aux2)*f23(i2(l))
            dux23=(aux1+aux2)*i2(l)*f23(i2(l)-1)
            dux1=i1(l)*f12(i1(l)-1)*f13(i3(l))
            dux2=i3(l)*f12(i3(l)-1)*f13(i1(l))
            dux12=(dux1+dux2)*f23(i2(l))
            dux1=i3(l)*f12(i1(l))*f13(i3(l)-1)
            dux2=i1(l)*f12(i3(l))*f13(i1(l)-1)
            dux13=(dux1+dux2)*f23(i2(l))
         endif
         ener=ener+cf(l)*aux
         der12=der12+cf(l)*dux12
         der13=der13+cf(l)*dux13
         der23=der23+cf(l)*dux23
    2 continue
      der(1)=der12*(1.d0-vex1*r12)*dexp(-vex1*r12)
      der(2)=der13*(1.d0-vex1*r13)*dexp(-vex1*r13)
      der(3)=der23*(1.d0-vex2*r23)*dexp(-vex2*r23)
      return
      end
