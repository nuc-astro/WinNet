!
! file contains fermi function routines:
!
! function zfermim12 does a rational function fit for the order -1/2 integral
! function zfermi12 does a rational function fit for the order 1/2 integral
! function zfermi1 does a rational function fit for the order 1 integral
! function zfermi32 does a rational function fit for the order 3/2 integral
! function zfermi2 does a rational function fit for the order 2 integral
! function zfermi52 does a rational function fit for the order 5/2 integral
! function zfermi3 does a rational function fit for the order 3 integral
!
! function ifermim12 is a rational function fit for the inverse of order -1/2
! function ifermi12 is a rational function fit for the inverse of order 1/2
! function ifermi32 is a rational function fit for the inverse of order 3/2
! function ifermi52 is a rational function fit for the inverse of order 5/2





      real*8 function zfermim12(x)
      include 'implno.dek'
!
! this routine applies a rational function expansion to get the fermi function
! of order -1/2 evaluated at x. maximum error is 1.23d-12.
! reference: antia apjs 84,101 1993
!
! declare
      integer :: i,m1,k1,m2,k2
      real*8  :: x,an,a1(12),b1(12),a2(12),b2(12),rn,den,xx

! load the coefficients of the expansion
      data  an,m1,k1,m2,k2 /-0.5d0, 7, 7, 11, 11/
      data  (a1(i),i=1,8)/ 1.71446374704454d7,    3.88148302324068d7, &
                           3.16743385304962d7,    1.14587609192151d7, &
                           1.83696370756153d6,    1.14980998186874d5, &
                           1.98276889924768d3,    1.0d0/
      data  (b1(i),i=1,8)/ 9.67282587452899d6,    2.87386436731785d7, &
                           3.26070130734158d7,    1.77657027846367d7, &
                           4.81648022267831d6,    6.13709569333207d5, &
                           3.13595854332114d4,    4.35061725080755d2/
      data (a2(i),i=1,12)/-4.46620341924942d-15, -1.58654991146236d-12, &
                          -4.44467627042232d-10, -6.84738791621745d-8, &
                          -6.64932238528105d-6,  -3.69976170193942d-4, &
                          -1.12295393687006d-2,  -1.60926102124442d-1, &
                          -8.52408612877447d-1,  -7.45519953763928d-1, &
                           2.98435207466372d0,    1.0d0/
      data (b2(i),i=1,12)/-2.23310170962369d-15, -7.94193282071464d-13, &
                          -2.22564376956228d-10, -3.43299431079845d-8, &
                          -3.33919612678907d-6,  -1.86432212187088d-4, &
                          -5.69764436880529d-3,  -8.34904593067194d-2, &
                          -4.78770844009440d-1,  -4.99759250374148d-1, &
                           1.86795964993052d0,    4.16485970495288d-1/


      if (x .lt. 2.0d0) then
       xx = exp(x)
       rn = xx + a1(m1)
       do i=m1-1,1,-1
        rn = rn*xx + a1(i)
       enddo
       den = b1(k1+1)
       do i=k1,1,-1
        den = den*xx + b1(i)
       enddo
       zfermim12 = xx * rn/den
!
      else
       xx = 1.0d0/(x*x)
       rn = xx + a2(m2)
       do i=m2-1,1,-1
        rn = rn*xx + a2(i)
       enddo
       den = b2(k2+1)
       do i=k2,1,-1
        den = den*xx + b2(i)
       enddo
       zfermim12 = sqrt(x)*rn/den
      end if
      return
      end function zfermim12






      real*8 function zfermi12(x)
      include 'implno.dek'
!
! this routine applies a rational function expansion to get the fermi
! function of order 1/2 evaluated at x. maximum error is 5.47d-13.
! reference: antia apjs 84,101 1993
!
! declare
      integer :: i,m1,k1,m2,k2
      real*8  :: x,an,a1(12),b1(12),a2(12),b2(12),rn,den,xx


! load the coefficients of the expansion
      data  an,m1,k1,m2,k2 /0.5d0, 7, 7, 10, 11/
      data  (a1(i),i=1,8)/5.75834152995465d6,   1.30964880355883d7, &
                          1.07608632249013d7,   3.93536421893014d6, &
                          6.42493233715640d5,   4.16031909245777d4, &
                          7.77238678539648d2,   1.0d0/
      data  (b1(i),i=1,8)/6.49759261942269d6,   1.70750501625775d7, &
                          1.69288134856160d7,   7.95192647756086d6, &
                          1.83167424554505d6,   1.95155948326832d5, &
                          8.17922106644547d3,   9.02129136642157d1/
      data (a2(i),i=1,11)/4.85378381173415d-14, 1.64429113030738d-11, &
                          3.76794942277806d-9,  4.69233883900644d-7, &
                          3.40679845803144d-5,  1.32212995937796d-3, &
                          2.60768398973913d-2,  2.48653216266227d-1, &
                          1.08037861921488d0,   1.91247528779676d0, &
                          1.0d0/
      data (b2(i),i=1,12)/7.28067571760518d-14, 2.45745452167585d-11, &
                          5.62152894375277d-9,  6.96888634549649d-7, &
                          5.02360015186394d-5,  1.92040136756592d-3, &
                          3.66887808002874d-2,  3.24095226486468d-1, &
                          1.16434871200131d0,   1.34981244060549d0, &
                          2.01311836975930d-1, -2.14562434782759d-2/


      if (x .lt. 2.0d0) then
       xx = exp(x)
       rn = xx + a1(m1)
       do i=m1-1,1,-1
        rn = rn*xx + a1(i)
       enddo
       den = b1(k1+1)
       do i=k1,1,-1
        den = den*xx + b1(i)
       enddo
       zfermi12 = xx * rn/den

      else
       xx = 1.0d0/(x*x)
       rn = xx + a2(m2)
       do i=m2-1,1,-1
        rn = rn*xx + a2(i)
       enddo
       den = b2(k2+1)
       do i=k2,1,-1
        den = den*xx + b2(i)
       enddo
       zfermi12 = x*sqrt(x)*rn/den
      end if
      return
      end function zfermi12





      real*8 function zfermi1(x)
      include 'implno.dek'
!
! this routine applies a rational function expansion to get the fermi
! function of order 1 evaluated at x. maximum error is 1.0e-8.
! reference: antia  priv comm. 11sep94
!
! declare
      integer :: i,m1,k1,m2,k2
      real*8  :: x,an,a1(12),b1(12),a2(12),b2(12),rn,den,xx

! load the coefficients of the expansion
      data  an,m1,k1,m2,k2 /1.0, 7, 4, 9, 5/
      data  (a1(i),i=1,8)/-7.606458638543d7,  -1.143519707857d8, &
                          -5.167289383236d7,  -7.304766495775d6, &
                          -1.630563622280d5,   3.145920924780d3, &
                          -7.156354090495d1,   1.0d0/
      data  (b1(i),i=1,5)/-7.606458639561d7,  -1.333681162517d8, &
                          -7.656332234147d7,  -1.638081306504d7, &
                          -1.044683266663d6/
      data (a2(i),i=1,10)/-3.493105157219d-7, -5.628286279892d-5, &
                          -5.188757767899d-3, -2.097205947730d-1, &
                          -3.353243201574d0,  -1.682094530855d1, &
                          -2.042542575231d1,   3.551366939795d0, &
                          -2.400826804233d0,   1.0d0/
      data  (b2(i),i=1,6)/-6.986210315105d-7, -1.102673536040d-4, &
                          -1.001475250797d-2, -3.864923270059d-1, &
                          -5.435619477378d0,  -1.563274262745d1/


      if (x .lt. 2.0d0) then
       xx = exp(x)
       rn = xx + a1(m1)
       do i=m1-1,1,-1
        rn = rn*xx + a1(i)
       enddo
       den = b1(k1+1)
       do i=k1,1,-1
        den = den*xx + b1(i)
       enddo
       zfermi1 = xx * rn/den

      else
       xx = 1.0d0/(x*x)
       rn = xx + a2(m2)
       do i=m2-1,1,-1
        rn = rn*xx + a2(i)
       enddo
       den = b2(k2+1)
       do i=k2,1,-1
        den = den*xx + b2(i)
       enddo
       zfermi1 = x*x*rn/den
      end if
      return
      end function zfermi1





      real*8 function zfermi32(x)
      include 'implno.dek'
!
! this routine applies a rational function expansion to get the fermi
! function of order 3/2 evaluated at x. maximum error is 5.07d-13.
! reference: antia apjs 84,101 1993
!
! declare
      integer :: i,m1,k1,m2,k2
      real*8  :: x,an,a1(12),b1(12),a2(12),b2(12),rn,den,xx

! load the coefficients of the expansion
      data  an,m1,k1,m2,k2 /1.5d0, 6, 7, 9, 10/
      data  (a1(i),i=1,7)/4.32326386604283d4,   8.55472308218786d4, &
                          5.95275291210962d4,   1.77294861572005d4, &
                          2.21876607796460d3,   9.90562948053193d1, &
                          1.0d0/
      data  (b1(i),i=1,8)/3.25218725353467d4,   7.01022511904373d4, &
                          5.50859144223638d4,   1.95942074576400d4, &
                          3.20803912586318d3,   2.20853967067789d2, &
                          5.05580641737527d0,   1.99507945223266d-2/
      data (a2(i),i=1,10)/2.80452693148553d-13, 8.60096863656367d-11, &
                          1.62974620742993d-8,  1.63598843752050d-6, &
                          9.12915407846722d-5,  2.62988766922117d-3, &
                          3.85682997219346d-2,  2.78383256609605d-1, &
                          9.02250179334496d-1,  1.0d0/
      data (b2(i),i=1,11)/7.01131732871184d-13, 2.10699282897576d-10, &
                          3.94452010378723d-8,  3.84703231868724d-6, &
                          2.04569943213216d-4,  5.31999109566385d-3, &
                          6.39899717779153d-2,  3.14236143831882d-1, &
                          4.70252591891375d-1, -2.15540156936373d-2, &
                          2.34829436438087d-3/


      if (x .lt. 2.0d0) then
       xx = exp(x)
       rn = xx + a1(m1)
       do i=m1-1,1,-1
        rn = rn*xx + a1(i)
       enddo
       den = b1(k1+1)
       do i=k1,1,-1
        den = den*xx + b1(i)
       enddo
       zfermi32 = xx * rn/den

      else
       xx = 1.0d0/(x*x)
       rn = xx + a2(m2)
       do i=m2-1,1,-1
        rn = rn*xx + a2(i)
       enddo
       den = b2(k2+1)
       do i=k2,1,-1
        den = den*xx + b2(i)
       enddo
       zfermi32 = x*x*sqrt(x)*rn/den
      end if
      return
      end function zfermi32




      real*8 function zfermi2(x)
      include 'implno.dek'
!
! this routine applies a rational function expansion to get the fermi
! function of order 2 evaluated at x. maximum error is 1.0e-8.
! reference: antia  priv comm. 11sep94
!
! declare
      integer :: i,m1,k1,m2,k2
      real*8  :: x,an,a1(12),b1(12),a2(12),b2(12),rn,den,xx

! load the coefficients of the expansion
      data  an,m1,k1,m2,k2 /2.0, 7, 4, 5, 9/
      data  (a1(i),i=1,8)/-1.434885992395d8,  -2.001711155617d8, &
                          -8.507067153428d7,  -1.175118281976d7, &
                          -3.145120854293d5,   4.275771034579d3, &
                          -8.069902926891d1,   1.0d0/
      data  (b1(i),i=1,5)/-7.174429962316d7,  -1.090535948744d8, &
                          -5.350984486022d7,  -9.646265123816d6, &
                          -5.113415562845d5/
      data  (a2(i),i=1,6)/ 6.919705180051d-8,  1.134026972699d-5, &
                           7.967092675369d-4,  2.432500578301d-2, &
                           2.784751844942d-1,  1.0d0/
      data (b2(i),i=1,10)/ 2.075911553728d-7,  3.197196691324d-5, &
                           2.074576609543d-3,  5.250009686722d-2, &
                           3.171705130118d-1, -1.147237720706d-1, &
                           6.638430718056d-2, -1.356814647640d-2, &
                          -3.648576227388d-2,  3.621098757460d-2/


      if (x .lt. 2.0d0) then
       xx = exp(x)
       rn = xx + a1(m1)
       do i=m1-1,1,-1
        rn = rn*xx + a1(i)
       enddo
       den = b1(k1+1)
       do i=k1,1,-1
        den = den*xx + b1(i)
       enddo
       zfermi2 = xx * rn/den

      else
       xx = 1.0d0/(x*x)
       rn = xx + a2(m2)
       do i=m2-1,1,-1
        rn = rn*xx + a2(i)
       enddo
       den = b2(k2+1)
       do i=k2,1,-1
        den = den*xx + b2(i)
       enddo
       zfermi2 = x*x*x*rn/den
      end if
      return
      end function zfermi2






      real*8 function zfermi52(x)
      include 'implno.dek'
!
! this routine applies a rational function expansion to get the fermi
! function of order 5/2 evaluated at x. maximum error is 2.47d-13.
! reference: antia apjs 84,101 1993
!
! declare
      integer :: i,m1,k1,m2,k2
      real*8  :: x,an,a1(12),b1(12),a2(12),b2(12),rn,den,xx

! load the coefficients of the expansion
      data  an,m1,k1,m2,k2 /2.5d0, 6, 7, 10, 9/
      data  (a1(i),i=1,7)/6.61606300631656d4,   1.20132462801652d5, &
                          7.67255995316812d4,   2.10427138842443d4, &
                          2.44325236813275d3,   1.02589947781696d2, &
                          1.0d0/
      data  (b1(i),i=1,8)/1.99078071053871d4,   3.79076097261066d4, &
                          2.60117136841197d4,   7.97584657659364d3, &
                          1.10886130159658d3,   6.35483623268093d1, &
                          1.16951072617142d0,   3.31482978240026d-3/
      data (a2(i),i=1,11)/8.42667076131315d-12, 2.31618876821567d-9, &
                          3.54323824923987d-7,  2.77981736000034d-5, &
                          1.14008027400645d-3,  2.32779790773633d-2, &
                          2.39564845938301d-1,  1.24415366126179d0, &
                          3.18831203950106d0,   3.42040216997894d0, &
                          1.0d0/
      data (b2(i),i=1,10)/2.94933476646033d-11, 7.68215783076936d-9, &
                          1.12919616415947d-6,  8.09451165406274d-5, &
                          2.81111224925648d-3,  3.99937801931919d-2, &
                          2.27132567866839d-1,  5.31886045222680d-1, &
                          3.70866321410385d-1,  2.27326643192516d-2/


      if (x .lt. 2.0d0) then
       xx = exp(x)
       rn = xx + a1(m1)
       do i=m1-1,1,-1
        rn = rn*xx + a1(i)
       enddo
       den = b1(k1+1)
       do i=k1,1,-1
        den = den*xx + b1(i)
       enddo
       zfermi52 = xx * rn/den

      else
       xx = 1.0d0/(x*x)
       rn = xx + a2(m2)
       do i=m2-1,1,-1
        rn = rn*xx + a2(i)
       enddo
       den = b2(k2+1)
       do i=k2,1,-1
        den = den*xx + b2(i)
       enddo
       zfermi52 = x*x*x*sqrt(x)*rn/den
      end if
      return
      end function zfermi52





      real*8 function zfermi3(x)
      include 'implno.dek'
!
! this routine applies a rational function expansion to get the fermi
! function of order 3 evaluated at x. maximum error is 1.0e-8.
! reference: antia  priv comm. 11sep94
!
! declare
      integer  :: i,m1,k1,m2,k2
      real*8   :: x,an,a1(12),b1(12),a2(12),b2(12),rn,den,xx

! load the coefficients of the expansion
      data  an,m1,k1,m2,k2 /3.0, 4, 6, 7, 7/
      data  (a1(i),i=1,5)/ 6.317036716422d2,    7.514163924637d2, &
                           2.711961035750d2,    3.274540902317d1, &
                           1.0d0/
      data  (b1(i),i=1,7)/ 1.052839452797d2,    1.318163114785d2, &
                           5.213807524405d1,    7.500064111991d0, &
                           3.383020205492d-1,   2.342176749453d-3, &
                          -8.445226098359d-6/
      data  (a2(i),i=1,8)/ 1.360999428425d-8,   1.651419468084d-6, &
                           1.021455604288d-4,   3.041270709839d-3, &
                           4.584298418374d-2,   3.440523212512d-1, &
                           1.077505444383d0,    1.0d0/
      data  (b2(i),i=1,8)/ 5.443997714076d-8,   5.531075760054d-6, &
                           2.969285281294d-4,   6.052488134435d-3, &
                           5.041144894964d-2,   1.048282487684d-1, &
                           1.280969214096d-2,  -2.851555446444d-3/


      if (x .lt. 2.0d0) then
       xx = exp(x)
       rn = xx + a1(m1)
       do i=m1-1,1,-1
        rn = rn*xx + a1(i)
       enddo
       den = b1(k1+1)
       do i=k1,1,-1
        den = den*xx + b1(i)
       enddo
       zfermi3 = xx * rn/den

      else
       xx = 1.0d0/(x*x)
       rn = xx + a2(m2)
       do i=m2-1,1,-1
        rn = rn*xx + a2(i)
       enddo
       den = b2(k2+1)
       do i=k2,1,-1
        den = den*xx + b2(i)
       enddo
       zfermi3 = x*x*x*x*rn/den
      end if
      return
      end function zfermi3





      real*8 function ifermim12(f)
      include 'implno.dek'
!
! this routine applies a rational function expansion to get the inverse
! fermi function of order -1/2 when it is equal to f.
! maximum error is 3.03d-9.   reference: antia apjs 84,101 1993
!
! declare
      integer  :: i,m1,k1,m2,k2
      real*8   :: f,an,a1(12),b1(12),a2(12),b2(12),rn,den,ff

! load the coefficients of the expansion
      data  an,m1,k1,m2,k2 /-0.5d0, 5, 6, 6, 6/
      data  (a1(i),i=1,6)/-1.570044577033d4,   1.001958278442d4, &
                          -2.805343454951d3,   4.121170498099d2, &
                          -3.174780572961d1,   1.0d0/
      data  (b1(i),i=1,7)/-2.782831558471d4,   2.886114034012d4, &
                          -1.274243093149d4,   3.063252215963d3, &
                          -4.225615045074d2,   3.168918168284d1, &
                          -1.008561571363d0/
      data  (a2(i),i=1,7)/ 2.206779160034d-8,  -1.437701234283d-6, &
                           6.103116850636d-5,  -1.169411057416d-3, &
                           1.814141021608d-2,  -9.588603457639d-2, &
                           1.0d0/
      data  (b2(i),i=1,7)/ 8.827116613576d-8,  -5.750804196059d-6, &
                           2.429627688357d-4,  -4.601959491394d-3, &
                           6.932122275919d-2,  -3.217372489776d-1, &
                           3.124344749296d0/

      if (f .lt. 4.0d0) then
       rn = f + a1(m1)
       do i=m1-1,1,-1
        rn = rn*f + a1(i)
       enddo
       den = b1(k1+1)
       do i=k1,1,-1
        den = den*f + b1(i)
       enddo
       ifermim12 = log(f * rn/den)

      else
       ff = 1.0d0/f**(1.0d0/(1.0d0 + an))
       rn = ff + a2(m2)
       do i=m2-1,1,-1
        rn = rn*ff + a2(i)
       enddo
       den = b2(k2+1)
       do i=k2,1,-1
        den = den*ff + b2(i)
       enddo
       ifermim12 = rn/(den*ff)
      end if
      return
      end function ifermim12






      real*8 function ifermi12(f)
      include 'implno.dek'
!
! this routine applies a rational function expansion to get the inverse
! fermi function of order 1/2 when it is equal to f.
! maximum error is 4.19d-9.   reference: antia apjs 84,101 1993
!
! declare
      integer  :: i,m1,k1,m2,k2
      real*8   :: f,an,a1(12),b1(12),a2(12),b2(12),rn,den,ff

! load the coefficients of the expansion
      data  an,m1,k1,m2,k2 /0.5d0, 4, 3, 6, 5/
      data  (a1(i),i=1,5)/ 1.999266880833d4,   5.702479099336d3, &
                           6.610132843877d2,   3.818838129486d1, &
                           1.0d0/
      data  (b1(i),i=1,4)/ 1.771804140488d4,  -2.014785161019d3, &
                           9.130355392717d1,  -1.670718177489d0/
      data  (a2(i),i=1,7)/-1.277060388085d-2,  7.187946804945d-2, &
                          -4.262314235106d-1,  4.997559426872d-1, &
                          -1.285579118012d0,  -3.930805454272d-1, &
                           1.0d0/
      data  (b2(i),i=1,6)/-9.745794806288d-3,  5.485432756838d-2, &
                          -3.299466243260d-1,  4.077841975923d-1, &
                          -1.145531476975d0,  -6.067091689181d-2/


      if (f .lt. 4.0d0) then
       rn = f + a1(m1)
       do i=m1-1,1,-1
        rn = rn*f + a1(i)
       enddo
       den = b1(k1+1)
       do i=k1,1,-1
        den = den*f + b1(i)
       enddo
       ifermi12 = log(f * rn/den)

      else
       ff = 1.0d0/f**(1.0d0/(1.0d0 + an))
       rn = ff + a2(m2)
       do i=m2-1,1,-1
        rn = rn*ff + a2(i)
       enddo
       den = b2(k2+1)
       do i=k2,1,-1
        den = den*ff + b2(i)
       enddo
       ifermi12 = rn/(den*ff)
      end if
      return
      end function ifermi12






      real*8 function ifermi32(f)
      include 'implno.dek'
!
! this routine applies a rational function expansion to get the inverse
! fermi function of order 3/2 when it is equal to f.
! maximum error is 2.26d-9.   reference: antia apjs 84,101 1993
!
! declare
      integer  :: i,m1,k1,m2,k2
      real*8   :: f,an,a1(12),b1(12),a2(12),b2(12),rn,den,ff

! load the coefficients of the expansion
      data  an,m1,k1,m2,k2 /1.5d0, 3, 4, 6, 5/
      data  (a1(i),i=1,4)/ 1.715627994191d2,   1.125926232897d2, &
                           2.056296753055d1,   1.0d0/
      data  (b1(i),i=1,5)/ 2.280653583157d2,   1.193456203021d2, &
                           1.167743113540d1,  -3.226808804038d-1, &
                           3.519268762788d-3/
      data  (a2(i),i=1,7)/-6.321828169799d-3, -2.183147266896d-2, &
                          -1.057562799320d-1, -4.657944387545d-1, &
                          -5.951932864088d-1,  3.684471177100d-1, &
                           1.0d0/
      data  (b2(i),i=1,6)/-4.381942605018d-3, -1.513236504100d-2, &
                          -7.850001283886d-2, -3.407561772612d-1, &
                          -5.074812565486d-1, -1.387107009074d-1/


      if (f .lt. 4.0d0) then
       rn = f + a1(m1)
       do i=m1-1,1,-1
        rn = rn*f + a1(i)
       enddo
       den = b1(k1+1)
       do i=k1,1,-1
        den = den*f + b1(i)
       enddo
       ifermi32 = log(f * rn/den)

      else
       ff = 1.0d0/f**(1.0d0/(1.0d0 + an))
       rn = ff + a2(m2)
       do i=m2-1,1,-1
        rn = rn*ff + a2(i)
       enddo
       den = b2(k2+1)
       do i=k2,1,-1
        den = den*ff + b2(i)
       enddo
       ifermi32 = rn/(den*ff)
      end if
      return
      end function ifermi32





      real*8 function ifermi52(f)
      include 'implno.dek'
!
! this routine applies a rational function expansion to get the inverse
! fermi function of order 5/2 when it is equal to f.
! maximum error is 6.17d-9.   reference: antia apjs 84,101 1993
!
! declare
      integer  :: i,m1,k1,m2,k2
      real*8   :: f,an,a1(12),b1(12),a2(12),b2(12),rn,den,ff

! load the coefficients of the expansion
      data  an,m1,k1,m2,k2 /2.5d0, 2, 3, 6, 6/
      data  (a1(i),i=1,3)/ 2.138969250409d2,   3.539903493971d1, &
                           1.0d0/
      data  (b1(i),i=1,4)/ 7.108545512710d2,   9.873746988121d1, &
                           1.067755522895d0,  -1.182798726503d-2/
      data  (a2(i),i=1,7)/-3.312041011227d-2,  1.315763372315d-1, &
                          -4.820942898296d-1,  5.099038074944d-1, &
                           5.495613498630d-1, -1.498867562255d0, &
                           1.0d0/
      data  (b2(i),i=1,7)/-2.315515517515d-2,  9.198776585252d-2, &
                          -3.835879295548d-1,  5.415026856351d-1, &
                          -3.847241692193d-1,  3.739781456585d-2, &
                          -3.008504449098d-2/


      if (f .lt. 4.0d0) then
       rn = f + a1(m1)
       do i=m1-1,1,-1
        rn = rn*f + a1(i)
       enddo
       den = b1(k1+1)
       do i=k1,1,-1
        den = den*f + b1(i)
       enddo
       ifermi52 = log(f * rn/den)

      else
       ff = 1.0d0/f**(1.0d0/(1.0d0 + an))
       rn = ff + a2(m2)
       do i=m2-1,1,-1
        rn = rn*ff + a2(i)
       enddo
       den = b2(k2+1)
       do i=k2,1,-1
        den = den*ff + b2(i)
       enddo
       ifermi52 = rn/(den*ff)
      end if
      return
      end function ifermi52
