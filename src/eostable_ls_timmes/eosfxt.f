
c..routine eosfxt computes a stellar eos assuming complete ionozation
c..routine etages makes a good guess for the chemical potential
c..routine xneroot gets the thermodynamics of the electrons and positrons







      subroutine eosfxt
      include 'implno.dek'
      include 'const.dek'
      include 'vector_eos.dek'

c..given a temperature temp [K], density den [g/cm**3], and a
composition
c..characterized by abar (average weight) and zbar (average charge),
c..this routine returns all the other thermodynamic quantities.

c..of interest is the pressure [erg/cm**3], specific thermal energy [erg/gr],
c..the entropy [erg/g/K], with their derivatives with respect to temperature,
c..density, abar, and zbar.

c..other quantites such the normalized chemical potential eta (plus its
c..derivatives), number density of electrons and positron pair (along
c..with their derivatives), adiabatic indices, specific heats, and
c..relativistically correct sound speed are also returned.

c..this routine assumes planckian photons, an ideal gas of ions,
c..and an electron-positron gas with an arbitrary degree of relativity
c..and degeneracy. the full fermi-dirac integrals and their derivatives
c..with respect to eta and beta are computed to machine precision, and
c..all other derivatives are analytic.

c..references: cox & giuli (c&g) chapter 24,
c..            timmes & arnett, apj supp. 125, 277, 1999
c..            timmes & swesty, apj supp. 126, 501, 2000




c..a dictionary of terms used:
c..this routine has now been pipelined.
c..all the input and output variables are in the file vector_eos.dek.
c..the vector name is the scaler name appended with an "_row",
c..for example, temp_row(i), den_row(i), and so on.



c..input:
c..temp     = temperature
c..den      = density
c..abar     = average number of nucleons per nuclei
c..zbar     = average number of protons per nuclei


c..output:

c..pres     = total pressure
c..dpresdd  = derivative of total pressure with respect to density
c..dpresdt  = derivative of total pressure with respect to temperature
c..dpresda  = derivative of total pressure with respect to abar
c..dpresdz  = derivative of total pressure with respect to zbar

c..ener     = total internal energy
c..denerdd  = derivative of total energy with respect to density
c..denerdt  = derivative of total energy with respect to temperature
c..denerda  = derivative of total energy with respect to abar
c..denerdz  = derivative of total energy with respect to zbar

c..entr     = total entropy
c..dentrdd  = derivative of total entropy with respect to density
c..dentrdt  = derivative of total entropy with respect to temperature
c..dentrda  = derivative of total entropy with respect to abar
c..dentrdz  = derivative of total entropy with respect to zbar



c..prad     = radiation pressure
c..dpraddd  = derivative of the radiation pressure with density
c..dpraddt  = derivative of the radiation pressure with temperature
c..dpradda  = derivative of the radiation pressure with abar
c..dpraddz  = derivative of the radiation pressure with zbar

c..erad     = radiation energy
c..deraddd  = derivative of the radiation energy with density
c..deraddt  = derivative of the radiation energy with temperature
c..deradda  = derivative of the radiation energy with abar
c..deraddz  = derivative of the radiation energy with zbar

c..srad     = radiation entropy
c..dsraddd  = derivative of the radiation entropy with density
c..dsraddt  = derivative of the radiation entropy with temperature
c..dsradda  = derivative of the radiation entropy with abar
c..dsraddz  = derivative of the radiation entropy with zbar

c..radmult  = radiation multiplier (useful for turning radiation off/on)




c..xni      = number density of ions
c..dxnidd   = derivative of the ion number density with density
c..dxnidt   = derivative of the ion number density with temperature
c..dxnida   = derivative of the ion number density with abar
c..dxnidz   = derivative of the ion number density with zbar

c..pion     = ion pressure
c..dpiondd  = derivative of the ion pressure with density
c..dpiondt  = derivative of the ion pressure with temperature
c..dpionda  = derivative of the ion pressure with abar
c..dpiondz  = derivative of the ion pressure with zbar

c..eion     = ion energy
c..deiondd  = derivative of the ion energy with density
c..deiondt  = derivative of the ion energy with temperature
c..deionda  = derivative of the ion energy with abar
c..deiondz  = derivative of the ion energy with zbar

c..sion     = ion entropy
c..dsiondd  = derivative of the ion entropy with density
c..dsiondt  = derivative of the ion entropy with temperature
c..dsionda  = derivative of the ion entropy with abar
c..dsiondz  = derivative of the ion entropy with zbar

c..ionmult  = ion multiplier (useful for turning ions off/on)


c..etaele   = electron chemical potential
c..detadd   = derivative of the electron chem potential with density
c..detadt   = derivative of the electron chem potential with temperature
c..detada   = derivative of the electron chem potential with abar
c..detadz   = derivative of the electron chem potential with zbar

c..etapos   = positron degeneracy parameter




c..xne       = number density of electrons
c..dxnedd    = derivative of the electron number density with density
c..dxnedt    = derivative of the electron number density with temperature
c..dxneda    = derivative of the electron number density with abar
c..dxnedz    = derivative of the electron number density with zbar

c..xnefer    = fermi integral electron number density
c..dxneferdd = derivative of the fermi electron number density with density
c..dxneferdt = derivative of the fermi electron number density with temperature
c..dxneferda = derivative of the fermi electron number density with abar
c..dxneferdz = derivative of the fermi electron number density with zbar

c..xnpfer    = fermi integral positron number density
c..dxnpferdd = derivative of the fermi positron number density with density
c..dxnpferdt = derivative of the fermi positron number density with temperature
c..dxnpferda = derivative of the fermi positron number density with abar
c..dxnpferdz = derivative of the fermi positron number density with zbar

c..pele      = electron pressure
c..dpeledd   = derivative of the electron pressure with density
c..dpeledt   = derivative of the electron pressure with temperature
c..dpeleda   = derivative of the electron pressure with abar
c..dpeledz   = derivative of the electron pressure with zbar

c..eele     = electron energy
c..deeledd   = derivative of the electron energy with density
c..deeledt   = derivative of the electron energy with temperature
c..deeleda   = derivative of the electron energy with abar
c..deeledz   = derivative of the electron energy with zbar

c..sele     = electron entropy
c..dseledd   = derivative of the electron entropy with density
c..dseledt   = derivative of the electron entropy with temperature
c..dseleda   = derivative of the electron entropy with abar
c..dseledz   = derivative of the electron entropy with zbar


c..ppos     = positron pressure
c..dpposdd   = derivative of the positron pressure with density
c..dpposdt   = derivative of the positron pressure with temperature
c..dpposda   = derivative of the positron pressure with abar
c..dpposdz   = derivative of the positron pressure with zbar

c..epos     = electron energy
c..deposdd   = derivative of the positron energy with density
c..deposdt   = derivative of the positron energy with temperature
c..deposda   = derivative of the positron energy with abar
c..deposdz   = derivative of the positron energy with zbar

c..spos     = electron entropy
c..dsposdd   = derivative of the positron entropy with density
c..dsposdt   = derivative of the positron entropy with temperature
c..dsposda   = derivative of the positron entropy with abar
c..dsposdz   = derivative of the positron entropy with zbar

c..pep      = electron + positron pressure
c..dpepdd   = derivative of the electron+positron pressure with density
c..dpepdt   = derivative of the electron+positron pressure with temperature
c..dpepda   = derivative of the electron+positron pressure with abar
c..dpepdz   = derivative of the electron+positron pressure with zbar

c..eep      = electron + ositron energy
c..deepdd   = derivative of the electron+positron energy with density
c..deepdt   = derivative of the electron+positron energy with temperature
c..deepda   = derivative of the electron+positron energy with abar
c..deepdz   = derivative of the electron+positron energy with zbar

c..sep      = electron + positron entropy
c..dsepdd   = derivative of the electron+positron entropy with density
c..dsepdt   = derivative of the electron+positron entropy with temperature
c..dsepda   = derivative of the electron+positron entropy with abar
c..dsepdz   = derivative of the electron+positron entropy with zbar

c..elemult  = electron multiplier (useful for turning e-e+ off/on)


c..eip      = ionization potential ennergy
c..deipdd   = derivative of ionization energy with density
c..deipdt   = derivative of ionization energy with temperature
c..deipda   = derivative of ionization energy with abar
c..deipdz   = derivative of ionization energy with zbar


c..sip      = ionization potential ennergy
c..dsipdd   = derivative of ionization energy with density
c..dsipdt   = derivative of ionization energy with temperature
c..dsipda   = derivative of ionization energy with abar
c..dsipdz   = derivative of ionization energy with zbar

c..potmult  = ionization energy multiplier (useful for turning off ionization additions)



c..pcoul    = coulomb pressure correction
c..coulmult = coulomb component multiplier
c..dpcouldd = derivative of the coulomb pressure with density
c..dpcouldt = derivative of the coulomb pressure with temperature
c..dpcoulda = derivative of the coulomb pressure with abar
c..dpcouldz = derivative of the coulomb pressure with zbar

c..ecoul    = coulomb energy correction
c..decouldd = derivative of the coulomb energy with density
c..decouldt = derivative of the coulomb energy with temperature
c..decoulda = derivative of the coulomb energy with abar
c..decouldz = derivative of the coulomb energy with zbar

c..scoul    = coulomb entropy correction
c..dscouldd = derivative of the coulomb entropy with density
c..dscouldt = derivative of the coulomb entropy with temperature
c..dscoulda = derivative of the coulomb entropy with abar
c..dscouldz = derivative of the coulomb entropy with zbar


c..kt       = kerg * temperature
c..beta     = dimensionless ratio of kerg*temp/me*c^2

c..chit     = temperature exponent in the pressure equation of state
c..chid     = density exponent in the pressure equation of state
c..cv       = specific heat at constant volume
c..cp       = specific heat at constant pressure
c..gam1     = first adiabatic exponent
c..gam2     = second adiabatic exponent
c..gam3     = third adiabatic exponent
c..nabad    = adiabatic gradient
c..sound    = relativistically correct adiabatic sound speed
c..plasg    = ratio of electrostatic to thermal energy


c..dse      = thermodynamic consistency check de/dt = t*ds/dt
c..dpe      = thermodynamic consistency check p = d**2 de/dd + t*dpdt
c..dsp      = thermodynamic consistency check dp/dt = - d**2 ds/dd





c..declare the input
      double precision temp,den,zbar,abar


c..declare everything else
c..totals
      double precision pres,ener,entr,
     1                 dpresdd,dpresdt,dpresda,dpresdz,
     2                 denerdd,denerdt,denerda,denerdz,
     3                 dentrdd,dentrdt,dentrda,dentrdz


c..radiation
      integer          radmult
      double precision prad,erad,srad,
     1                 dpraddd,dpraddt,dpradda,dpraddz,
     2                 deraddd,deraddt,deradda,deraddz,
     3                 dsraddd,dsraddt,dsradda,dsraddz


c..ions
      integer          ionmult
      double precision pion,eion,sion,xni,etaion,
     1                 dpiondd,dpiondt,dpionda,dpiondz,
     2                 deiondd,deiondt,deionda,deiondz,
     3                 dsiondd,dsiondt,dsionda,dsiondz,
     4                 dxnidd,dxnidt,dxnida,dxnidz,
     5                 detaiondd,detaiondt,detaionda,detaiondz


c..electron-positrons
      integer          elemult
      double precision etaele,detadd,detadt,detada,detadz,
     1                 etapos,zeff,
     2                 xne,dxnedd,dxnedt,dxneda,dxnedz,
     3                 xnefer,dxneferdd,dxneferdt,dxneferda,dxneferdz,
     4                 xnpfer,dxnpferdd,dxnpferdt,dxnpferda,dxnpferdz,
     5                 pele,dpeledd,dpeledt,dpeleda,dpeledz,
     6                 ppos,dpposdd,dpposdt,dpposda,dpposdz,
     7                 pep,dpepdd,dpepdt,dpepda,dpepdz,
     8                 eele,deeledd,deeledt,deeleda,deeledz,
     9                 epos,deposdd,deposdt,deposda,deposdz,
     1                 eep,deepdd,deepdt,deepda,deepdz,
     2                 sele,dseledd,dseledt,dseleda,dseledz,
     3                 spos,dsposdd,dsposdt,dsposda,dsposdz,
     4                 sep,dsepdd,dsepdt,dsepda,dsepdz


c..ionization contributions
      integer          ionized,potmult
      double precision eip,deipdd,deipdt,deipda,deipdz,
     1                 sip,dsipdd,dsipdt,dsipda,dsipdz


c..coulomb corrections
      integer          coulmult
      double precision plasg,
     1                 pcoul,dpcouldd,dpcouldt,dpcoulda,dpcouldz,
     2                 ecoul,decouldd,decouldt,decoulda,decouldz,
     3                 scoul,dscouldd,dscouldt,dscoulda,dscouldz



c..various physical quantities based on derivatives
      double precision chit,chid,cv,cp,gam1,gam2,gam3,nabad,sound


c..for the maxwell relations
      double precision dse,dpe,dsp


c..miscelaneous local variables
      integer          i,j,niter,mode
      double precision kt,ktinv,x,y,z,xx,yy,zz,ages,agesav,agesnew,
     1                 ratio,ytot1,f,df,deninv,tempinv


c..various derived constants
      double precision third,sifac,eostol,fpmin
      parameter        (third  = 1.0d0/3.0d0,
     1                  sifac  = 8.6322745944370191d-45,
     2                  eostol = 1.0d-13,
     3                  fpmin  = 1.0d-14)

c..note: sifac = h**3/(2.0d0*pi*amu)**1.5d0





c..popular format statements for debugging
01    format(1x,5(a,1pe24.16))
02    format(1x,5(a,1pe16.8))
03    format(1x,1p5e16.8)



c..set the on/off switches
      radmult  = 1
      ionmult  = 1
      ionized  = 1
      elemult  = 1
      coulmult = 1
      potmult  = 1



c..start pipeline loop
      do j=jlo_eos,jhi_eos

       if (temp_row(j) .le. 0.0) stop 'temp less than 0 in eosfxt'
       if (den_row(j)  .le. 0.0) stop 'den less than 0 in eosfxt'

       temp  = temp_row(j)
       den   = den_row(j)
       abar  = abar_row(j)
       zbar  = zbar_row(j)
       ytot1 = 1.0d0/abar


c..initialize
       prad     = 0.0d0
       dpraddd  = 0.0d0
       dpraddt  = 0.0d0
       dpradda  = 0.0d0
       dpraddz  = 0.0d0

       erad     = 0.0d0
       deraddd  = 0.0d0
       deraddt  = 0.0d0
       deradda  = 0.0d0
       deraddz  = 0.0d0

       srad     = 0.0d0
       dsraddd  = 0.0d0
       dsraddt  = 0.0d0
       dsradda  = 0.0d0
       dsraddz  = 0.0d0

       xni      = 0.0d0
       dxnidd   = 0.0d0
       dxnidt   = 0.0d0
       dxnida   = 0.0d0
       dxnidz   = 0.0d0

       pion     = 0.0d0
       dpiondd  = 0.0d0
       dpiondt  = 0.0d0
       dpionda  = 0.0d0
       dpiondz  = 0.0d0

       eion     = 0.0d0
       deiondd  = 0.0d0
       deiondt  = 0.0d0
       deionda  = 0.0d0
       deiondz  = 0.0d0

       sion     = 0.0d0
       dsiondd  = 0.0d0
       dsiondt  = 0.0d0
       dsionda  = 0.0d0
       dsiondz  = 0.0d0

       xne      = 0.0d0
       dxnedd   = 0.0d0
       dxnedt   = 0.0d0
       dxneda   = 0.0d0
       dxnedz   = 0.0d0

       etaele   = 0.0d0
       detadd   = 0.0d0
       detadt   = 0.0d0
       detada   = 0.0d0
       detadz   = 0.0d0
       etapos   = 0.0d0

       xnefer    = 0.0d0
       dxneferdd = 0.0d0
       dxneferdt = 0.0d0
       dxneferda = 0.0d0
       dxneferdz = 0.0d0

       xnpfer    = 0.0d0
       dxnpferdd = 0.0d0
       dxnpferdt = 0.0d0
       dxnpferda = 0.0d0
       dxnpferdz = 0.0d0

       pele     = 0.0d0
       dpeledd  = 0.0d0
       dpeledt  = 0.0d0
       dpeleda  = 0.0d0
       dpeledz  = 0.0d0

       eele     = 0.0d0
       deeledd  = 0.0d0
       deeledt  = 0.0d0
       deeleda  = 0.0d0
       deeledz  = 0.0d0

       sele     = 0.0d0
       dseledd  = 0.0d0
       dseledt  = 0.0d0
       dseleda  = 0.0d0
       dseledz  = 0.0d0

       ppos     = 0.0d0
       dpposdd  = 0.0d0
       dpposdt  = 0.0d0
       dpposda  = 0.0d0
       dpeledz  = 0.0d0

       epos     = 0.0d0
       deposdd  = 0.0d0
       deposdt  = 0.0d0
       deposda  = 0.0d0
       deeledz  = 0.0d0

       spos     = 0.0d0
       dsposdd  = 0.0d0
       dsposdt  = 0.0d0
       dsposda  = 0.0d0
       dseledz  = 0.0d0

       pep      = 0.0d0
       dpepdd   = 0.0d0
       dpepdt   = 0.0d0
       dpepda   = 0.0d0
       dpepdz   = 0.0d0

       eep      = 0.0d0
       deepdd   = 0.0d0
       deepdt   = 0.0d0
       deepda   = 0.0d0
       deepdz   = 0.0d0

       sep      = 0.0d0
       dsepdd   = 0.0d0
       dsepdt   = 0.0d0
       dsepda   = 0.0d0
       dsepdz   = 0.0d0

       eip      = 0.0d0
       deipdd   = 0.0d0
       deipdt   = 0.0d0
       deipda   = 0.0d0
       deipdz   = 0.0d0

       sip      = 0.0d0
       dsipdd   = 0.0d0
       dsipdt   = 0.0d0
       dsipda   = 0.0d0
       dsipdz   = 0.0d0

       pcoul    = 0.0d0
       dpcouldd = 0.0d0
       dpcouldt = 0.0d0
       dpcoulda = 0.0d0
       dpcouldz = 0.0d0

       ecoul    = 0.0d0
       decouldd = 0.0d0
       decouldt = 0.0d0
       decoulda = 0.0d0
       decouldz = 0.0d0

       scoul    = 0.0d0
       dscouldd = 0.0d0
       dscouldt = 0.0d0
       dscoulda = 0.0d0
       dscouldz = 0.0d0


       kt      = kerg * temp
       ktinv   = 1.0d0/kt
       deninv  = 1.0d0/den
       tempinv = 1.0d0/temp




c..radiation section:
       if (radmult .ne. 0) then

c..pressure in erg/cm**3
        prad    = asol * third * temp * temp * temp * temp
        dpraddd = 0.0d0
        dpraddt = 4.0d0 * prad/temp
        dpradda = 0.0d0
        dpraddz = 0.0d0

c..energy in erg/gr
        erad    = 3.0d0 * prad * deninv
        deraddd = -erad * deninv
        deraddt = 3.0d0 * dpraddt * deninv
        deradda = 0.0d0
        deraddz = 0.0d0

c..entropy in erg/g/kelvin
        srad    = (prad*deninv + erad) * tempinv
        dsraddd = (dpraddd*deninv - prad*deninv**2 + deraddd) * tempinv
        dsraddt = (dpraddt*deninv + deraddt - srad)  * tempinv
        dsradda = 0.0d0
        dsraddz = 0.0d0
       end if




c..ion section:

c..number density in 1/cm**3,
        xni     = avo * ytot1 * den
        dxnidd  = avo * ytot1
        dxnidt  = 0.0d0
        dxnida  = -xni * ytot1
        dxnidz  = 0.0d0

       if (ionmult .ne. 0) then

c..pressure in erg/cm**3
        pion    = xni * kt
        dpiondd = dxnidd * kt
        dpiondt = xni * kerg
        dpionda = -pion * ytot1
        dpiondz = 0.0d0

c.. energy in erg/gr
        eion    = 1.5d0 * pion*deninv
        deiondd = (1.5d0 * dpiondd - eion)*deninv
        deiondt = 1.5d0 * dpiondt*deninv
        deionda = 1.5d0 * dpionda*deninv
        deiondz = 0.0d0


c..ion degeneracy parameter (c&g 9.60)
        y         = 1.0d0/(abar*kt)
        yy        = y * sqrt(y)
        z         = xni * sifac * yy
        etaion    = log(z)
        xx        = 1.0d0/xni
        detaiondd = dxnidd*xx
        detaiondt = dxnidt*xx - 1.5d0*tempinv
        detaionda = dxnida*xx - 1.5d0*ytot1
        detaiondz = dxnidz*xx


c..entropy in erg/gr/kelvin
c..the last term is the usual  etaion * kerg * xni/den
c..sometimes called the sacker-tetrode equation

        sion    = (eion + pion*deninv)*tempinv - etaion * kerg*avo*ytot1
        if (sion.lt.0.0d0) then
c~           write(*,*) "Warning: Negative entropy of ions!"
c~           write(*,*) sion
          sion = 0.0d0
        end if

        dsiondd = (deiondd + dpiondd*deninv - pion*deninv**2)*tempinv
     1            - detaiondd * kerg * avo*ytot1

        dsiondt = (deiondt + dpiondt*deninv)*tempinv
     1            - (eion + pion*deninv)*tempinv**2
     2            - detaiondt * kerg * avo*ytot1

        dsionda = (deionda + dpionda*deninv)*tempinv
     1            - detaionda * kerg * avo*ytot1
     2            + etaion * kerg * avo * ytot1**2

        dsiondz = 0.0d0
       end if









c..electron-positron section:
       if (elemult .ne. 0) then

c..make a good guess at the the electron degeneracy parameter eta
        call etages(xni,zbar,temp,ages)
        agesav = ages


c..newton-raphson to get the electron/positron quantities
        eosfail = .false.
        mode   = 0
        do i=1,100

         call xneroot(mode,den,temp,abar,zbar,ionized,
     1                ages,f,df,
     2                etaele,detadd,detadt,detada,detadz,
     3                etapos,zeff,
     4                xne,dxnedd,dxnedt,dxneda,dxnedz,
     5                xnefer,dxneferdd,dxneferdt,dxneferda,dxneferdz,
     6                xnpfer,dxnpferdd,dxnpferdt,dxnpferda,dxnpferdz,
     7                pele,dpeledd,dpeledt,dpeleda,dpeledz,
     8                ppos,dpposdd,dpposdt,dpposda,dpposdz,
     9                pep,dpepdd,dpepdt,dpepda,dpepdz,
     &                eele,deeledd,deeledt,deeleda,deeledz,
     1                epos,deposdd,deposdt,deposda,deposdz,
     2                eep,deepdd,deepdt,deepda,deepdz,
     3                sele,dseledd,dseledt,dseleda,dseledz,
     4                spos,dsposdd,dsposdt,dsposda,dsposdz,
     5                sep,dsepdd,dsepdt,dsepda,dsepdz,
     6                potmult,
     7                eip,deipdd,deipdt,deipda,deipdz,
     8                sip,dsipdd,dsipdt,dsipda,dsipdz)


         if (df .eq. 0.0) goto 11
         ratio   = f/df
         agesnew = ages - ratio
         z       = abs((agesnew - ages)/ages)
         ages    = agesnew
         niter   = i
         if (z .lt. eostol .or. abs(ratio) .le. fpmin) goto 20
        enddo

11      write(6,*)
        write(6,*) 'newton-raphson failed in routine eosfxt'
        write(6,01) 'temp  =',temp,' den =',den
        write(6,01) 'z     =',z,' ages=',ages, ' agesav=',agesav
        write(6,01) 'eostol=',eostol
        write(6,01) 'f/df  =',f/df,' f   =',f,    ' df    =',df
        write(6,01) 'fpmin =',fpmin
        write(6,*)
        call flush(6)
        eosfail = .true.
        return
20      continue



c..with the converged values, get the energy, pressure, and entropy
         mode = 1
         call xneroot(mode,den,temp,abar,zbar,ionized,
     1                ages,f,df,
     2                etaele,detadd,detadt,detada,detadz,
     3                etapos,zeff,
     4                xne,dxnedd,dxnedt,dxneda,dxnedz,
     5                xnefer,dxneferdd,dxneferdt,dxneferda,dxneferdz,
     6                xnpfer,dxnpferdd,dxnpferdt,dxnpferda,dxnpferdz,
     7                pele,dpeledd,dpeledt,dpeleda,dpeledz,
     8                ppos,dpposdd,dpposdt,dpposda,dpposdz,
     9                pep,dpepdd,dpepdt,dpepda,dpepdz,
     &                eele,deeledd,deeledt,deeleda,deeledz,
     1                epos,deposdd,deposdt,deposda,deposdz,
     2                eep,deepdd,deepdt,deepda,deepdz,
     3                sele,dseledd,dseledt,dseleda,dseledz,
     4                spos,dsposdd,dsposdt,dsposda,dsposdz,
     5                sep,dsepdd,dsepdt,dsepda,dsepdz,
     6                potmult,
     7                eip,deipdd,deipdt,deipda,deipdz,
     8                sip,dsipdd,dsipdt,dsipda,dsipdz)

       end if





c..coulomb corretions section:
       if (coulmult .ne. 0) then

        call coulomb(den,temp,abar,zbar,
     1               pion,dpiondd,dpiondt,dpionda,dpiondz,
     2               xne,dxnedd,dxnedt,dxneda,dxnedz,
     3               plasg,
     4               pcoul,dpcouldd,dpcouldt,dpcoulda,dpcouldz,
     5               ecoul,decouldd,decouldt,decoulda,decouldz,
     6               scoul,dscouldd,dscouldt,dscoulda,dscouldz)


c..bomb proof the cooulomb correctins
        x   = prad + pion + pele + ppos + pcoul
        if (x .le. 0.0) then

c         write(6,*)
c         write(6,*) 'coulomb corrections are causing a negative pressure'
c         write(6,*) 'setting all coulomb corrections to zero'
c         write(6,*)

         pcoul    = 0.0d0
         dpcouldd = 0.0d0
         dpcouldt = 0.0d0
         dpcoulda = 0.0d0
         dpcouldz = 0.0d0
         ecoul    = 0.0d0
         decouldd = 0.0d0
         decouldt = 0.0d0
         decoulda = 0.0d0
         decouldz = 0.0d0
         scoul    = 0.0d0
         dscouldd = 0.0d0
         dscouldt = 0.0d0
         dscoulda = 0.0d0
         dscouldz = 0.0d0
        end if
       end if



c..sum all the components
       pres    = prad + pion + pele + ppos + pcoul
       ener    = erad + eion + eele + epos + ecoul + eip
       entr    = srad + sion + sele + spos + scoul + sip

       dpresdd = dpraddd + dpiondd + dpepdd + dpcouldd
       dpresdt = dpraddt + dpiondt + dpepdt + dpcouldt
       dpresda = dpradda + dpionda + dpepda + dpcoulda
       dpresdz = dpraddz + dpiondz + dpepdz + dpcouldz

       denerdd = deraddd + deiondd + deepdd + decouldd + deipdd
       denerdt = deraddt + deiondt + deepdt + decouldt + deipdt
       denerda = deradda + deionda + deepda + decoulda + deipda
       denerdz = deraddz + deiondz + deepdz + decouldz + deipdz

       dentrdd = dsraddd + dsiondd + dsepdd + dscouldd + dsipdd
       dentrdt = dsraddt + dsiondt + dsepdt + dscouldt + dsipdt
       dentrda = dsradda + dsionda + dsepda + dscoulda + dsipda
       dentrdz = dsraddz + dsiondz + dsepdz + dscouldz + dsipdz





c..the temperature and density exponents (c&g 9.81 9.82)
c..the specific heat at constant volume (c&g 9.92)
c..the third adiabatic exponent (c&g 9.93)
c..the first adiabatic exponent (c&g 9.97)
c..the second adiabatic exponent (c&g 9.105)
c..the specific heat at constant pressure (c&g 9.98)
c..and relativistic formula for the sound speed (c&g 14.29)

       zz    = pres/den
       chit  = temp/pres * dpresdt
       chid  = dpresdd/zz
       cv    = denerdt
       x     = zz * chit/(temp * cv)
       gam3  = x + 1.0d0
       gam1  = chit*x + chid
       nabad = x/gam1
       gam2  = 1.0d0/(1.0d0 - nabad)
       cp    = cv * gam1/chid
       z     = 1.0d0 + (ener + clight*clight)/zz
       sound = clight * sqrt(gam1/z)




c..maxwell relations; each is zero if the consistency is perfect
c..delicate subtraction in very degenerate regions causes roundoff error

       dse = temp*dentrdt/denerdt - 1.0d0

       dpe = (denerdd*den**2 + temp*dpresdt)/pres - 1.0d0

       dsp = -dentrdd*(den**2/dpresdt) - 1.0d0




c..store this row
        ptot_row(j)   = pres
        dpt_row(j)    = dpresdt
        dpd_row(j)    = dpresdd
        dpa_row(j)    = dpresda
        dpz_row(j)    = dpresdz

        etot_row(j)   = ener
        det_row(j)    = denerdt
        ded_row(j)    = denerdd
        dea_row(j)    = denerda
        dez_row(j)    = denerdz

        stot_row(j)   = entr
        dst_row(j)    = dentrdt
        dsd_row(j)    = dentrdd
        dsa_row(j)    = dentrda
        dsz_row(j)    = dentrdz

        prad_row(j)   = prad
        erad_row(j)   = erad
        srad_row(j)   = srad

        pion_row(j)   = pion
        dpiont_row(j) = dpiondt
        eion_row(j)   = eion
        sion_row(j)   = sion

        xni_row(j)    = xni

        pele_row(j)   = pele
        ppos_row(j)   = ppos
        dpept_row(j)  = dpepdt
        dpepd_row(j)  = dpepdd
        dpepa_row(j)  = dpepda
        dpepz_row(j)  = dpepdz

        eele_row(j)   = eele
        epos_row(j)   = epos
        deept_row(j)  = deepdt
        deepd_row(j)  = deepdd
        deepa_row(j)  = deepda
        deepz_row(j)  = deepdz

        sele_row(j)   = sele
        spos_row(j)   = spos
        dsept_row(j)  = dsepdt
        dsepd_row(j)  = dsepdd
        dsepa_row(j)  = dsepda
        dsepz_row(j)  = dsepdz

        xnem_row(j)   = xne
        xne_row(j)    = xnefer
        dxnet_row(j)  = dxneferdt + dxnpferdt
        dxned_row(j)  = dxneferdd + dxnpferdd
        dxnea_row(j)  = dxneferda + dxnpferda
        dxnez_row(j)  = dxneferdz + dxnpferdz
        xnp_row(j)    = xnpfer
        zeff_row(j)   = zeff

        etaele_row(j) = etaele
        detat_row(j)  = detadt
        detad_row(j)  = detadd
        detaa_row(j)  = detada
        detaz_row(j)  = detadz
        etapos_row(j) = etapos

        eip_row(j)    = eip
        sip_row(j)    = sip

        pcou_row(j)   = pcoul
        ecou_row(j)   = ecoul
        scou_row(j)   = scoul
        plasg_row(j)  = plasg

        dse_row(j)    = dse
        dpe_row(j)    = dpe
        dsp_row(j)    = dsp

        cv_row(j)     = cv
        cp_row(j)     = cp
        gam1_row(j)   = gam1
        gam2_row(j)   = gam2
        gam3_row(j)   = gam3
        cs_row(j)     = sound


c..for debugging
c        crap1_row(j)   = etaele
c        dcrap1d_row(j) = detadd
c        dcrap1t_row(j) = detadt
c        dcrap1a_row(j) = detada
c        dcrap1z_row(j) = detadz


c..end of pipeline loop
      enddo
      return
      end subroutine








      subroutine xneroot(mode,den,temp,abar,zbar,ionized,
     1                  aa,f,df,
     2                  etaele,detadd,detadt,detada,detadz,
     3                  etapos,zeff,
     4                  xne,dxnedd,dxnedt,dxneda,dxnedz,
     5                  xnefer,dxneferdd,dxneferdt,dxneferda,dxneferdz,
     6                  xnpfer,dxnpferdd,dxnpferdt,dxnpferda,dxnpferdz,
     7                  pele,dpeledd,dpeledt,dpeleda,dpeledz,
     8                  ppos,dpposdd,dpposdt,dpposda,dpposdz,
     9                  pep,dpepdd,dpepdt,dpepda,dpepdz,
     &                  eele,deeledd,deeledt,deeleda,deeledz,
     1                  epos,deposdd,deposdt,deposda,deposdz,
     2                  eep,deepdd,deepdt,deepda,deepdz,
     3                  sele,dseledd,dseledt,dseleda,dseledz,
     4                  spos,dsposdd,dsposdt,dsposda,dsposdz,
     5                  sep,dsepdd,dsepdt,dsepda,dsepdz,
     6                  potmult,
     7                  eip,deipdd,deipdt,deipda,deipdz,
     8                  sip,dsipdd,dsipdt,dsipda,dsipdz)

      include 'implno.dek'
      include 'const.dek'


c..this routine is called by a root finder to find the degeneracy parameter aa
c..where the number density from a saha equation equals the number
c..density as computed by the fermi-dirac integrals.

c..input is the mode (0 = root find mode, 1 = full calculation),
c..temperature temp, density dem, average weight abar, average charge zbar,
c..and the degeneracy parameter (chemical potential/kerg*temp) aa.

c..everything else is output. see the calling routine for the definitions
c..of the variables



c..declare the pass
      integer          mode,ionized,potmult
      double precision den,temp,abar,zbar,
     1                 aa,f,df,
     2                 etaele,detadd,detadt,detada,detadz,
     3                 etapos,zeff,
     4                 xne,dxnedd,dxnedt,dxneda,dxnedz,
     5                 xnefer,dxneferdd,dxneferdt,dxneferda,dxneferdz,
     6                 xnpfer,dxnpferdd,dxnpferdt,dxnpferda,dxnpferdz,
     7                 pele,dpeledd,dpeledt,dpeleda,dpeledz,
     8                 ppos,dpposdd,dpposdt,dpposda,dpposdz,
     9                 pep,dpepdd,dpepdt,dpepda,dpepdz,
     &                 eele,deeledd,deeledt,deeleda,deeledz,
     1                 epos,deposdd,deposdt,deposda,deposdz,
     2                 eep,deepdd,deepdt,deepda,deepdz,
     3                 sele,dseledd,dseledt,dseleda,dseledz,
     4                 spos,dsposdd,dsposdt,dsposda,dsposdz,
     5                 sep,dsepdd,dsepdt,dsepda,dsepdz,
     6                 eip,deipdd,deipdt,deipda,deipdz,
     7                 sip,dsipdd,dsipdt,dsipda,dsipdz




c..local variables
      double precision kt,beta,beta12,beta32,beta52,
     1                 chi,chifac,saha,sfac,
     2                 xni,dxnidd,dxnidt,dxnida,dxnidz,
     3                 f12,f12eta,f12beta,
     4                 f32,f32eta,f32beta,
     5                 f52,f52eta,f52beta,
     6                 dzeff_deta,dzeffdd,dzeffdt,dzeffda,dzeffdz,
     7                 dxne_deta,dxnefer_deta,dxnefer_dbeta,
     8                 detap_deta,detap_dbeta,
     9                 dxnpfer_detap,dxnpfer_deta,dxnpfer_dbeta,
     &                 dxep_deta,dxep_dbeta,
     1                 dpele_deta,dpele_dbeta,deele_deta,deele_dbeta,
     2                 dppos_detap,dppos_deta,dppos_dbeta,
     3                 depos_detap,depos_deta,depos_dbeta,
     4                 dsfac_deta,ytot1,zz,y,yy,ww,denion



      double precision xconst,pconst,econst,mecc,dbetadt,safe
      parameter        (xconst  = 2.4883740912221807d30,
     1                  pconst  = 1.3581730208282635d24,
     2                  econst  = 2.0372595312423953d24,
     3                  mecc    = me * clight * clight,
     4                  dbetadt = kerg/mecc,
     5                  safe    = 0.005d0)


c..note:
c..xconst = 8.0d0 * pi * sqrt(2.0d0) * (me/h)**3 * c**3
c..pconst = xconst * 2.0d0/3.0d0 * me * clight**2
c..econst = xconst * me * clight**2




c..some common factors
      ytot1   = 1.0d0/abar
      kt      = kerg * temp
      beta    = kt/mecc
      beta12  = sqrt(beta)
      beta32  = beta * beta12
      xni     = avo * ytot1 * den
      etaele  = aa


c..ion number density in 1/cm**3,
      xni     = avo * ytot1 * den
      dxnidd  = avo * ytot1
      dxnidt  = 0.0d0
      dxnida  = -xni * ytot1
      dxnidz  = 0.0d0





c..get the number density of free electrons
c..saha is the ratio of the ground state to the ionized state
c..this model is exact for a pure hydrogen composition
c..denion is a crude pressure ionization model

      if (ionized .eq. 0) then

       denion     = 0.1d0
       chi        = hion * ev2erg * zbar

       chifac     = chi/kt
       yy         = chifac - den/denion

       if (yy .gt. 200.0) then
        saha   = 1.0d90
        f      = 0.0d0
        df     = 1.0d0
        etaele = -100.0d0
        if (mode .eq. 0) return

       else if (yy .lt. -200.0) then
        saha = 0.0d0

       else
        ww   = min(200.0d0,chifac + etaele - den/denion)
        saha = 2.0d0 * exp(ww)
       end if


c..assume fully ionized
      else
       denion     = 1.0e30
       chi        = 0.0d0
       chifac     = 0.0d0
       saha       = 0.0d0
      end if



c..the saha factor, effective charge, and
c..the number density of free electrons

      sfac       = 1.0d0/(1.0d0 + saha)
      dsfac_deta = -sfac*sfac*saha

      zeff       = zbar * sfac
      dzeff_deta = zbar * dsfac_deta

      xne        = xni * zeff
      dxne_deta  = xni * dzeff_deta





c..get the fermi-dirac integral electron contribution
      call dfermi(0.5d0, etaele, beta, f12, f12eta, f12beta)
      call dfermi(1.5d0, etaele, beta, f32, f32eta, f32beta)

      zz            = xconst * beta32
      yy            = f12 + beta * f32
      xnefer        = zz * yy
      dxnefer_deta  = zz * (f12eta + beta * f32eta)
      dxnefer_dbeta = xconst * beta12 * (1.5d0 * yy
     1                       +  beta * (f12beta + f32 + beta * f32beta))



c..if the temperature is not too low, get the positron contributions
c..chemical equilibrium means etaele + etapos = eta_photon = 0.
      etapos        = 0.0d0
      detap_deta    = 0.0d0
      detap_dbeta   = 0.0d0
      xnpfer        = 0.0d0
      dxnpfer_detap = 0.0d0
      dxnpfer_dbeta = 0.0d0

      if (beta .gt. 0.02) then
       etapos      = -aa - 2.0d0/beta
       detap_deta  = -1.0d0
       detap_dbeta = 2.0d0/beta**2
       call dfermi(0.5d0, etapos, beta, f12, f12eta, f12beta)
       call dfermi(1.5d0, etapos, beta, f32, f32eta, f32beta)
       xnpfer        = zz * (f12 + beta * f32)
       dxnpfer_detap = zz * (f12eta + beta * f32eta)
       dxnpfer_dbeta = xconst * beta12 * (1.5d0 * (f12 + beta * f32)
     1                 +  beta * (f12beta + f32 + beta * f32beta))
      end if



c..charge neutrality means ne_ionizat = ne_elect - ne_posit
      f  = xnefer - xnpfer - xne


c..derivative of f with eta for newton-like root finders
      df = dxnefer_deta  - dxnpfer_detap * detap_deta  - dxne_deta



c..if we are in root finder mode, return


      if (mode .eq. 0) return





c..if we are not in root finder mode, polish off the calculation


c..all the derivatives are in terms of eta and beta.
c..we want to convert to temperature, density, abar and zbar derivatives.
c..so, after the root find above on eta we have: xne = xnefer - xnpfer
c..taking the derivative of this and solving for the unknown eta derivatives
c..leads to these expressions:


      dxnpfer_deta  = dxnpfer_detap * detap_deta
      dxnpfer_dbeta = dxnpfer_dbeta + dxnpfer_detap * detap_dbeta
      dxep_deta     = dxnefer_deta  - dxnpfer_deta
      dxep_dbeta    = dxnefer_dbeta - dxnpfer_dbeta


      y      = 1.0d0/(dxep_deta - dxne_deta)

      detadd = (xne/den - dxne_deta/denion)  * y
      detadt = -(dxne_deta*chifac/temp + dxep_dbeta*dbetadt) * y
      detada = -xne/abar * y
      detadz = (xne/zbar + dxne_deta * chifac/zbar)*y



c..derivatives of the effective charge
      dzeffdd = dzeff_deta * (detadd - 1.0d0/denion)
      dzeffdt = dzeff_deta * (detadt - chifac/temp)
      dzeffda = dzeff_deta * detada
      dzeffdz = sfac + dzeff_deta * (detadz + chifac/zbar)


c..derivatives of the electron number density
      dxnedd = dxnidd * zeff + xni * dzeffdd
      dxnedt = dxnidt * zeff + xni * dzeffdt
      dxneda = dxnida * zeff + xni * dzeffda
      dxnedz = dxnidz * zeff + xni * dzeffdz



c..derivatives of the fermi integral electron number densities
      dxneferdd = dxnefer_deta * detadd
      dxneferdt = dxnefer_deta * detadt + dxnefer_dbeta * dbetadt
      dxneferda = dxnefer_deta * detada
      dxneferdz = dxnefer_deta * detadz


c..derivatives of the fermi integral positron number densities
      dxnpferdd = dxnpfer_deta * detadd
      dxnpferdt = dxnpfer_deta * detadt + dxnpfer_dbeta * dbetadt
      dxnpferda = dxnpfer_deta * detada
      dxnpferdz = dxnpfer_deta * detadz





c..now get the pressure and energy
c..for the electrons

      beta52  = beta * beta32
      yy      = pconst * beta52
      zz      = econst * beta52

      call dfermi(1.5d0, etaele, beta, f32, f32eta, f32beta)
      call dfermi(2.5d0, etaele, beta, f52, f52eta, f52beta)

      pele        = yy * (f32 + 0.5d0 * beta * f52)
      dpele_deta  = yy * (f32eta + 0.5d0 * beta * f52eta)
      dpele_dbeta = pconst * beta32 * (2.5d0 * (f32 + 0.5d0*beta* f52)
     1              + beta* (f32beta + 0.5d0*f52 + 0.5d0*beta*f52beta))

      eele        = zz * (f32 + beta * f52)
      deele_deta  = zz * (f32eta + beta * f52eta)
      deele_dbeta = econst * beta32 * (2.5d0 * (f32 + beta * f52)
     1              + beta * (f32beta + f52 + beta * f52beta))


c..for the positrons
      ppos        = 0.0d0
      dppos_detap = 0.0d0
      dppos_dbeta = 0.0d0
      epos        = 0.0d0
      depos_detap = 0.0d0
      depos_dbeta = 0.0d0

      if (beta .gt. 0.02) then
       call dfermi(1.5d0, etapos, beta, f32, f32eta, f32beta)
       call dfermi(2.5d0, etapos, beta, f52, f52eta, f52beta)

       ppos        = yy * (f32 + 0.5d0 * beta * f52)
       dppos_detap = yy * (f32eta + 0.5d0*beta *f52eta)
       dppos_dbeta = pconst * beta32 * (2.5d0 * (f32 + 0.5d0*beta*f52)
     1               + beta*(f32beta + 0.5d0*f52 + 0.5d0*beta*f52beta))

       epos        = zz * (f32 + beta * f52)
       depos_detap = zz * (f32eta + beta * f52eta)
       depos_dbeta = econst * beta32 * (2.5d0 * (f32 + beta * f52)
     1               + beta * (f32beta + f52 + beta * f52beta))
      end if


c..derivatives of the electron pressure
      dpeledd = dpele_deta * detadd
      dpeledt = dpele_deta * detadt + dpele_dbeta * dbetadt
      dpeleda = dpele_deta * detada
      dpeledz = dpele_deta * detadz


c..derivatives of the electron energy
      deeledd = deele_deta * detadd
      deeledt = deele_deta * detadt + deele_dbeta * dbetadt
      deeleda = deele_deta * detada
      deeledz = deele_deta * detadz


c..derivatives of the positron pressure
      dppos_deta  = dppos_detap * detap_deta
      dppos_dbeta = dppos_dbeta + dppos_detap * detap_dbeta
      dpposdd     = dppos_deta * detadd
      dpposdt     = dppos_deta * detadt + dppos_dbeta * dbetadt
      dpposda     = dppos_deta * detada
      dpposdz     = dppos_deta * detadz


c..derivatives of the positron energy
      depos_deta  = depos_detap * detap_deta
      depos_dbeta = depos_dbeta + depos_detap * detap_dbeta
      deposdd     = depos_deta * detadd
      deposdt     = depos_deta * detadt + depos_dbeta * dbetadt
      deposda     = depos_deta * detada
      deposdz     = depos_deta * detadz



c..electron+positron pressure and its derivatives
c..note: at high temperatures and low densities, dpepdd is very small
c..and can go negative, so limit it to be positive definite
      pep    = pele    + ppos
      dpepdd = max(dpeledd + dpposdd, 1.0d-30)
      dpepdt = dpeledt + dpposdt
      dpepda = dpeleda + dpposda
      dpepdz = dpeledz + dpposdz


c..electron+positron thermal energy and its derivatives
      eep    = eele    + epos
      deepdd = deeledd + deposdd
      deepdt = deeledt + deposdt
      deepda = deeleda + deposda
      deepdz = deeledz + deposdz




 114  format(1x,1p5e24.16)


c..electron entropy in erg/gr/kelvin and its derivatives
      y       = kerg/den

      sele    = ((pele + eele)/kt - etaele*xnefer) * y

      dseledd = ((dpeledd + deeledd)/kt
     1            - detadd*xnefer)*y
     2            - etaele*dxneferdd*y
     3            - sele/den

      dseledt = ((dpeledt + deeledt)/kt
     1             - detadt*xnefer
     2             - etaele*dxneferdt
     3             - (pele + eele)/(kt*temp))*y

      dseleda = ((dpeleda + deeleda)/kt - detada*xnefer
     1             - etaele*dxneferda)*y

      dseledz = ((dpeledz + deeledz)/kt - detadz*xnefer
     1             - etaele*dxneferdz)*y



c..positron entropy in erg/gr/kelvin and its derivatives
      spos    = ((ppos + epos)/kt - etapos*xnpfer) * y

      dsposdd = ((dpposdd + deposdd)/kt
     1           - detap_deta*detadd*xnpfer
     2           - etapos*dxnpferdd)*y - spos/den

      dsposdt = ((dpposdt + deposdt)/kt
     1           - (detap_deta*detadt + detap_dbeta*dbetadt)*xnpfer
     2           - etapos*dxnpferdt
     3           - (ppos + epos)/(kt*temp))*y

      dsposda = ((dpposda + deposda)/kt
     1           - detap_deta*detada*xnpfer
     2           - etapos*dxnpferda)*y

      dsposdz = ((dpposdz + deposdz)/kt - detap_deta*detadz*xnpfer
     1             - etapos*dxnpferdz)*y


c..and their sum
      sep     = sele + spos
      dsepdd  = dseledd + dsposdd
      dsepdt  = dseledt + dsposdt
      dsepda  = dseleda + dsposda
      dsepdz  = dseledz + dsposdz



c..adjust for the rest mass energy of the positrons
      y       = 2.0d0 * mecc
      epos    = epos    + y * xnpfer
      deposdd = deposdd + y * dxnpferdd
      deposdt = deposdt + y * dxnpferdt
      deposda = deposda + y * dxnpferda
      deposdz = deposdz + y * dxnpferdz


c..and resum
      deepdd = deeledd + deposdd
      deepdt = deeledt + deposdt
      deepda = deeleda + deposda
      deepdz = deeledz + deposdz


c..convert the electron-positron thermal energy in erg/cm**3 to
c..a specific thermal energy in erg/gr

      eele    = eele/den
      deeledd = deeledd/den - eele/den
      deeledt = deeledt/den
      deeleda = deeleda/den
      deeledz = deeledz/den

      epos    = epos/den
      deposdd = deposdd/den - epos/den
      deposdt = deposdt/den
      deposda = deposda/den
      deposdz = deposdz/den


c..and resum
      deepdd = deeledd + deposdd
      deepdt = deeledt + deposdt
      deepda = deeleda + deposda
      deepdz = deeledz + deposdz




c..and take care of the ionization potential contributions
      if (potmult .eq. 0) then
       eip    = 0.0d0
       deipdd = 0.0d0
       deipdt = 0.0d0
       deipda = 0.0d0
       deipdz = 0.0d0
       sip    = 0.0d0
       dsipdd = 0.0d0
       dsipdt = 0.0d0
       dsipda = 0.0d0
       dsipdz = 0.0d0
      else

       eip    = chi * xne
       deipdd = chi * dxnedd
       deipdt = chi * dxnedt
       deipda = chi * dxneda
       deipdz = chi * dxnedz + hion*ev2erg*xne

c..the ionization entropy in erg/gr/kelvin and its derivatives
       y      = kerg/den

c       sip    = (eip/kt - etaele*xne) * y

c       dsipdd = (deipdd/kt
c     1            - detadd*xne)*y
c     2            - etaele*dxnedd*y
c     3            - sip/den

c       dsipdt = (deipdt/kt
c     1             - detadt*xne
c     2             - etaele*dxnedt
c     3             - eip/(kt*temp))*y


       sip    = eip/kt * y

       dsipdd = deipdd/kt*y - sip/den

       dsipdt = (deipdt/kt - eip/(kt*temp))*y

       dsipda = deipda/kt*y

       dsipdz = deipdz/kt*y


c..convert the ionization energy from erg/cm**3 to  erg/gr

       eip    = eip/den
       deipdd = deipdd/den - eip/den
       deipdt = deipdt/den
       deipda = deipda/den
       deipdz = deipdz/den

      end if

      return
      end subroutine xneroot





      subroutine coulomb(den,temp,abar,zbar,
     1                   pion,dpiondd,dpiondt,dpionda,dpiondz,
     2                   xne,dxnedd,dxnedt,dxneda,dxnedz,
     3                   plasg,
     4                   pcoul,dpcouldd,dpcouldt,dpcoulda,dpcouldz,
     5                   ecoul,decouldd,decouldt,decoulda,decouldz,
     6                   scoul,dscouldd,dscouldt,dscoulda,dscouldz)
      include 'implno.dek'
      include 'const.dek'


c..this routine implments coulomb corrections
c..see yakovlev & shalybkov 1989, uniform background corrections

c..input


c..declare the pass
      double precision den,temp,abar,zbar,
     1                 pion,dpiondd,dpiondt,dpionda,dpiondz,
     2                 xne,dxnedd,dxnedt,dxneda,dxnedz,
     3                 plasg,
     4                 pcoul,dpcouldd,dpcouldt,dpcoulda,dpcouldz,
     5                 ecoul,decouldd,decouldt,decoulda,decouldz,
     6                 scoul,dscouldd,dscouldt,dscoulda,dscouldz


c..local variables
c..for the uniform background coulomb correction
      double precision ytot1,kt,ktinv,
     1                 s,dsdd,dsdt,dsda,dsdz,sinv,
     2                 aele,daeledd,daeledt,daeleda,daeledz,aeleinv,
     3                 eplasg,eplasgdd,eplasgdt,eplasgda,eplasgdz,
     4                 aion,
c     5                 lami,inv_lami,lamidd,lamida,lamidz,
     6                 plasgdd,plasgdt,plasgda,plasgdz,
     7                 x,y,z


      double precision u0,a1,b1,c1,d1,e1,a2,b2,c2
      parameter        (a1 = -0.898004d0,
     1                  b1 =  0.96786d0,
     2                  c1 =  0.220703d0,
     3                  d1 = -0.86097d0,
     4                  e1 =  2.5269d0,
     5                  a2 =  0.29561d0,
     6                  b2 =  1.9885d0,
     7                  c2 =  0.288675d0)


c..various derived constants
      double precision third,forth,fiveth,esqu,forthpi
      parameter        (third   = 1.0d0/3.0d0,
     1                  forth   = 4.0d0/3.0d0,
     2                  fiveth  = 5.0d0/3.0d0,
     3                  esqu    = qe*qe,
     4                  forthpi = forth * pi)



c..common variables
       ytot1   = 1.0d0/abar
       kt      = kerg * temp
       ktinv   = 1.0d0/kt



c..yakovlev & shalybkov eqs 5, 9 and 10
       s        = forthpi * xne
       dsdd     = forthpi * dxnedd
       dsdt     = forthpi * dxnedt
       dsda     = forthpi * dxneda
       dsdz     = forthpi * dxnedz
       sinv     = 1.0d0/s

c..electron-sphere radius aele
       aele     = sinv**third
       z        = -third * aele * sinv
       daeledd  = z * dsdd
       daeledt  = z * dsdt
       daeleda  = z * dsda
       daeledz  = z * dsdz
       aeleinv  = 1.0d0/aele

c..electron coupling parameter eplasg
       eplasg   = esqu * ktinv * aeleinv
       z        = -eplasg * aeleinv
       eplasgdd = z * daeledd
       eplasgdt = z * daeledt - eplasg*ktinv*kerg
       eplasgda = z * daeleda
       eplasgdz = z * daeledz

c..ion-sphere radius aion
       x        = zbar**third
       aion     = x * aele

c..ion coupling parameter plasg
       z        = x*x*x*x*x
       plasg    = z * eplasg
       plasgdd  = z * eplasgdd
       plasgdt  = z * eplasgdt
       plasgda  = z * eplasgda
       plasgdz  = z * eplasgdz + fiveth*x*x * eplasg

c       write(6,*)
c       write(6,112) aion,aele
c       write(6,112) plasg,plasgdd,plasgdt,plasgda,plasgdz
c       write(6,*)



c..yakovlev & shalybkov 1989 equations 82, 85, 86, 87
       if (plasg .ge. 1.0) then
        x        = plasg**(0.25d0)
        u0       = a1*plasg + b1*x + c1/x + d1
        ecoul    = pion/den * u0
        pcoul    = third * ecoul * den
        scoul    = -avo*ytot1*kerg *
     1              (3.0d0*b1*x - 5.0d0*c1/x
     1             + d1 * (log(plasg) - 1.0d0) - e1)

        y        = avo/abar*kt*(a1 + 0.25d0/plasg*(b1*x - c1/x))
        decouldd = y * plasgdd
        decouldt = y * plasgdt + ecoul/temp
        decoulda = y * plasgda - ecoul/abar
        decouldz = y * plasgdz

        y        = third * den
        dpcouldd = third * ecoul + y*decouldd
        dpcouldt = y * decouldt
        dpcoulda = y * decoulda
        dpcouldz = y * decouldz


        y        = -avo*kerg/(abar*plasg)*(0.75d0*b1*x +1.25d0*c1/x +d1)
        dscouldd = y * plasgdd
        dscouldt = y * plasgdt
        dscoulda = y * plasgda - scoul/abar
        dscouldz = y * plasgdz


c..yakovlev & shalybkov 1989 equations 102, 103, 104
       else if (plasg .lt. 1.0) then
        x        = plasg*sqrt(plasg)
        y        = plasg**b2
        z        = c2 * x - third * a2 * y
        pcoul    = -pion * z
        ecoul    = 3.0d0 * pcoul/den
        scoul    = -avo/abar*kerg*(c2*x -a2*(b2-1.0d0)/b2*y)

        s        = 1.5d0*c2*x/plasg - third*a2*b2*y/plasg
        dpcouldd = -dpiondd*z - pion*s*plasgdd
        dpcouldt = -dpiondt*z - pion*s*plasgdt
        dpcoulda = -dpionda*z - pion*s*plasgda
        dpcouldz = -dpiondz*z - pion*s*plasgdz

        s        = 3.0d0/den
        decouldd = s * dpcouldd - ecoul/den
        decouldt = s * dpcouldt
        decoulda = s * dpcoulda
        decouldz = s * dpcouldz

        s        = -avo*kerg/(abar*plasg)*(1.5d0*c2*x -a2*(b2-1.0d0)*y)
        dscouldd = s * plasgdd
        dscouldt = s * plasgdt
        dscoulda = s * plasgda - scoul/abar
        dscouldz = s * plasgdz
       end if




      return
      end subroutine coulomb




      subroutine etages(xni,zbar,temp,eta)
      include 'implno.dek'
      include 'const.dek'

c..this routine makes a damn good guess for the electron degeneracy
c..parameter eta.
c..input is the ion number density xni, average charge zbar,
c..average atomic weigt abar, and temperature temp.
c..output is a guess at the chemical potential eta


c..declare the pass
      double precision  xni,zbar,temp,eta

c..declare
      double precision xne,x,y,z,kt,beta,tmkt,xnefac

      double precision rt2,rt3,rtpi,cpf0,cpf1,cpf2,cpf3,
     1                 twoth,fa0,forpi,mecc
      parameter        (rt2     = 1.4142135623730951d0,
     1                  rt3     = 1.7320508075688772d0,
     2                  rtpi    = 1.7724538509055159d0,
     3                  cpf0    = h/(me*clight),
     4                  cpf1    = 3.0d0/(8.0d0*pi) * cpf0**3,
     5                  cpf2    = 4.0d0/cpf1,
     6                  cpf3    = 2.0d0*rt3*rtpi/(rt2*cpf1),
     7                  twoth   = 2.0d0/3.0d0,
     8                  fa0     = 64.0d0/(9.0d0*pi),
     9                  forpi   = 4.0d0 * pi,
     &                  mecc    = me * clight * clight)

c..notes: rt2=sqrt(2)  rt3=sqrt(3)  rtpi=sqrt(pi)


c..for the purposes of guessing eta, assume full ionization
      xne   = xni * zbar
      kt    = kerg * temp
      beta  = kt/mecc


c..number density of ionized electrons (c&g 24.354k) and number density at
c..turning point (c&g 24.354i). if either of these exceed the number density
c..as given by a saha equation, then pairs are important. set alfa = 1/2.

      if (beta .ge. 1.0) then
       x = cpf2 * beta * beta
      else
       x = cpf3 * beta * (1.0d0 + 0.75d0*beta) * exp(-1.0d0/beta)
      end if
      if (x .ge. xne) then
       eta = -0.5d0


c..get the dimensionless number density (c&g 24.313), if it is large apply the
c..formula (c&g 24.309) to get a possible alfa, if not large do a two term
c..binomial expansion on (c&g 24.309) to estimate alfa.

      else
       z = (xne*cpf1)**twoth
       if (z .ge. 1.0e-6) then
        y = (sqrt(z + 1.0d0) - 1.0d0)/beta
       else
        y = z * (1.0d0 - z * 0.25d0) * 0.5d0/beta
       end if


c..isolate the constant in front of the number density integral. if it is
c..small enough run the divine approximation backwards with c&g 24.43. then
c..join it smoothly with the lower limit.

       x = log10(xne**0.6d0/temp)
       if (x .le. 9.5) then
        z = ((1.0d0 + fa0*beta)*sqrt(1.0d0 + fa0*beta*0.5) - 1.0d0)/fa0
        tmkt    = 2.0d0 * me/h * kt/h
        xnefac  = forpi * tmkt * sqrt(tmkt)
        eta = -log(xnefac*rtpi*(0.5d0+0.75d0*z)/xne)
        if (x .ge. 8.5) eta = eta*(9.5d0-x) + y * (1.0d0 - (9.5d0-x))
       else
        eta = y
       end if
      end if

      return
      end subroutine etages






c..routine dfermi gets the fermi-dirac functions and their derivaties
c..routine fdfunc1 forms the integrand of the fermi-dirac functions
c..routine fdfunc2 same as fdfunc but with the change of variable z**2=x
c..routine dqleg010 does 10 point gauss-legendre integration  9 fig accuracy
c..routine dqleg020 does 20 point gauss-legendre integration 14 fig accuracy
c..routine dqleg040 does 40 point gauss-legendre integration 18 fig accuracy
c..routine dqleg080 does 80 point gauss-legendre integration 32 fig accuracy
c..routine dqlag010 does 10 point gauss-laguerre integration  9 fig accuracy
c..routine dqlag020 does 20 point gauss-laguerre integration 14 fig accuracy
c..routine dqlag040 does 40 point gauss-laguerre integration 18 fig accuracy
c..routine dqlag080 does 80 point gauss-laguerre integration 32 fig accuracy



      subroutine dfermi(dk,eta,theta,fd,fdeta,fdtheta)
      include 'implno.dek'
c..
c..this routine computes the fermi-dirac integrals of
c..index dk, with degeneracy parameter eta and relativity parameter theta.
c..input is dk the double precision index of the fermi-dirac function,
c..eta the degeneracy parameter, and theta the relativity parameter.
c..the output is fd is computed by applying three 10-point
c..gauss-legendre and one 10-point gauss-laguerre rules over
c..four appropriate subintervals. the derivative with respect to eta is
c..output in fdeta, and the derivative with respct to theta is in fdtheta.
c..within each subinterval the fd kernel.
c..
c..this routine delivers at least 9 figures of accuracy
c..
c..reference: j.m. aparicio, apjs 117, 632 1998
c..
c..declare
c..declare
!debug: needed to comment line below for undefined references reasons...
!      external         fdfunc1,fdfunc2
      double precision dk,eta,theta,fd,fdeta,fdtheta,
     1                 d,sg,a1,b1,c1,a2,b2,c2,d2,e2,a3,b3,c3,d3,e3,
     2                 eta1,xi,xi2,x1,x2,x3,s1,s2,s3,s12,par(3),
     3                 res1,dres1,ddres1,res2,dres2,ddres2,
     4                 res3,dres3,ddres3,res4,dres4,ddres4


c   parameters defining the location of the breakpoints for the
c   subintervals of integration:
      data d   / 3.3609d 0 /
      data sg  / 9.1186d-2 /
      data a1  / 6.7774d 0 /
      data b1  / 1.1418d 0 /
      data c1  / 2.9826d 0 /
      data a2  / 3.7601d 0 /
      data b2  / 9.3719d-2 /
      data c2  / 2.1063d-2 /
      data d2  / 3.1084d 1 /
      data e2  / 1.0056d 0 /
      data a3  / 7.5669d 0 /
      data b3  / 1.1695d 0 /
      data c3  / 7.5416d-1 /
      data d3  / 6.6558d 0 /
      data e3  /-1.2819d-1 /


c   integrand parameters:
      par(1)=dk
      par(2)=eta
      par(3)=theta


c   definition of xi:
      eta1=sg*(eta-d)
      if (eta1.le.5.d1) then
        xi=log(1.d0+exp(eta1))/sg
      else
        xi=eta-d
      endif
      xi2=xi*xi

c   definition of the x_i:
      x1=(a1  +b1*xi+c1*   xi2)
     +  /(1.d0+c1*xi)
      x2=(a2  +b2*xi+c2*d2*xi2)
     +  /(1.d0+e2*xi+c2*   xi2)
      x3=(a3  +b3*xi+c3*d3*xi2)
     +  /(1.d0+e3*xi+c3*   xi2)

c   breakpoints:
      s1=x1-x2
      s2=x1
      s3=x1+x3
      s12=sqrt(s1)

c   quadrature integrations:

c 9 significant figure accuracy
c      call dqleg010(fdfunc2, 0.d0,  s12, res1, dres1, ddres1, par,3)
c      call dqleg010(fdfunc1,   s1,   s2, res2, dres2, ddres2, par,3)
c      call dqleg010(fdfunc1,   s2,   s3, res3, dres3, ddres3, par,3)
c      call dqlag010(fdfunc1,   s3, 1.d0, res4, dres4, ddres4, par,3)

c 14 significant figure accuracy
      call dqleg020(fdfunc2, 0.d0,  s12, res1, dres1, ddres1, par,3)
      call dqleg020(fdfunc1,   s1,   s2, res2, dres2, ddres2, par,3)
      call dqleg020(fdfunc1,   s2,   s3, res3, dres3, ddres3, par,3)
      call dqlag020(fdfunc1,   s3, 1.d0, res4, dres4, ddres4, par,3)

c 18 significant figure accuracy
c      call dqleg040(fdfunc2, 0.d0,  s12, res1, dres1, ddres1, par,3)
c      call dqleg040(fdfunc1,   s1,   s2, res2, dres2, ddres2, par,3)
c     call dqleg040(fdfunc1,   s2,   s3, res3, dres3, ddres3, par,3)
c     call dqlag040(fdfunc1,   s3, 1.d0, res4, dres4, ddres4, par,3)

c 32 significant figure accuracy
c      call dqleg080(fdfunc2, 0.d0,  s12, res1, dres1, ddres1, par,3)
c      call dqleg080(fdfunc1,   s1,   s2, res2, dres2, ddres2, par,3)
c      call dqleg080(fdfunc1,   s2,   s3, res3, dres3, ddres3, par,3)
c      call dqlag080(fdfunc1,   s3, 1.d0, res4, dres4, ddres4, par,3)


c..sum the contributions
      fd      = res1 + res2 + res3 + res4
      fdeta   = dres1 + dres2 + dres3 + dres4
      fdtheta = ddres1 + ddres2 + ddres3 + ddres4
      return
      end subroutine dfermi




      subroutine fdfunc1(x,par,n,fd,fdeta,fdtheta)
      include 'implno.dek'
c..
c..forms the fermi-dirac integrand and its derivatives with eta and theta.
c..on input x is the integration variable, par(1) is the double precision
c..index, par(2) is the degeneravy parameter, and par(3) is the relativity
c..parameter. on output fd is the integrand, fdeta is the derivative
c..with respect to eta, and fdtheta is the derivative with respect to theta.
c..
c..declare
      integer          n
      double precision x,par(n),dk,eta,theta,fd,fdeta,fdtheta,
     1                 factor,dxst,denom,denom2,xdk,xdkp1

c..initialize
      dk    = par(1)
      eta   = par(2)
      theta = par(3)
      xdk   = x**dk
      xdkp1 = x * xdk
      dxst  = sqrt(1.0d0 + 0.5d0*x*theta)

c   avoid overflow in the exponentials at large x
      if ((x-eta) .lt. 1.0d2) then
       factor  = exp(x-eta)
       denom   = factor + 1.0d0
       fd      = xdk * dxst / denom
       fdeta   = fd * factor / denom
       denom2  = 4.0d0 * dxst * denom
       fdtheta = xdkp1 / denom2

      else
       factor   = exp(eta-x)
       fd       = xdk * dxst * factor
       fdeta    = fd
       denom2   = 4.0d0 * dxst
       fdtheta  = xdkp1/denom2 * factor
      endif

      return
      end subroutine fdfunc1




      subroutine fdfunc2(x,par,n,fd,fdeta,fdtheta)
      include 'implno.dek'
c..
c..forms the fermi-dirac integrand and its derivatives with eta and theta,
c..when the z**2=x variable change has been made.
c..on input x is the integration variable, par(1) is the double precision
c..index, par(2) is the degeneravy parameter, and par(3) is the relativity
c..parameter. on output fd is the integrand, fdeta is the derivative
c..with respect to eta, and fdtheta is the derivative with respect to theta.
c..
c..declare
      integer          n
      double precision x,par(n),dk,eta,theta,fd,fdeta,fdtheta,
     1                 factor,dxst,denom,denom2,xdk,xdkp1,xsq

      dk    = par(1)
      eta   = par(2)
      theta = par(3)
      xsq   = x * x
      xdk   = x**(2.0d0 * dk + 1.0d0)
      xdkp1 = xsq * xdk
      dxst  = sqrt(1.0d0 + 0.5d0 * xsq * theta)

c   avoid an overflow in the denominator at large x:
      if ((xsq-eta) .lt. 1.d2) then
       factor  = exp(xsq - eta)
       denom   = factor + 1.0d0
       fd      = 2.0d0 * xdk * dxst/denom
       fdeta   = fd * factor/denom
       denom2  = 4.0d0 * dxst * denom
       fdtheta = 2.0d0 * xdkp1/denom2

      else
       factor  = exp(eta - xsq)
       fd      = 2.0d0 * xdk * dxst * factor
       fdeta   = fd
       denom2  = 4.0d0 * dxst
       fdtheta = 2.0d0 * xdkp1/denom2 * factor
      endif

      return
      end subroutine fdfunc2





      subroutine dqleg010(f,a,b,result,dresult,ddresult,par,n)
      include 'implno.dek'
c..
c..10 point gauss-legendre rule for the fermi-dirac function and
c..its derivatives with respect to eta and theta.
c..on input f is the name of the subroutine containing the integrand,
c..a is the lower end point of the interval, b is the higher end point,
c..par is an array of constant parameters to be passed to subroutine f,
c..and n is the length of the par array. on output result is the
c..approximation from applying the 10-point gauss-legendre rule,
c..dresult is the derivative with respect to eta, and ddresult is the
c..derivative with respect to theta.
c..
c..note: since the number of nodes is even, zero is not an abscissa.
c..
c..declare
      external         f
      integer          j,n
      double precision a,b,result,dresult,ddresult,par(n),
     1                 absc1,absc2,center,hlfrun,wg(5),xg(5),
     2                 fval1,dfval1,ddfval1,fval2,dfval2,ddfval2

c the abscissae and weights are given for the interval (-1,1).
c xg     - abscissae of the 20-point gauss-legendre rule
c          for half of the usual run (-1,1), i.e.
c          the positive nodes of the 20-point rule
c wg     - weights of the 20-point gauss rule.
c
c abscissae and weights were evaluated with 100 decimal digit arithmetic.

      data xg (  1) /   1.4887433898 1631210884 8260011297 19984 d -1 /
      data xg (  2) /   4.3339539412 9247190799 2659431657 84162 d -1 /
      data xg (  3) /   6.7940956829 9024406234 3273651148 73575 d -1 /
      data xg (  4) /   8.6506336668 8984510732 0966884234 93048 d -1 /
      data xg (  5) /   9.7390652851 7171720077 9640120844 52053 d -1 /

      data wg (  1) /   2.9552422471 4752870173 8929946513 38329 d -1 /
      data wg (  2) /   2.6926671930 9996355091 2269215694 69352 d -1 /
      data wg (  3) /   2.1908636251 5982043995 5349342281 63192 d -1 /
      data wg (  4) /   1.4945134915 0580593145 7763396576 97332 d -1 /
      data wg (  5) /   6.6671344308 6881375935 6880989333 17928 d -2 /


c           list of major variables
c           -----------------------
c
c           absc   - abscissa
c           fval*  - function value
c           result - result of the 10-point gauss formula

      center   = 0.5d0 * (a+b)
      hlfrun   = 0.5d0 * (b-a)
      result   = 0.0d0
      dresult  = 0.0d0
      ddresult = 0.0d0
      do j=1,5
        absc1 = center + hlfrun*xg(j)
        absc2 = center - hlfrun*xg(j)
        call f(absc1, par, n, fval1, dfval1, ddfval1)
        call f(absc2, par, n, fval2, dfval2, ddfval2)
        result   = result + (fval1 + fval2)*wg(j)
        dresult  = dresult + (dfval1 + dfval2)*wg(j)
        ddresult = ddresult + (ddfval1 + ddfval2)*wg(j)
      enddo
      result   = result * hlfrun
      dresult  = dresult * hlfrun
      ddresult = ddresult * hlfrun
      return
      end subroutine dqleg010




      subroutine dqleg020(f,a,b,result,dresult,ddresult,par,n)
      include 'implno.dek'
c..
c..20 point gauss-legendre rule for the fermi-dirac function and
c..its derivatives with respect to eta and theta.
c..on input f is the name of the subroutine containing the integrand,
c..a is the lower end point of the interval, b is the higher end point,
c..par is an array of constant parameters to be passed to subroutine f,
c..and n is the length of the par array. on output result is the
c..approximation from applying the 20-point gauss-legendre rule,
c..dresult is the derivative with respect to eta, and ddresult is the
c..derivative with respect to theta.
c..
c..note: since the number of nodes is even, zero is not an abscissa.
c..
c..declare
      external         f
      integer          j,n
      double precision a,b,result,dresult,ddresult,par(n),
     1                 absc1,absc2,center,hlfrun,wg(10),xg(10),
     2                 fval1,dfval1,ddfval1,fval2,dfval2,ddfval2


c the abscissae and weights are given for the interval (-1,1).
c xg     - abscissae of the 20-point gauss-legendre rule
c          for half of the usual run (-1,1), i.e.
c          the positive nodes of the 20-point rule
c wg     - weights of the 20-point gauss rule.
c
c abscissae and weights were evaluated with 100 decimal digit arithmetic.


      data xg (  1) /   7.6526521133 4973337546 4040939883 82110 d -2 /
      data xg (  2) /   2.2778585114 1645078080 4961953685 74624 d -1 /
      data xg (  3) /   3.7370608871 5419560672 5481770249 27237 d -1 /
      data xg (  4) /   5.1086700195 0827098004 3640509552 50998 d -1 /
      data xg (  5) /   6.3605368072 6515025452 8366962262 85936 d -1 /
      data xg (  6) /   7.4633190646 0150792614 3050703556 41590 d -1 /
      data xg (  7) /   8.3911697182 2218823394 5290617015 20685 d -1 /
      data xg (  8) /   9.1223442825 1325905867 7524412032 98113 d -1 /
      data xg (  9) /   9.6397192727 7913791267 6661311972 77221 d -1 /
      data xg ( 10) /   9.9312859918 5094924786 1223884713 20278 d -1 /

      data wg (  1) /   1.5275338713 0725850698 0843319550 97593 d -1 /
      data wg (  2) /   1.4917298647 2603746787 8287370019 69436 d -1 /
      data wg (  3) /   1.4209610931 8382051329 2983250671 64933 d -1 /
      data wg (  4) /   1.3168863844 9176626898 4944997481 63134 d -1 /
      data wg (  5) /   1.1819453196 1518417312 3773777113 82287 d -1 /
      data wg (  6) /   1.0193011981 7240435036 7501354803 49876 d -1 /
      data wg (  7) /   8.3276741576 7047487247 5814322204 62061 d -2 /
      data wg (  8) /   6.2672048334 1090635695 0653518704 16063 d -2 /
      data wg (  9) /   4.0601429800 3869413310 3995227493 21098 d -2 /
      data wg ( 10) /   1.7614007139 1521183118 6196235185 28163 d -2 /


c           list of major variables
c           -----------------------
c
c           absc   - abscissa
c           fval*  - function value
c           result - result of the 20-point gauss formula


      center   = 0.5d0 * (a+b)
      hlfrun   = 0.5d0 * (b-a)
      result   = 0.0d0
      dresult  = 0.0d0
      ddresult = 0.0d0
      do j=1,10
        absc1 = center + hlfrun*xg(j)
        absc2 = center - hlfrun*xg(j)
        call f(absc1, par, n, fval1, dfval1, ddfval1)
        call f(absc2, par, n, fval2, dfval2, ddfval2)
        result   = result + (fval1 + fval2)*wg(j)
        dresult  = dresult + (dfval1 + dfval2)*wg(j)
        ddresult = ddresult + (ddfval1 + ddfval2)*wg(j)
      enddo
      result   = result * hlfrun
      dresult  = dresult * hlfrun
      ddresult = ddresult * hlfrun
      return
      end subroutine dqleg020





      subroutine dqleg040(f,a,b,result,dresult,ddresult,par,n)
      include 'implno.dek'
c..
c..40 point gauss-legendre rule for the fermi-dirac function and
c..its derivatives with respect to eta and theta.
c..on input f is the name of the subroutine containing the integrand,
c..a is the lower end point of the interval, b is the higher end point,
c..par is an array of constant parameters to be passed to subroutine f,
c..and n is the length of the par array. on output result is the
c..approximation from applying the 40-point gauss-legendre rule,
c..dresult is the derivative with respect to eta, and ddresult is the
c..derivative with respect to theta.
c..
c..note: since the number of nodes is even, zero is not an abscissa.
c..
c..declare
      external         f
      integer          j,n
      double precision a,b,result,dresult,ddresult,par(n),
     1                 absc1,absc2,center,hlfrun,wg(20),xg(20),
     2                 fval1,dfval1,ddfval1,fval2,dfval2,ddfval2


c the abscissae and weights are given for the interval (-1,1).
c xg     - abscissae of the 40-point gauss-legendre rule
c          for half of the usual run (-1,1), i.e.
c          the positive nodes of the 40-point rule
c wg     - weights of the 40-point gauss rule.
c
c abscissae and weights were evaluated with 100 decimal digit arithmetic.

      data xg (  1) /   3.8772417506 0508219331 9344402462 32946 d -2 /
      data xg (  2) /   1.1608407067 5255208483 4512844080 24113 d -1 /
      data xg (  3) /   1.9269758070 1371099715 5168520651 49894 d -1 /
      data xg (  4) /   2.6815218500 7253681141 1843448085 96183 d -1 /
      data xg (  5) /   3.4199409082 5758473007 4924811791 94310 d -1 /
      data xg (  6) /   4.1377920437 1605001524 8797458037 13682 d -1 /
      data xg (  7) /   4.8307580168 6178712908 5665742448 23004 d -1 /
      data xg (  8) /   5.4946712509 5128202075 9313055295 17970 d -1 /
      data xg (  9) /   6.1255388966 7980237952 6124502306 94877 d -1 /
      data xg ( 10) /   6.7195668461 4179548379 3545149614 94109 d -1 /
      data xg ( 11) /   7.2731825518 9927103280 9964517549 30548 d -1 /
      data xg ( 12) /   7.7830565142 6519387694 9715455064 94848 d -1 /
      data xg ( 13) /   8.2461223083 3311663196 3202306660 98773 d -1 /
      data xg ( 14) /   8.6595950321 2259503820 7818083546 19963 d -1 /
      data xg ( 15) /   9.0209880696 8874296728 2533308684 93103 d -1 /
      data xg ( 16) /   9.3281280827 8676533360 8521668452 05716 d -1 /
      data xg ( 17) /   9.5791681921 3791655804 5409994527 59285 d -1 /
      data xg ( 18) /   9.7725994998 3774262663 3702837129 03806 d -1 /
      data xg ( 19) /   9.9072623869 9457006453 0543522213 72154 d -1 /
      data xg ( 20) /   9.9823770971 0559200349 6227024205 86492 d -1 /

      data wg (  1) /   7.7505947978 4248112637 2396295832 63269 d -2 /
      data wg (  2) /   7.7039818164 2479655883 0753428381 02485 d -2 /
      data wg (  3) /   7.6110361900 6262423715 5807592249 48230 d -2 /
      data wg (  4) /   7.4723169057 9682642001 8933626132 46731 d -2 /
      data wg (  5) /   7.2886582395 8040590605 1068344251 78358 d -2 /
      data wg (  6) /   7.0611647391 2867796954 8363085528 68323 d -2 /
      data wg (  7) /   6.7912045815 2339038256 9010823192 39859 d -2 /
      data wg (  8) /   6.4804013456 6010380745 5452956675 27300 d -2 /
      data wg (  9) /   6.1306242492 9289391665 3799640839 85959 d -2 /
      data wg ( 10) /   5.7439769099 3915513666 1773091042 59856 d -2 /
      data wg ( 11) /   5.3227846983 9368243549 9647977226 05045 d -2 /
      data wg ( 12) /   4.8695807635 0722320614 3416044814 63880 d -2 /
      data wg ( 13) /   4.3870908185 6732719916 7468604171 54958 d -2 /
      data wg ( 14) /   3.8782167974 4720176399 7203129044 61622 d -2 /
      data wg ( 15) /   3.3460195282 5478473926 7818308641 08489 d -2 /
      data wg ( 16) /   2.7937006980 0234010984 8915750772 10773 d -2 /
      data wg ( 17) /   2.2245849194 1669572615 0432418420 85732 d -2 /
      data wg ( 18) /   1.6421058381 9078887128 6348488236 39272 d -2 /
      data wg ( 19) /   1.0498284531 1528136147 4217106727 96523 d -2 /
      data wg ( 20) /   4.5212770985 3319125847 1732878185 33272 d -3 /


c           list of major variables
c           -----------------------
c
c           absc   - abscissa
c           fval*  - function value
c           result - result of the 20-point gauss formula


      center   = 0.5d0 * (a+b)
      hlfrun   = 0.5d0 * (b-a)
      result   = 0.0d0
      dresult  = 0.0d0
      ddresult = 0.0d0
      do j=1,20
        absc1 = center + hlfrun*xg(j)
        absc2 = center - hlfrun*xg(j)
        call f(absc1, par, n, fval1, dfval1, ddfval1)
        call f(absc2, par, n, fval2, dfval2, ddfval2)
        result   = result + (fval1 + fval2)*wg(j)
        dresult  = dresult + (dfval1 + dfval2)*wg(j)
        ddresult = ddresult + (ddfval1 + ddfval2)*wg(j)
      enddo
      result   = result * hlfrun
      dresult  = dresult * hlfrun
      ddresult = ddresult * hlfrun
      return
      end subroutine dqleg040





      subroutine dqleg080(f,a,b,result,dresult,ddresult,par,n)
      include 'implno.dek'
c..
c..80 point gauss-legendre rule for the fermi-dirac function and
c..its derivatives with respect to eta and theta.
c..on input f is the name of the subroutine containing the integrand,
c..a is the lower end point of the interval, b is the higher end point,
c..par is an array of constant parameters to be passed to subroutine f,
c..and n is the length of the par array. on output result is the
c..approximation from applying the 80-point gauss-legendre rule,
c..dresult is the derivative with respect to eta, and ddresult is the
c..derivative with respect to theta.
c..
c..note: since the number of nodes is even, zero is not an abscissa.
c..
c..declare
      external         f
      integer          j,n
      double precision a,b,result,dresult,ddresult,par(n),
     1                 absc1,absc2,center,hlfrun,wg(40),xg(40),
     2                 fval1,dfval1,ddfval1,fval2,dfval2,ddfval2


c the abscissae and weights are given for the interval (-1,1).
c xg     - abscissae of the 80-point gauss-legendre rule
c          for half of the usual run (-1,1), i.e.
c          the positive nodes of the 80-point rule
c wg     - weights of the 80-point gauss rule.
c
c abscissae and weights were evaluated with 100 decimal digit arithmetic.


      data xg (  1) /   1.9511383256 7939976543 5123410745 45479 d -2 /
      data xg (  2) /   5.8504437152 4206686289 9332188341 77944 d -2 /
      data xg (  3) /   9.7408398441 5845990632 7845010493 69020 d -2 /
      data xg (  4) /   1.3616402280 9143886559 2410780007 17067 d -1 /
      data xg (  5) /   1.7471229183 2646812559 3390480112 86195 d -1 /
      data xg (  6) /   2.1299450285 7666132572 3885386663 21823 d -1 /
      data xg (  7) /   2.5095235839 2272120493 1588160350 04797 d -1 /
      data xg (  8) /   2.8852805488 4511853109 1393014347 13898 d -1 /
      data xg (  9) /   3.2566437074 7701914619 1129436273 58695 d -1 /
      data xg ( 10) /   3.6230475349 9487315619 0432863589 63588 d -1 /
      data xg ( 11) /   3.9839340588 1969227024 3796425175 33757 d -1 /
      data xg ( 12) /   4.3387537083 1756093062 3867003631 81958 d -1 /
      data xg ( 13) /   4.6869661517 0544477036 0783649358 08657 d -1 /
      data xg ( 14) /   5.0280411188 8784987593 6727503675 68003 d -1 /
      data xg ( 15) /   5.3614592089 7131932019 8572531254 00904 d -1 /
      data xg ( 16) /   5.6867126812 2709784725 4857866248 27158 d -1 /
      data xg ( 17) /   6.0033062282 9751743154 7462991640 06848 d -1 /
      data xg ( 18) /   6.3107577304 6871966247 9283872893 36863 d -1 /
      data xg ( 19) /   6.6085989898 6119801735 9671228443 17234 d -1 /
      data xg ( 20) /   6.8963764434 2027600771 2076124389 35266 d -1 /
      data xg ( 21) /   7.1736518536 2099880254 0682582938 15278 d -1 /
      data xg ( 22) /   7.4400029758 3597272316 5405279309 13673 d -1 /
      data xg ( 23) /   7.6950242013 5041373865 6160687490 26083 d -1 /
      data xg ( 24) /   7.9383271750 4605449948 6393117384 54358 d -1 /
      data xg ( 25) /   8.1695413868 1463470371 1249940122 95707 d -1 /
      data xg ( 26) /   8.3883147358 0255275616 6230439028 67064 d -1 /
      data xg ( 27) /   8.5943140666 3111096977 1921234916 56492 d -1 /
      data xg ( 28) /   8.7872256767 8213828703 7733436391 24407 d -1 /
      data xg ( 29) /   8.9667557943 8770683194 3240719673 95986 d -1 /
      data xg ( 30) /   9.1326310257 1757654164 7336561509 47478 d -1 /
      data xg ( 31) /   9.2845987717 2445795953 0459590754 53133 d -1 /
      data xg ( 32) /   9.4224276130 9872674752 2660045000 01735 d -1 /
      data xg ( 33) /   9.5459076634 3634905493 4815170210 29508 d -1 /
      data xg ( 34) /   9.6548508904 3799251452 2731556714 54998 d -1 /
      data xg ( 35) /   9.7490914058 5727793385 6452300691 36276 d -1 /
      data xg ( 36) /   9.8284857273 8629070418 2880277091 16473 d -1 /
      data xg ( 37) /   9.8929130249 9755531026 5031671366 31385 d -1 /
      data xg ( 38) /   9.9422754096 5688277892 0635036649 11698 d -1 /
      data xg ( 39) /   9.9764986439 8237688899 4942081831 22985 d -1 /
      data xg ( 40) /   9.9955382265 1630629880 0804990945 67184 d -1 /

      data wg (  1) /   3.9017813656 3066548112 8043925275 40483 d -2 /
      data wg (  2) /   3.8958395962 7695311986 2552477226 08223 d -2 /
      data wg (  3) /   3.8839651059 0519689317 7418266878 71658 d -2 /
      data wg (  4) /   3.8661759774 0764633270 7711026715 66912 d -2 /
      data wg (  5) /   3.8424993006 9594231852 1243632949 01384 d -2 /
      data wg (  6) /   3.8129711314 4776383442 0679156573 62019 d -2 /
      data wg (  7) /   3.7776364362 0013974897 7497642632 10547 d -2 /
      data wg (  8) /   3.7365490238 7304900267 0537705783 86691 d -2 /
      data wg (  9) /   3.6897714638 2760088391 5099657340 52192 d -2 /
      data wg ( 10) /   3.6373749905 8359780439 6499104652 28136 d -2 /
      data wg ( 11) /   3.5794393953 4160546028 6158881615 44542 d -2 /
      data wg ( 12) /   3.5160529044 7475934955 2659238869 68812 d -2 /
      data wg ( 13) /   3.4473120451 7539287943 6422673102 98320 d -2 /
      data wg ( 14) /   3.3733214984 6115228166 7516306423 87284 d -2 /
      data wg ( 15) /   3.2941939397 6454013828 3618090195 95361 d -2 /
      data wg ( 16) /   3.2100498673 4877731480 5649028725 06960 d -2 /
      data wg ( 17) /   3.1210174188 1147016424 4286672060 35518 d -2 /
      data wg ( 18) /   3.0272321759 5579806612 2001009090 11747 d -2 /
      data wg ( 19) /   2.9288369583 2678476927 6758601957 91396 d -2 /
      data wg ( 20) /   2.8259816057 2768623967 5319796501 45302 d -2 /
      data wg ( 21) /   2.7188227500 4863806744 1870668054 42598 d -2 /
      data wg ( 22) /   2.6075235767 5651179029 6874360026 92871 d -2 /
      data wg ( 23) /   2.4922535764 1154911051 1784700321 98023 d -2 /
      data wg ( 24) /   2.3731882865 9301012931 9252461356 84162 d -2 /
      data wg ( 25) /   2.2505090246 3324619262 2158968616 87390 d -2 /
      data wg ( 26) /   2.1244026115 7820063887 1073725061 31285 d -2 /
      data wg ( 27) /   1.9950610878 1419989288 9192871511 35633 d -2 /
      data wg ( 28) /   1.8626814208 2990314287 3541415215 72090 d -2 /
      data wg ( 29) /   1.7274652056 2693063585 8420713129 09998 d -2 /
      data wg ( 30) /   1.5896183583 7256880449 0290922917 85257 d -2 /
      data wg ( 31) /   1.4493508040 5090761169 6207458346 05500 d -2 /
      data wg ( 32) /   1.3068761592 4013392937 8682589705 63403 d -2 /
      data wg ( 33) /   1.1624114120 7978269164 6676999543 26348 d -2 /
      data wg ( 34) /   1.0161766041 1030645208 3185035240 69436 d -2 /
      data wg ( 35) /   8.6839452692 6085842640 9452204034 28135 d -3 /
      data wg ( 36) /   7.1929047681 1731275267 5570867956 50747 d -3 /
      data wg ( 37) /   5.6909224514 0319864926 9107117162 01847 d -3 /
      data wg ( 38) /   4.1803131246 9489523673 9304201681 35132 d -3 /
      data wg ( 39) /   2.6635335895 1268166929 3535831668 45546 d -3 /
      data wg ( 40) /   1.1449500031 8694153454 4171941315 63611 d -3 /


c           list of major variables
c           -----------------------
c
c           absc   - abscissa
c           fval*  - function value
c           result - result of the 20-point gauss formula


      center   = 0.5d0 * (a+b)
      hlfrun   = 0.5d0 * (b-a)
      result   = 0.0d0
      dresult  = 0.0d0
      ddresult = 0.0d0
      do j=1,40
        absc1 = center + hlfrun*xg(j)
        absc2 = center - hlfrun*xg(j)
        call f(absc1, par, n, fval1, dfval1, ddfval1)
        call f(absc2, par, n, fval2, dfval2, ddfval2)
        result   = result + (fval1 + fval2)*wg(j)
        dresult  = dresult + (dfval1 + dfval2)*wg(j)
        ddresult = ddresult + (ddfval1 + ddfval2)*wg(j)
      enddo
      result   = result * hlfrun
      dresult  = dresult * hlfrun
      ddresult = ddresult * hlfrun
      return
      end subroutine dqleg080





      subroutine dqlag010(f,a,b,result,dresult,ddresult,par,n)
      include 'implno.dek'
c..
c..10 point gauss-laguerre rule for the fermi-dirac function.
c..on input f is the external function defining the integrand
c..f(x)=g(x)*w(x), where w(x) is the gaussian weight
c..w(x)=exp(-(x-a)/b) and g(x) a smooth function,
c..a is the lower end point of the interval, b is the higher end point,
c..par is an array of constant parameters to be passed to the function f,
c..and n is the length of the par array. on output result is the
c..approximation from applying the 10-point gauss-laguerre rule.
c..since the number of nodes is even, zero is not an abscissa.
c..
c..declare
      external         f
      integer          j,n
      double precision a,b,result,dresult,ddresult,par(n),
     1                 absc,wg(10),xg(10),fval,dfval,ddfval


c the abscissae and weights are given for the interval (0,+inf).
c xg     - abscissae of the 10-point gauss-laguerre rule
c wg     - weights of the 10-point gauss rule. since f yet
c          includes the weight function, the values in wg
c          are actually exp(xg) times the standard
c          gauss-laguerre weights
c
c abscissae and weights were evaluated with 100 decimal digit arithmetic.

      data xg (  1) /   1.3779347054 0492430830 7725056527 11188 d -1 /
      data xg (  2) /   7.2945454950 3170498160 3731216760 78781 d -1 /
      data xg (  3) /   1.8083429017 4031604823 2920075750 60883 d  0 /
      data xg (  4) /   3.4014336978 5489951448 2532221408 39067 d  0 /
      data xg (  5) /   5.5524961400 6380363241 7558486868 76285 d  0 /
      data xg (  6) /   8.3301527467 6449670023 8767197274 52218 d  0 /
      data xg (  7) /   1.1843785837 9000655649 1853891914 16139 d  1 /
      data xg (  8) /   1.6279257831 3781020995 3265393583 36223 d  1 /
      data xg (  9) /   2.1996585811 9807619512 7709019559 44939 d  1 /
      data xg ( 10) /   2.9920697012 2738915599 0879334079 91951 d  1 /

      data wg (  1) /   3.5400973860 6996308762 2268914420 67608 d -1 /
      data wg (  2) /   8.3190230104 3580738109 8296581278 49577 d -1 /
      data wg (  3) /   1.3302885617 4932817875 2792194393 99369 d  0 /
      data wg (  4) /   1.8630639031 1113098976 3988735482 46693 d  0 /
      data wg (  5) /   2.4502555580 8301016607 2693731657 52256 d  0 /
      data wg (  6) /   3.1227641551 3518249615 0818263314 55472 d  0 /
      data wg (  7) /   3.9341526955 6152109865 5812459248 23077 d  0 /
      data wg (  8) /   4.9924148721 9302310201 1485652433 15445 d  0 /
      data wg (  9) /   6.5722024851 3080297518 7668710376 11234 d  0 /
      data wg ( 10) /   9.7846958403 7463069477 0086638718 59813 d  0 /


c           list of major variables
c           -----------------------
c           absc   - abscissa
c           fval*  - function value
c           result - result of the 10-point gauss formula

      result   = 0.0d0
      dresult  = 0.0d0
      ddresult = 0.0d0
      do j=1,10
       absc = a+b*xg(j)
       call f(absc, par, n, fval, dfval, ddfval)
       result   = result + fval*wg(j)
       dresult  = dresult + dfval*wg(j)
       ddresult = ddresult + ddfval*wg(j)
      enddo
      result   = result*b
      dresult  = dresult*b
      ddresult = ddresult*b
      return
      end subroutine dqlag010





      subroutine dqlag020(f,a,b,result,dresult,ddresult,par,n)
      include 'implno.dek'
c..
c..20 point gauss-laguerre rule for the fermi-dirac function.
c..on input f is the external function defining the integrand
c..f(x)=g(x)*w(x), where w(x) is the gaussian weight
c..w(x)=dexp(-(x-a)/b) and g(x) a smooth function,
c..a is the lower end point of the interval, b is the higher end point,
c..par is an array of constant parameters to be passed to the function f,
c..and n is the length of the par array. on output result is the
c..approximation from applying the 20-point gauss-laguerre rule.
c..since the number of nodes is even, zero is not an abscissa.
c..
c..declare
      external         f
      integer          j,n
      double precision a,b,result,dresult,ddresult,par(n),
     1                 absc,wg(20),xg(20),fval,dfval,ddfval


c the abscissae and weights are given for the interval (0,+inf).
c xg     - abscissae of the 20-point gauss-laguerre rule
c wg     - weights of the 20-point gauss rule. since f yet
c          includes the weight function, the values in wg
c          are actually exp(xg) times the standard
c          gauss-laguerre weights
c
c abscissae and weights were evaluated with 100 decimal digit arithmetic.

      data xg (  1) /   7.0539889691 9887533666 8900458421 50958 d -2 /
      data xg (  2) /   3.7212681800 1611443794 2413887611 46636 d -1 /
      data xg (  3) /   9.1658210248 3273564667 7162770741 83187 d -1 /
      data xg (  4) /   1.7073065310 2834388068 7689667413 05070 d  0 /
      data xg (  5) /   2.7491992553 0943212964 5030460494 81338 d  0 /
      data xg (  6) /   4.0489253138 5088692237 4953369133 33219 d  0 /
      data xg (  7) /   5.6151749708 6161651410 4539885651 89234 d  0 /
      data xg (  8) /   7.4590174536 7106330976 8860218371 81759 d  0 /
      data xg (  9) /   9.5943928695 8109677247 3672734282 79837 d  0 /
      data xg ( 10) /   1.2038802546 9643163096 2340929886 55158 d  1 /
      data xg ( 11) /   1.4814293442 6307399785 1267971004 79756 d  1 /
      data xg ( 12) /   1.7948895520 5193760173 6579099261 25096 d  1 /
      data xg ( 13) /   2.1478788240 2850109757 3517036959 46692 d  1 /
      data xg ( 14) /   2.5451702793 1869055035 1867748464 15418 d  1 /
      data xg ( 15) /   2.9932554631 7006120067 1365613516 58232 d  1 /
      data xg ( 16) /   3.5013434240 4790000062 8493590668 81395 d  1 /
      data xg ( 17) /   4.0833057056 7285710620 2956770780 75526 d  1 /
      data xg ( 18) /   4.7619994047 3465021399 4162715285 11211 d  1 /
      data xg ( 19) /   5.5810795750 0638988907 5077344449 72356 d  1 /
      data xg ( 20) /   6.6524416525 6157538186 4031879146 06659 d  1 /

      data wg (  1) /   1.8108006241 8989255451 6754059131 10644 d -1 /
      data wg (  2) /   4.2255676787 8563974520 3441725664 58197 d -1 /
      data wg (  3) /   6.6690954670 1848150373 4821149925 15927 d -1 /
      data wg (  4) /   9.1535237278 3073672670 6046847718 68067 d -1 /
      data wg (  5) /   1.1695397071 9554597380 1478222395 77476 d  0 /
      data wg (  6) /   1.4313549859 2820598636 8449948915 14331 d  0 /
      data wg (  7) /   1.7029811379 8502272402 5332616332 06720 d  0 /
      data wg (  8) /   1.9870158907 9274721410 9218392751 29020 d  0 /
      data wg (  9) /   2.2866357812 5343078546 2228546814 95651 d  0 /
      data wg ( 10) /   2.6058347275 5383333269 4989509540 33323 d  0 /
      data wg ( 11) /   2.9497837342 1395086600 2354168272 85951 d  0 /
      data wg ( 12) /   3.3253957820 0931955236 9519374217 51118 d  0 /
      data wg ( 13) /   3.7422554705 8981092111 7072932653 77811 d  0 /
      data wg ( 14) /   4.2142367102 5188041986 8080637824 78746 d  0 /
      data wg ( 15) /   4.7625184614 9020929695 2921978390 96371 d  0 /
      data wg ( 16) /   5.4217260442 4557430380 3082979899 81779 d  0 /
      data wg ( 17) /   6.2540123569 3242129289 5184903007 07542 d  0 /
      data wg ( 18) /   7.3873143890 5443455194 0300191964 64791 d  0 /
      data wg ( 19) /   9.1513287309 8747960794 3482425529 50528 d  0 /
      data wg ( 20) /   1.2893388645 9399966710 2628712874 85278 d  1 /


c           list of major variables
c           -----------------------
c           absc   - abscissa
c           fval*  - function value
c           result - result of the 20-point gauss formula

      result   = 0.0d0
      dresult  = 0.0d0
      ddresult = 0.0d0
      do j=1,20
       absc = a+b*xg(j)
       call f(absc, par, n, fval, dfval, ddfval)
       result   = result + fval*wg(j)
       dresult  = dresult + dfval*wg(j)
       ddresult = ddresult + ddfval*wg(j)
      enddo
      result   = result*b
      dresult  = dresult*b
      ddresult = ddresult*b
      return
      end subroutine dqlag020




      subroutine dqlag040(f,a,b,result,dresult,ddresult,par,n)
      include 'implno.dek'
c..
c..20 point gauss-laguerre rule for the fermi-dirac function.
c..on input f is the external function defining the integrand
c..f(x)=g(x)*w(x), where w(x) is the gaussian weight
c..w(x)=dexp(-(x-a)/b) and g(x) a smooth function,
c..a is the lower end point of the interval, b is the higher end point,
c..par is an array of constant parameters to be passed to the function f,
c..and n is the length of the par array. on output result is the
c..approximation from applying the 20-point gauss-laguerre rule.
c..since the number of nodes is even, zero is not an abscissa.
c..
c..declare
      external         f
      integer          j,n
      double precision a,b,result,dresult,ddresult,par(n),
     1                 absc,wg(40),xg(40),fval,dfval,ddfval


c the abscissae and weights are given for the interval (0,+inf).
c xg     - abscissae of the 20-point gauss-laguerre rule
c wg     - weights of the 20-point gauss rule. since f yet
c          includes the weight function, the values in wg
c          are actually exp(xg) times the standard
c          gauss-laguerre weights
c
c abscissae and weights were evaluated with 100 decimal digit arithmetic.

      data xg (  1) /   3.5700394308 8883851220 8447128660 08554 d -2 /
      data xg (  2) /   1.8816228315 8698516003 5893462190 95913 d -1 /
      data xg (  3) /   4.6269428131 4576453564 9375245611 90364 d -1 /
      data xg (  4) /   8.5977296397 2934922257 2722246887 22412 d -1 /
      data xg (  5) /   1.3800108205 2733718649 8000329595 26559 d  0 /
      data xg (  6) /   2.0242091359 2282673344 2066002800 13075 d  0 /
      data xg (  7) /   2.7933693535 0681645765 3514486026 64039 d  0 /
      data xg (  8) /   3.6887026779 0827020959 1526351908 68698 d  0 /
      data xg (  9) /   4.7116411465 5497269361 8722836277 47369 d  0 /
      data xg ( 10) /   5.8638508783 4371811427 3164237995 82987 d  0 /
      data xg ( 11) /   7.1472479081 0228825068 5691951979 42362 d  0 /
      data xg ( 12) /   8.5640170175 8616376271 8522042088 13232 d  0 /
      data xg ( 13) /   1.0116634048 4519394068 4962965639 52448 d  1 /
      data xg ( 14) /   1.1807892294 0045848428 4158670436 06304 d  1 /
      data xg ( 15) /   1.3640933712 5370872283 7167636065 01202 d  1 /
      data xg ( 16) /   1.5619285893 3390738372 0196365218 80145 d  1 /
      data xg ( 17) /   1.7746905950 0956630425 7387749542 43772 d  1 /
      data xg ( 18) /   2.0028232834 5748905296 1261481017 51172 d  1 /
      data xg ( 19) /   2.2468249983 4984183513 7178622899 45366 d  1 /
      data xg ( 20) /   2.5072560772 4262037943 9608620940 09769 d  1 /
      data xg ( 21) /   2.7847480009 1688627207 5170414045 57997 d  1 /
      data xg ( 22) /   3.0800145739 4454627007 5438519619 11114 d  1 /
      data xg ( 23) /   3.3938657084 9137196090 9885858628 19990 d  1 /
      data xg ( 24) /   3.7272245880 4760043283 2076099060 74207 d  1 /
      data xg ( 25) /   4.0811492823 8869204661 5567558160 06426 d  1 /
      data xg ( 26) /   4.4568603175 3344627071 2302063449 83559 d  1 /
      data xg ( 27) /   4.8557763533 0599922809 6204880670 67936 d  1 /
      data xg ( 28) /   5.2795611187 2169329693 5202113739 17638 d  1 /
      data xg ( 29) /   5.7301863323 3936274950 3374699589 21651 d  1 /
      data xg ( 30) /   6.2100179072 7751116121 6819905789 89921 d  1 /
      data xg ( 31) /   6.7219370927 1269987990 8027755188 87054 d  1 /
      data xg ( 32) /   7.2695158847 6124621175 2192772426 19385 d  1 /
      data xg ( 33) /   7.8572802911 5713092805 4389683348 12596 d  1 /
      data xg ( 34) /   8.4911231135 7049845427 0156470966 63186 d  1 /
      data xg ( 35) /   9.1789874671 2363769923 3719348062 73153 d  1 /
      data xg ( 36) /   9.9320808717 4468082501 0905416548 68123 d  1 /
      data xg ( 37) /   1.0767244063 9388272520 7967676113 22664 d  2 /
      data xg ( 38) /   1.1712230951 2690688807 6506441235 50702 d  2 /
      data xg ( 39) /   1.2820184198 8255651192 5411043896 31263 d  2 /
      data xg ( 40) /   1.4228004446 9159997888 3488353595 41764 d  2 /

      data wg (  1) /   9.1625471157 4598973115 1169808013 74830 d -2 /
      data wg (  2) /   2.1342058490 5012080007 1933671215 12341 d -1 /
      data wg (  3) /   3.3571811668 0284673880 5107016162 92191 d -1 /
      data wg (  4) /   4.5854093503 3497560385 4323803764 52497 d -1 /
      data wg (  5) /   5.8206816577 9105168990 9963654015 43283 d -1 /
      data wg (  6) /   7.0649521636 7219392989 8300156730 16682 d -1 /
      data wg (  7) /   8.3202690300 3485238099 1129479783 49523 d -1 /
      data wg (  8) /   9.5887819879 4443111448 1226796760 28906 d -1 /
      data wg (  9) /   1.0872761620 3054971575 3869333172 02661 d  0 /
      data wg ( 10) /   1.2174623279 7778097895 4277850665 60948 d  0 /
      data wg ( 11) /   1.3496954913 5676530792 3938594423 94519 d  0 /
      data wg ( 12) /   1.4842549297 7684671120 5611786129 78719 d  0 /
      data wg ( 13) /   1.6214441628 1182197802 3168843164 54527 d  0 /
      data wg ( 14) /   1.7615953746 7676961118 4242204209 81598 d  0 /
      data wg ( 15) /   1.9050746658 9479967668 2993205972 79371 d  0 /
      data wg ( 16) /   2.0522883472 6171671760 1995822729 47454 d  0 /
      data wg ( 17) /   2.2036905532 4509588909 8283443281 40570 d  0 /
      data wg ( 18) /   2.3597925385 2320332354 0373753789 01497 d  0 /
      data wg ( 19) /   2.5211741403 7643299165 3136902874 22820 d  0 /
      data wg ( 20) /   2.6884980554 0884226415 9505447063 74659 d  0 /
      data wg ( 21) /   2.8625278132 1044881203 4763959831 04311 d  0 /
      data wg ( 22) /   3.0441506653 1151710041 0439679543 33670 d  0 /
      data wg ( 23) /   3.2344070972 6353194177 4902394288 67111 d  0 /
      data wg ( 24) /   3.4345293984 2774809220 3984818916 02464 d  0 /
      data wg ( 25) /   3.6459928249 9408907238 9656466994 90434 d  0 /
      data wg ( 26) /   3.8705845972 1651656808 4753202134 44338 d  0 /
      data wg ( 27) /   4.1104986804 3282265583 5822472639 51577 d  0 /
      data wg ( 28) /   4.3684687232 5406347450 8083382729 45025 d  0 /
      data wg ( 29) /   4.6479589840 7446688299 3033998838 83991 d  0 /
      data wg ( 30) /   4.9534461124 0989326218 6961507855 62721 d  0 /
      data wg ( 31) /   5.2908484059 0073657468 7373657188 58968 d  0 /
      data wg ( 32) /   5.6682046090 3297677000 7305290232 63795 d  0 /
      data wg ( 33) /   6.0967964147 4342030593 3760108591 98806 d  0 /
      data wg ( 34) /   6.5931088610 3999953794 4296642062 94899 d  0 /
      data wg ( 35) /   7.1824959955 3689315064 4298016266 99574 d  0 /
      data wg ( 36) /   7.9066663113 8422877369 3107423105 86595 d  0 /
      data wg ( 37) /   8.8408924928 1034652079 1255950630 26792 d  0 /
      data wg ( 38) /   1.0140899265 6211694839 0946003069 40468 d  1 /
      data wg ( 39) /   1.2210021299 2046038985 2264858758 81108 d  1 /
      data wg ( 40) /   1.6705520642 0242974052 4687743985 73553 d  1 /


c           list of major variables
c           -----------------------
c           absc   - abscissa
c           fval*  - function value
c           result - result of the 20-point gauss formula


      result   = 0.0d0
      dresult  = 0.0d0
      ddresult = 0.0d0
      do j=1,40
       absc = a+b*xg(j)
       call f(absc, par, n, fval, dfval, ddfval)
       result   = result + fval*wg(j)
       dresult  = dresult + dfval*wg(j)
       ddresult = ddresult + ddfval*wg(j)
      enddo
      result   = result*b
      dresult  = dresult*b
      ddresult = ddresult*b
      return
      end subroutine dqlag040





      subroutine dqlag080(f,a,b,result,dresult,ddresult,par,n)
      include 'implno.dek'
c..
c..20 point gauss-laguerre rule for the fermi-dirac function.
c..on input f is the external function defining the integrand
c..f(x)=g(x)*w(x), where w(x) is the gaussian weight
c..w(x)=dexp(-(x-a)/b) and g(x) a smooth function,
c..a is the lower end point of the interval, b is the higher end point,
c..par is an array of constant parameters to be passed to the function f,
c..and n is the length of the par array. on output result is the
c..approximation from applying the 20-point gauss-laguerre rule.
c..since the number of nodes is even, zero is not an abscissa.
c..
c..declare
      external         f
      integer          j,n
      double precision a,b,result,dresult,ddresult,par(n),
     1                 absc,wg(80),xg(80),fval,dfval,ddfval


c the abscissae and weights are given for the interval (0,+inf).
c xg     - abscissae of the 20-point gauss-laguerre rule
c wg     - weights of the 20-point gauss rule. since f yet
c          includes the weight function, the values in wg
c          are actually exp(xg) times the standard
c          gauss-laguerre weights
c
c abscissae and weights were evaluated with 100 decimal digit arithmetic.

      data xg (  1) /   1.7960423300 6983655540 1031924740 16803 d -2 /
      data xg (  2) /   9.4639912994 3539888113 9027246521 72943 d -2 /
      data xg (  3) /   2.3262286812 5867569207 7061572163 49831 d -1 /
      data xg (  4) /   4.3199254780 2387480255 7861724977 70411 d -1 /
      data xg (  5) /   6.9282886135 2021839905 7022136354 46867 d -1 /
      data xg (  6) /   1.0152325561 8947143744 6254368599 35350 d  0 /
      data xg (  7) /   1.3993276878 4287277414 4190514309 78382 d  0 /
      data xg (  8) /   1.8452623038 3584513811 1771177695 99966 d  0 /
      data xg (  9) /   2.3532088716 0926152447 2447080161 40181 d  0 /
      data xg ( 10) /   2.9233646865 5542632483 6912342597 32862 d  0 /
      data xg ( 11) /   3.5559523140 4613405944 9673083246 38370 d  0 /
      data xg ( 12) /   4.2512200823 0987808316 4857664485 77637 d  0 /
      data xg ( 13) /   5.0094426336 2016477243 3677068182 06389 d  0 /
      data xg ( 14) /   5.8309215386 0871901982 1271132956 05083 d  0 /
      data xg ( 15) /   6.7159859778 5131711156 5500876351 99430 d  0 /
      data xg ( 16) /   7.6649934948 9177306073 4189090478 23480 d  0 /
      data xg ( 17) /   8.6783308251 6770109543 4422555426 61083 d  0 /
      data xg ( 18) /   9.7564148057 4293071316 5509446173 66591 d  0 /
      data xg ( 19) /   1.0899693371 2878553774 3610010214 89406 d  1 /
      data xg ( 20) /   1.2108646642 3656999007 0548486983 15593 d  1 /
      data xg ( 21) /   1.3383788112 7786473701 6298406038 33297 d  1 /
      data xg ( 22) /   1.4725665943 5085855393 3580762618 38437 d  1 /
      data xg ( 23) /   1.6134864371 6624665791 6585454289 90907 d  1 /
      data xg ( 24) /   1.7612005243 8144378598 6356869435 86520 d  1 /
      data xg ( 25) /   1.9157749684 2412479221 7299702056 74985 d  1 /
      data xg ( 26) /   2.0772799909 7920960924 4193790104 89579 d  1 /
      data xg ( 27) /   2.2457901204 5404583114 0959169508 77516 d  1 /
      data xg ( 28) /   2.4213844068 9586473771 9224693924 47092 d  1 /
      data xg ( 29) /   2.6041466560 1655866929 3900535654 35682 d  1 /
      data xg ( 30) /   2.7941656841 8594655558 2330692936 92111 d  1 /
      data xg ( 31) /   2.9915355964 9009855011 2704121157 37715 d  1 /
      data xg ( 32) /   3.1963560902 2089207107 7488875426 36533 d  1 /
      data xg ( 33) /   3.4087327864 7261898749 8349473428 60505 d  1 /
      data xg ( 34) /   3.6287775928 7814544588 0319884362 16948 d  1 /
      data xg ( 35) /   3.8566091009 2922104582 5630521729 08535 d  1 /
      data xg ( 36) /   4.0923530218 0312671999 0958505955 44326 d  1 /
      data xg ( 37) /   4.3361426651 7312302957 8267604682 19500 d  1 /
      data xg ( 38) /   4.5881194661 2788863456 2664899748 78378 d  1 /
      data xg ( 39) /   4.8484335660 8331891358 7372733535 63006 d  1 /
      data xg ( 40) /   5.1172444544 6070105959 8894323349 07144 d  1 /
      data xg ( 41) /   5.3947216789 5544471206 2102787875 72430 d  1 /
      data xg ( 42) /   5.6810456334 6362231341 2485032441 02122 d  1 /
      data xg ( 43) /   5.9764084342 1099549427 2959612774 71927 d  1 /
      data xg ( 44) /   6.2810148963 9264772036 2729175902 88682 d  1 /
      data xg ( 45) /   6.5950836257 4560573434 6406271607 92248 d  1 /
      data xg ( 46) /   6.9188482420 2362773741 9802886482 37373 d  1 /
      data xg ( 47) /   7.2525587544 2633453588 3896526165 68450 d  1 /
      data xg ( 48) /   7.5964831127 8641748269 4497949747 96502 d  1 /
      data xg ( 49) /   7.9509089629 0888369620 5728262599 80809 d  1 /
      data xg ( 50) /   8.3161456401 0536896630 4295068758 48705 d  1 /
      data xg ( 51) /   8.6925264419 6156234481 1659260404 48396 d  1 /
      data xg ( 52) /   9.0804112300 9407559518 4117278203 18427 d  1 /
      data xg ( 53) /   9.4801894215 9474332072 0718891387 35302 d  1 /
      data xg ( 54) /   9.8922834446 9405791648 0193727380 36790 d  1 /
      data xg ( 55) /   1.0317152750 8039130233 0470941673 45654 d  2 /
      data xg ( 56) /   1.0755298497 7539906327 6078907989 75954 d  2 /
      data xg ( 57) /   1.1207269048 4128333623 9300461662 11013 d  2 /
      data xg ( 58) /   1.1673666467 3503666318 1578881308 01099 d  2 /
      data xg ( 59) /   1.2155154249 0952625566 8638957521 10813 d  2 /
      data xg ( 60) /   1.2652466579 6515540341 5702654316 53573 d  2 /
      data xg ( 61) /   1.3166419525 2120310870 0890863080 06192 d  2 /
      data xg ( 62) /   1.3697924668 6936973947 5706372894 63788 d  2 /
      data xg ( 63) /   1.4248005891 2161601930 8265692004 55232 d  2 /
      data xg ( 64) /   1.4817820245 5004441818 6523848360 07732 d  2 /
      data xg ( 65) /   1.5408684228 1798697859 4174252655 96259 d  2 /
      data xg ( 66) /   1.6022107287 0095715935 2684168930 10646 d  2 /
      data xg ( 67) /   1.6659835193 4053918744 5211797337 12213 d  2 /
      data xg ( 68) /   1.7323907133 4249503830 9065037750 56999 d  2 /
      data xg ( 69) /   1.8016732304 9032317982 4302089977 01523 d  2 /
      data xg ( 70) /   1.8741194967 6963772390 4901345880 21771 d  2 /
      data xg ( 71) /   1.9500802244 1532991450 3904796005 99643 d  2 /
      data xg ( 72) /   2.0299898419 5074937824 8076778237 14777 d  2 /
      data xg ( 73) /   2.1143987049 4836466691 4849046955 42608 d  2 /
      data xg ( 74) /   2.2040236815 1735739654 0442066777 63168 d  2 /
      data xg ( 75) /   2.2998320607 5680004348 4109696758 44754 d  2 /
      data xg ( 76) /   2.4031908705 5841540417 5974604792 19628 d  2 /
      data xg ( 77) /   2.5161587933 0499611167 4449393109 73194 d  2 /
      data xg ( 78) /   2.6421382388 3199102097 6961086914 35553 d  2 /
      data xg ( 79) /   2.7876673304 6004563652 0141725306 11597 d  2 /
      data xg ( 80) /   2.9696651199 5651345758 8528591557 03581 d  2 /

      data wg (  1) /   4.6093103133 0609664705 2513213955 10083 d -2 /
      data wg (  2) /   1.0731300778 3932752564 1503203043 98860 d -1 /
      data wg (  3) /   1.6866442954 7948111794 2204577827 02406 d -1 /
      data wg (  4) /   2.3008808938 4940054411 2571819781 93282 d -1 /
      data wg (  5) /   2.9160130250 2437964832 1693187729 43752 d -1 /
      data wg (  6) /   3.5322675357 5408236352 7231258056 47046 d -1 /
      data wg (  7) /   4.1498817755 0940466187 1976863112 80092 d -1 /
      data wg (  8) /   4.7690979230 2936241314 7770254185 05661 d -1 /
      data wg (  9) /   5.3901621847 4955374499 5076565223 27912 d -1 /
      data wg ( 10) /   6.0133249744 7190529086 7652488407 39512 d -1 /
      data wg ( 11) /   6.6388413639 6680571849 4422407272 99214 d -1 /
      data wg ( 12) /   7.2669716361 4156688973 5672962491 40514 d -1 /
      data wg ( 13) /   7.8979818942 8428531349 7930783987 88294 d -1 /
      data wg ( 14) /   8.5321447143 8152298354 5981624313 62968 d -1 /
      data wg ( 15) /   9.1697398383 3892698590 3429000315 53302 d -1 /
      data wg ( 16) /   9.8110549100 4005747195 0601559842 18607 d -1 /
      data wg ( 17) /   1.0456386258 0654218147 5684456631 76029 d  0 /
      data wg ( 18) /   1.1106039730 0025890771 1247632597 29371 d  0 /
      data wg ( 19) /   1.1760331584 1226175056 6510765192 08666 d  0 /
      data wg ( 20) /   1.2419589444 9809359279 3517618178 71338 d  0 /
      data wg ( 21) /   1.3084153330 3134064261 1885428459 54645 d  0 /
      data wg ( 22) /   1.3754376757 4892843813 1559170934 90796 d  0 /
      data wg ( 23) /   1.4430627938 7849270398 3124172072 47308 d  0 /
      data wg ( 24) /   1.5113291075 8830693847 6550205599 17703 d  0 /
      data wg ( 25) /   1.5802767765 3099415830 2018787231 21659 d  0 /
      data wg ( 26) /   1.6499478528 0267874116 0120428193 55036 d  0 /
      data wg ( 27) /   1.7203864478 1283277182 0042814522 90770 d  0 /
      data wg ( 28) /   1.7916389147 6093832891 4426205276 88915 d  0 /
      data wg ( 29) /   1.8637540486 4909708435 9257090286 88162 d  0 /
      data wg ( 30) /   1.9367833060 3070923513 9254343278 41646 d  0 /
      data wg ( 31) /   2.0107810470 1134222912 6149881755 55546 d  0 /
      data wg ( 32) /   2.0858048023 8741046429 3039785129 89079 d  0 /
      data wg ( 33) /   2.1619155692 4159897378 3163440488 27763 d  0 /
      data wg ( 34) /   2.2391781388 2364652373 4539974474 45645 d  0 /
      data wg ( 35) /   2.3176614611 4651854068 6060480434 96370 d  0 /
      data wg ( 36) /   2.3974390514 4001430514 1172386388 49980 d  0 /
      data wg ( 37) /   2.4785894444 4973417756 3691644552 22527 d  0 /
      data wg ( 38) /   2.5611967035 7790455335 1155092225 72643 d  0 /
      data wg ( 39) /   2.6453509930 6968892850 4634410003 67534 d  0 /
      data wg ( 40) /   2.7311492228 9915138861 4102871311 69260 d  0 /
      data wg ( 41) /   2.8186957777 5934171703 1418737478 11157 d  0 /
      data wg ( 42) /   2.9081033436 8223018934 5502767774 92687 d  0 /
      data wg ( 43) /   2.9994938483 9685626832 4124518299 68724 d  0 /
      data wg ( 44) /   3.0929995346 9357468116 6951083530 33660 d  0 /
      data wg ( 45) /   3.1887641899 4712376429 3652715016 23466 d  0 /
      data wg ( 46) /   3.2869445597 5337531998 3781070122 16956 d  0 /
      data wg ( 47) /   3.3877119796 0397652334 0549087621 54571 d  0 /
      data wg ( 48) /   3.4912542659 8732012281 7324237827 64895 d  0 /
      data wg ( 49) /   3.5977779176 9613046096 2947301749 02943 d  0 /
      data wg ( 50) /   3.7075106900 1745708341 0271556592 28179 d  0 /
      data wg ( 51) /   3.8207046196 5311695152 0299594304 67622 d  0 /
      data wg ( 52) /   3.9376395977 1430720676 8005406573 30923 d  0 /
      data wg ( 53) /   4.0586276133 8354481597 4201161879 88679 d  0 /
      data wg ( 54) /   4.1840178238 1424031850 6076923345 03121 d  0 /
      data wg ( 55) /   4.3142026492 9613425820 0845732179 87912 d  0 /
      data wg ( 56) /   4.4496251505 3655906604 9828201553 77774 d  0 /
      data wg ( 57) /   4.5907880226 3617511042 9598491489 29810 d  0 /
      data wg ( 58) /   4.7382646459 8929537394 7538735058 38770 d  0 /
      data wg ( 59) /   4.8927127796 6692168696 8869367432 83567 d  0 /
      data wg ( 60) /   5.0548916853 4039512820 5725071351 75938 d  0 /
      data wg ( 61) /   5.2256837559 4272391089 2780101660 22467 d  0 /
      data wg ( 62) /   5.4061221337 9727909853 3235123407 17863 d  0 /
      data wg ( 63) /   5.5974264018 4041404016 5536941589 80053 d  0 /
      data wg ( 64) /   5.8010493213 7643943530 6261624553 94841 d  0 /
      data wg ( 65) /   6.0187389387 8099108768 0151515140 26344 d  0 /
      data wg ( 66) /   6.2526224749 1437403092 9342134800 91928 d  0 /
      data wg ( 67) /   6.5053217351 7668675787 4827196636 96133 d  0 /
      data wg ( 68) /   6.7801152120 0777294201 2873479800 59368 d  0 /
      data wg ( 69) /   7.0811712202 5414518776 1743119167 59402 d  0 /
      data wg ( 70) /   7.4138924461 5305421421 6956062266 87752 d  0 /
      data wg ( 71) /   7.7854415484 1612700386 2327403392 30532 d  0 /
      data wg ( 72) /   8.2055734781 4596472333 9050861009 17119 d  0 /
      data wg ( 73) /   8.6880138399 6161871469 4199580582 55237 d  0 /
      data wg ( 74) /   9.2528697341 5578523923 5565062019 79918 d  0 /
      data wg ( 75) /   9.9311447184 0215736008 3709865340 09772 d  0 /
      data wg ( 76) /   1.0773973641 4646829405 7508435229 90655 d  1 /
      data wg ( 77) /   1.1873891246 5097447081 9508877108 77400 d  1 /
      data wg ( 78) /   1.3422885849 7264236139 7349401540 89734 d  1 /
      data wg ( 79) /   1.5919780161 6897924449 5542522001 85978 d  1 /
      data wg ( 80) /   2.1421454296 4372259537 5210361864 15127 d  1 /


c           list of major variables
c           -----------------------
c           absc   - abscissa
c           fval*  - function value
c           result - result of the 20-point gauss formula

      result   = 0.0d0
      dresult  = 0.0d0
      ddresult = 0.0d0
      do j=1,80
       absc = a+b*xg(j)
       call f(absc, par, n, fval, dfval, ddfval)
       result   = result + fval*wg(j)
       dresult  = dresult + dfval*wg(j)
       ddresult = ddresult + ddfval*wg(j)
      enddo
      result   = result*b
      dresult  = dresult*b
      ddresult = ddresult*b
      return
      end subroutine dqlag080
