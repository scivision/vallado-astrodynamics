      module ast2body
      use comm
      use astmath
      use asttime
      implicit none
*   this file contains fundamental astrodynamic procedures and functions
*   using 2-body dynamics. the routines span a wide range of material, and
*   they come from chapters 2, 3, 5, and 11.
*
*                            companion code for
*               fundamentals of astrodynamics and applications
*                                   2007
*                             by david vallado
*
*       (w) 719-573-2600, email dvallado@agi.com
*
*    current :
*              30 may 07  david vallado
*                           3rd edition baseline
*    changes :
*              21 jul 05  david vallado
*                           2nd printing baseline
*              14 May 01  David Vallado
*                           2nd edition baseline
*              23 Nov 87  David Vallado
*                           Original baseline
*
*     *****************************************************************
*
*     Uses object files:
*         Astmath,
*         AstTime
*     Uses common files:
*         Astmath.cmn
*         Astconst.cmn
*         Astreduc.cmn


*         Sinv= ( Ecc * sin(Nu) ) / sqrt( 1.0D0+ 2.0D0*Ecc*cos(Nu) + Ecc*Ecc )
*            he= 1.0D0+ 2.0D0*Ecc*cos(Nu) + Ecc*Ecc
*            IF ( He .gt. 0.000001D0 ) THEN
*                Cosv= ( 1.0D0 + Ecc * cos(Nu) ) / sqrt(he)
*              ELSE
*                Cosv= 0.0D0   ! Fix for special case (from testgau)
*
*            fpa= DACOS( Cosv )
*             ANGLE( r,v  fpa )

      contains
* ------------------------------------------------------------------------------
*
*                           SUBROUTINE rv2coe
*
*  this subroutine finds the classical orbital elements given the Geocentric
*    Equatorial Position and Velocity vectors.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    R           - IJK Position vector            km
*    V           - IJK Velocity vector            km / s
*
*  Outputs       :
*    P           - SemiLatus rectum               km
*    A           - semimajor axis                 km
*    Ecc         - Eccentricity
*    Incl        - inclination                    0.0D0 to Pi rad
*    Omega       - Longitude of Ascending Node    0.0D0 to 2Pi rad
*    Argp        - Argument of Perigee            0.0D0 to 2Pi rad
*    Nu          - True anomaly                   0.0D0 to 2Pi rad
*    M           - Mean anomaly                   0.0D0 to 2Pi rad
*    ArgLat      - Argument of Latitude      (CI) 0.0D0 to 2Pi rad
*    LamTrue     - True Longitude            (CE) 0.0D0 to 2Pi rad
*    LonPer      - Longitude of Periapsis    (EE) 0.0D0 to 2Pi rad
*
*  Locals        :
*    HBar        - Angular Momentum H Vector      km2 / s
*    EBar        - Eccentricity     E Vector
*    NBar        - Line of Nodes    N Vector
*    c1          - V**2 - u/R
*    RDotV       - R DOT V
*    Hk          - Hk norm vector
*    SME         - Specfic Mechanical Energy      km2 / s2
*    i           - index
*    E           - Eccentric, Parabolic,
*                  Hyperbolic Anomaly             rad
*    Temp        - Temporary variable
*    TypeOrbit   - Type of orbit                  EE, EI, CE, CI
*
*  Coupling      :
*    MAG         - Magnitude of a vector
*    CROSS       - CROSS product of two vectors
*    DOT         - DOT product of two vectors
*    ANGLE       - Find the ANGLE between two vectors
*    NEWTONNU    - Find the mean anomaly
*
*  References    :
*    Vallado       2007, 121, Alg 9, Ex 2-5
*
* ------------------------------------------------------------------------------

      SUBROUTINE rv2coe      ( R, V, P, A, Ecc, Incl, Omega, Argp, Nu,
     &                         M, ArgLat, TrueLon, LonPer )
        IMPLICIT NONE
        REAL*8 R(3), V(3), P, A, Ecc, Incl, Omega, Argp, Nu, M, ArgLat,
     &         TrueLon, LonPer

* -----------------------------  Locals  ------------------------------
        REAL*8 c1, RDotV, hk, SME, Hbar(3), Ebar(3), Nbar(3),
     &         E, Temp,  maghbar, magnbar, magr, magv
        INTEGER i
        CHARACTER*2 TypeOrbit

        ! --------------------  Implementation   ----------------------
        magr = norm2( R )
        magv = norm2( V )
        ! ------------------  Find H N and E vectors   ----------------
        CALL CROSS( R, V, HBar )
        maghbar = norm2(Hbar)
        IF ( maghbar .gt. Small ) THEN
            NBar(1)= -HBar(2)
            NBar(2)=  HBar(1)
            NBar(3)=   0.0D0
            magnbar = norm2( Nbar )
            c1 = magv**2 - mu/magr
            RDotV= dot_product( R, V )
            DO i= 1 , 3
                EBar(i)= (c1*R(i) - RDotV*V(i))/mu
              ENDDO

            Ecc = norm2( EBar )

            ! ------------  Find a e and semi-Latus rectum   ----------
            SME= ( magv*magv*0.5D0 ) - ( mu/magr )
            IF ( abs( SME ) .gt. Small ) THEN
                A= -mu / (2.0D0*SME)
              ELSE
                A= Infinite
              ENDIF
            P = maghbar*maghbar/mu

            ! -----------------  Find inclination   -------------------
            Hk= HBar(3)/maghbar
c            IF ( abs( abs(Hk) - 1.0D0 ) .lt. Small ) THEN
c                ! -------------  Equatorial Orbits   ------------------
c                IF ( abs(HBar(3)) .gt. 0.0D0 ) THEN
c                    Hk= DSIGN(1.0D0, HBar(3))
c                  ENDIF
c              ENDIF
            Incl= DACOS( Hk ) 

            ! --------  Determine type of orbit for Later use  --------
            ! ------ Elliptical, Parabolic, Hyperbolic Inclined -------
            TypeOrbit= 'EI' 
            IF ( Ecc .lt. Small ) THEN
                ! ----------------  Circular Equatorial ---------------
                IF ( (Incl.lt.Small).or.(abs(Incl-Pi).lt.Small) ) THEN
                    TypeOrbit= 'CE'
                  ELSE
                    ! --------------  Circular Inclined ---------------
                    TypeOrbit= 'CI'
                  ENDIF
              ELSE
                ! - Elliptical, Parabolic, Hyperbolic Equatorial --
                IF ( (Incl.lt.Small).or.(abs(Incl-Pi).lt.Small) ) THEN
                    TypeOrbit= 'EE'
                  ENDIF
              ENDIF

            ! ----------  Find Longitude of Ascending Node ------------
            IF ( magnbar .gt. Small ) THEN
                Temp= NBar(1) / magnbar
                IF ( abs(Temp) .gt. 1.0D0 ) THEN
                    Temp= DSIGN(1.0D0, Temp)
                  ENDIF
                Omega= DACOS( Temp ) 
                IF ( NBar(2) .lt. 0.0D0 ) THEN
                    Omega= TwoPi - Omega
                  ENDIF
              ELSE
                Omega= Undefined 
              ENDIF

            ! ---------------- Find Argument of perigee ---------------
            IF ( TypeOrbit .eq. 'EI' ) THEN
                CALL ANGLE( NBar, EBar, Argp )
                IF ( EBar(3) .lt. 0.0D0 ) THEN
                    Argp= TwoPi - Argp 
                  ENDIF
              ELSE
                Argp= Undefined 
              ENDIF

            ! ------------  Find True Anomaly at Epoch    -------------
            IF ( TypeOrbit(1:1) .eq. 'E' ) THEN
                CALL ANGLE( EBar, r, Nu )
                IF ( RDotV .lt. 0.0D0 ) THEN
                    Nu= TwoPi - Nu 
                  ENDIF
              ELSE
                Nu= Undefined 
              ENDIF

            ! ----  Find Argument of Latitude - Circular Inclined -----
            IF ( TypeOrbit .eq. 'CI' ) THEN
                CALL ANGLE( NBar, R, ArgLat )
                IF ( R(3) .lt. 0.0D0 ) THEN
                    ArgLat= TwoPi - ArgLat
                  ENDIF
              ELSE
                ArgLat= Undefined 
              ENDIF

            ! -- Find Longitude of Perigee - Elliptical Equatorial ----
            IF ( ( Ecc.gt.Small ) .and. (TypeOrbit.eq.'EE') ) THEN
                Temp= EBar(1)/Ecc
                IF ( abs(Temp) .gt. 1.0D0 ) THEN
                    Temp= DSIGN(1.0D0, Temp)
                  ENDIF
                LonPer= DACOS( Temp ) 
                IF ( EBar(2) .lt. 0.0D0 ) THEN
                    LonPer= TwoPi - LonPer 
                  ENDIF
                IF ( Incl .gt. HalfPi ) THEN
                    LonPer= TwoPi - LonPer
                  ENDIF
              ELSE
                LonPer= Undefined
              ENDIF

            ! -------- Find True Longitude - Circular Equatorial ------
            IF ( ( magr.gt.Small ) .and. ( TypeOrbit.eq.'CE' ) ) THEN
                Temp= R(1)/magr
                IF ( abs(Temp) .gt. 1.0D0 ) THEN
                    Temp= DSIGN(1.0D0, Temp)
                  ENDIF
                TrueLon= DACOS( Temp )
                IF ( R(2) .lt. 0.0D0 ) THEN
                    TrueLon= TwoPi - TrueLon
                  ENDIF
                IF ( Incl .gt. HalfPi ) THEN
                    TrueLon= TwoPi - TrueLon
                  ENDIF
              ELSE
                TrueLon= Undefined
              ENDIF

            ! ------------ Find Mean Anomaly for all orbits -----------
            CALL NEWTONNU(Ecc, Nu, E, M )

         ELSE
           P    = Undefined
           A    = Undefined
           Ecc  = Undefined
           Incl = Undefined
           Omega= Undefined 
           Argp = Undefined 
           Nu   = Undefined 
           M    = Undefined 
           ArgLat  = Undefined 
           TrueLon= Undefined 
           LonPer = Undefined 
         ENDIF 

      RETURN
      END  ! end rv2coe

* ------------------------------------------------------------------------------
*
*                           SUBROUTINE coe2rv
*
*  this subroutine finds the position and velocity vectors in Geocentric
*    Equatorial (IJK) system given the classical orbit elements.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    P           - SemiLatus rectum               km
*    Ecc         - Eccentricity
*    Incl        - inclination                    0.0D0 to Pi rad
*    Omega       - Longitude of Ascending Node    0.0D0 to 2Pi rad
*    Argp        - Argument of Perigee            0.0D0 to 2Pi rad
*    Nu          - True anomaly                   0.0D0 to 2Pi rad
*    ArgLat      - Argument of Latitude      (CI) 0.0D0 to 2Pi rad
*    LamTrue     - True Longitude            (CE) 0.0D0 to 2Pi rad
*    LonPer      - Longitude of Periapsis    (EE) 0.0D0 to 2Pi rad
*
*  Outputs       :
*    R           - IJK Position vector            km
*    V           - IJK Velocity vector            km / s
*
*  Locals        :
*    Temp        - Temporary REAL*8 value
*    Rpqw        - PQW Position vector            km
*    Vpqw        - PQW Velocity vector            km / s
*    SinNu       - Sine of Nu
*    CosNu       - Cosine of Nu
*    TempVec     - PQW Velocity vector
*
*  Coupling      :
*    ROT3        - Rotation about the 3rd axis
*    ROT1        - Rotation about the 1st axis
*
*  References    :
*    Vallado       2007, 126, Alg 10, Ex 2-5
*
* ------------------------------------------------------------------------------
*
      SUBROUTINE coe2rv      ( P, Ecc, Incl, Omega, Argp, Nu, ArgLat,
     &                         TrueLon, LonPer, R, V )
        IMPLICIT NONE
        REAL*8 R(3), V(3), P, Ecc, Incl, Omega, Argp, Nu, ArgLat,
     &         TrueLon, LonPer
* -----------------------------  Locals  ------------------------------
        REAL*8 Rpqw(3), Vpqw(3), TempVec(3), Temp, SinNu, CosNu

        
        

        ! --------------------  Implementation   ----------------------
*       Determine what type of orbit is involved and set up the
*       set up angles for the special cases.
        ! -------------------------------------------------------------
        IF ( Ecc .lt. Small ) THEN
            ! ----------------  Circular Equatorial  ------------------
            IF ( (Incl.lt.Small).or.( abs(Incl-Pi).lt. Small ) ) THEN
                Argp = 0.0D0
                Omega= 0.0D0 
                Nu   = TrueLon 
              ELSE
                ! --------------  Circular Inclined  ------------------
                Argp= 0.0D0
                Nu  = ArgLat 
              ENDIF
          ELSE
            ! ---------------  Elliptical Equatorial  -----------------
            IF ( ( Incl.lt.Small) .or. (abs(Incl-Pi).lt.Small) ) THEN
                Argp = LonPer
                Omega= 0.0D0 
              ENDIF 
          ENDIF

        ! ----------  Form PQW position and velocity vectors ----------
        CosNu= cos(Nu)
        SinNu= sin(Nu)
        Temp = P / (1.0D0 + Ecc*CosNu)
        Rpqw(1)= Temp*CosNu
        Rpqw(2)= Temp*SinNu
        Rpqw(3)=     0.0D0
        IF ( abs(p) .lt. 0.00000001D0 ) THEN
            p= 0.00000001D0
          ENDIF
        Vpqw(1)=    -SinNu    * sqrt(mu/P)
        Vpqw(2)=  (Ecc + CosNu) * sqrt(mu/P)
        Vpqw(3)=      0.0D0

        ! ----------------  Perform transformation to IJK  ------------
        CALL ROT3( Rpqw   , -Argp , TempVec )
        CALL ROT1( TempVec, -Incl , TempVec )
        CALL ROT3( TempVec, -Omega,  R     )

        CALL ROT3( Vpqw   , -Argp , TempVec )
        CALL ROT1( TempVec, -Incl , TempVec )
        CALL ROT3( TempVec, -Omega, V     )

      RETURN
      END
*
* ----------------------------------------------------------------------------
*
*                           function flt2rv.m
*
*  this function transforms  the flight elements - latgc, lon, fpav, az,
*    position and velocity magnitude into an eci position and velocity vector.
*
*  author        : david vallado                  719-573-2600   17 jun 2002
*
*  revisions
*    vallado     - fix extra terms in rtasc calc                  8 oct 2002
*
*  inputs          description                    range / units
*    rmag        - eci position vector magnitude  km
*    vmag        - eci velocity vector magnitude  km/sec
*    latgc       - geocentric latitude            rad
*    lon         - longitude                      rad
*    fpa         - sat flight path angle          rad
*    az          - sat flight path az             rad
*    ttt         - julian centuries of tt         centuries
*    jdut1       - julian date of ut1             days from 4713 bc
*    lod         - excess length of day           sec
*    xp          - polar motion coefficient       arc sec
*    yp          - polar motion coefficient       arc sec
*    terms       - number of terms for ast calculation 0,2
*
*  outputs       :
*    r           - eci position vector            km
*    v           - eci velocity vector            km/s
*
*  locals        :
*    fpav        - sat flight path anglefrom vert rad
*
*  coupling      :
*    none        -
*
*  references    :
*    vallado       2007, 117
*    chobotov            67
*
* ----------------------------------------------------------------------------

      SUBROUTINE flt2rv      ( rmag, vmag, latgc, lon, fpa, az,
     &                         ttt, jdut1, lod, xp, yp, terms, r, v )
        IMPLICIT NONE
        REAL*8 R(3), V(3), rmag, vmag, latgc, lon, fpa, az, ttt, jdut1,
     &         lod, xp, yp
        INTEGER terms

* -----------------------------  Locals  ------------------------------
        REAL*8 recef(3), vecef(3), rtasc, decl, temp, fpav

        
        

        ! --------------------  Implementation   ----------------------
        ! -------- form position vector
        recef(1) = rmag*cos(latgc)*cos(lon)
        recef(2) = rmag*cos(latgc)*sin(lon)
        recef(3) = rmag*sin(latgc)

        ! -------- convert r to eci
        vecef(1) = 0.0D0
        vecef(2) = 0.0D0
        vecef(3) = 0.0D0
c        CALL ECI_ECEF( r,v, 'FROM', rECEF,vECEF,
c     &                 TTT, JDUT1, LOD, xp, yp, terms )

        ! ------------- calculate rtasc and decl ------------------
        temp= sqrt( r(1)*r(1) + r(2)*r(2) )

* v needs to be defined herexxxxxxxxx
        if ( temp .lt. small ) THEN
            rtasc= atan2( v(2) , v(1) )
          else
            rtasc= atan2( r(2) , r(1) )
          ENDIF
        decl= asin( r(3)/rmag )

        ! -------- form velocity vector
        fpav = halfpi - fpa
        v(1)= vmag*( cos(rtasc)*(-cos(az)*sin(fpav)*sin(decl) +
     &                cos(fpav)*cos(decl)) - sin(az)*sin(fpav)*
     &                sin(rtasc) )
        v(2)= vmag*( sin(rtasc)*(-cos(az)*sin(fpav)*sin(decl) +
     &                cos(fpav)*cos(decl)) + sin(az)*sin(fpav)*
     &                cos(rtasc) )
        v(3)= vmag*( cos(az)*cos(decl)*sin(fpav) + cos(fpav)*
     &                sin(decl) )

      RETURN
      END
*
* ----------------------------------------------------------------------------
*
*                           function rv2flt.m
*
*  this function transforms a position and velocity vector into the flight
*    elements - latgc, lon, fpa, az, position and velocity magnitude.
*
*  author        : david vallado                  719-573-2600   17 jun 2002
*
*  revisions
*    vallado     - add terms for ast calculation                 30 sep 2002
*
*  inputs          description                    range / units
*    r           - eci position vector            km
*    v           - eci velocity vector            km/s
*    ttt         - julian centuries of tt         centuries
*    jdut1       - julian date of ut1             days from 4713 bc
*    lod         - excess length of day           sec
*    xp          - polar motion coefficient       arc sec
*    yp          - polar motion coefficient       arc sec
*    terms       - number of terms for ast calculation 0,2
*
*  outputs       :
*    rmag        - eci position vector magnitude  km
*    vmag        - eci velocity vector magnitude  km/sec
*    latgc       - geocentric latitude            rad
*    lon         - longitude                      rad
*    fpa         - sat flight path angle          rad
*    az          - sat flight path az             rad
*
*  locals        :
*    fpav        - sat flight path anglefrom vert rad
*
*  coupling      :
*    none        -
*
*  references    :
*    vallado       2007, 117
*
* ----------------------------------------------------------------------------

      SUBROUTINE rv2flt      ( R, V, ttt, jdut1, lod, xp, yp, terms,
     &                         rmag, vmag, latgc, lon, fpa, az)
        IMPLICIT NONE
        REAL*8 R(3), V(3), rmag, vmag, latgc, lon, fpa, az, ttt, jdut1,
     &         lod, xp, yp, hmag
        INTEGER terms

* -----------------------------  Locals  ------------------------------
        REAL*8 recef(3), vecef(3), temp, fpav, rdotv, hcrossr(3), h(3)

        
        

        ! --------------------  Implementation   ----------------------
        rmag = norm2(r)
        vmag = norm2(v)

        ! -------- convert r to ecef for lat/lon calculation
c        CALL ECI_ECEF( r,v, 'TOO ', rECEF,vECEF,
c     &                 TTT, JDUT1, LOD, xp, yp, terms )

        ! ----------------- find longitude value  ----------------- uses ecef
        temp = sqrt( recef(1)*recef(1) + recef(2)*recef(2) )
        if ( temp .lt. small ) THEN
            lon= atan2( vecef(2), vecef(1) )
          else
            lon= atan2( recef(2), recef(1) )
          ENDIF

*        latgc = atan2( recef(3) , sqrt(recef(1)**2 + recef(2)**2) )
        latgc = asin( recef(3) / rmag )

        CALL cross(r, v, h)
        hmag = norm2(h)
        rdotv= dot_product(r, v)
        fpav= atan2(hmag, rdotv)
        fpa = halfpi - fpav

        CALL cross(h, r, hcrossr)

        az = atan2( r(1)*hcrossr(2) - r(2)*hcrossr(1),
     &               hcrossr(3)*rmag )

      RETURN
      END
*
* ------------------------------------------------------------------------------
*
*                           function eq2rv
*
*  this function finds the classical orbital elements given the equinoctial
*    elements.
*
*  author        : david vallado                  719-573-2600    9 jun 2002
*
*  revisions
*    vallado     - fix elliptical equatorial orbits case         19 oct 2002
*
*  inputs          description                    range / units
*    af          -
*    ag          -
*    n           - mean motion                    rad
*    meanlon     - mean longitude                 rad
*    chi         -
*    psi         -
*
*  outputs       :
*    r           - eci position vector            km
*    v           - eci velocity vector            km/s
*
*  locals        :
*    temp        - temporary variable
*    p           - semilatus rectum               km
*    a           - semimajor axis                 km
*    ecc         - eccentricity
*    incl        - inclination                    0.0  to pi rad
*    omega       - longitude of ascending node    0.0  to 2pi rad
*    argp        - argument of perigee            0.0  to 2pi rad
*    nu          - true anomaly                   0.0  to 2pi rad
*    m           - mean anomaly                   0.0  to 2pi rad
*    arglat      - argument of latitude      (ci) 0.0  to 2pi rad
*    truelon     - true longitude            (ce) 0.0  to 2pi rad
*    lonper      - longitude of periapsis    (ee) 0.0  to 2pi rad
*
*  coupling      :
*
*  references    :
*    vallado       2007, 116
*    chobotov            30
*
* ------------------------------------------------------------------------------

      SUBROUTINE eq2rv       ( af, ag, meanlon, n, chi, psi, r, v)
        IMPLICIT NONE
        REAL*8 R(3), V(3), af, ag, meanlon, n, chi, psi

* -----------------------------  Locals  ------------------------------
        REAL*8 p, a, ecc, incl, omega, argp, nu, m,
     &         arglat, truelon, lonper, e0
        INTEGER fr

        
        

        ! --------------------  Implementation   ----------------------
        arglat  = 999999.1D0
        lonper  = 999999.1D0
        truelon = 999999.1D0

        a = (mu/n**2)**(1.0D0/3.0D0)

        ecc = sqrt (af**2 + ag**2)

        p = a * (1.0D0 - ecc*ecc)

        incl = 2.0D0 * atan( sqrt(chi**2 + psi**2) )

        ! -------- setup retrograde factor ----------------------------
        fr = 1
        ! -------- set this so it only affects i = 180 deg orbits!! ---
        if (abs(incl-pi) .lt. small) THEN
            fr = -1
          ENDIF

        omega = atan2( chi, psi)

        argp = atan2( fr*ag, af ) - atan2( chi, psi )

        if ( ecc .lt. small ) THEN
            ! ----------------  circular equatorial  ------------------
            if ( (incl.lt.small) .or. ( abs(incl-pi).lt.small) ) THEN
                argp = 0.0D0
                omega= 0.0D0
              else
                ! --------------  circular inclined  ------------------
                argp= 0.0D0
              ENDIF
          else
            ! ---------------  elliptical equatorial  -----------------
            if ( (incl.lt.small) .or. (abs(incl-pi).lt.small) ) THEN
                omega= 0.0D0
              ENDIF
          ENDIF

        m = meanlon - omega - argp
        m = dmod(m+twopi, twopi)

        CALL newtonm ( ecc, m, e0, nu )

        ! ----------  fix for elliptical equatorial orbits ------------
        if ( ecc .lt. small ) THEN
            ! ----------------  circular equatorial  ------------------
            if ((incl.lt.small) .or. ( abs(incl-pi).lt. small )) THEN
                argp    = undefined
                omega   = undefined
                truelon = nu
              else
                ! --------------  circular inclined  ------------------
                argp  = undefined
                arglat= nu
              ENDIF
            nu   = undefined
          else
            ! ---------------  elliptical equatorial  -----------------
            if ( ( incl.lt.small) .or. (abs(incl-pi).lt.small) ) THEN
                lonper = argp
                argp    = undefined
                omega   = undefined
              ENDIF
          ENDIF

        ! -------- now convert back to position and velocity vectors
        CALL coe2rv(p, ecc, incl, omega, argp, nu, arglat, truelon,
     &       lonper, r, v)

      RETURN
      END

* ----------------------------------------------------------------------------
*
*                           function rv2eq.m
*
*  this function transforms a position and velocity vector into the flight
*    elements - latgc, lon, fpa, az, position and velocity magnitude.
*
*  author        : david vallado                  719-573-2600    7 jun 2002
*
*  revisions
*    vallado     - fix special orbit types (ee, hyper)           15 mar 2003
*
*  inputs          description                    range / units
*    r           - eci position vector            km
*    v           - eci velocity vector            km/s
*
*  outputs       :
*    af          -
*    ag          -
*    meanlon     - mean longitude                 rad
*    n           - mean motion                    rad/s
*    chi         -
*    psi         -
*
*  locals        :
*    none        -
*
*  coupling      :
*    none        -
*
*  references    :
*    vallado       2007, 116
*    chobotov            30
*
* ----------------------------------------------------------------------------

      SUBROUTINE rv2eq       ( R, V, af, ag, meanlon, n, chi, psi)
        IMPLICIT NONE
        REAL*8 R(3), V(3), af, ag, meanlon, n, chi, psi

* -----------------------------  Locals  ------------------------------
        REAL*8 p, a, ecc, incl, omega, argp, nu, m,
     &         arglat, truelon, lonper
        INTEGER fr

        
        

        ! --------------------  Implementation   ----------------------
        ! -------- convert to classical elements ----------------------
        CALL rv2coe( r, v,
     &       p, a, ecc, incl, omega, argp, nu, m, arglat, truelon,
     &       lonper )

        ! -------- setup retrograde factor ----------------------------
        fr = 1
        ! -------- set this so it only affects i = 180 deg orbits!! ---
        if (abs(incl-pi) .lt. small) THEN
            fr = -1
          ENDIF

        if ( ecc .lt. small ) THEN
            ! ----------------  circular equatorial  ------------------
            if ((incl.lt.small) .or. ( abs(incl-pi).lt. small )) THEN
                argp = 0.0D0
                omega= 0.0D0
              else
                ! --------------  circular inclined  ------------------
                argp= 0.0D0
              ENDIF
          else
            ! ---------------  elliptical equatorial  -----------------
            if ( ( incl.lt.small) .or. (abs(incl-pi).lt.small) ) THEN
                argp = lonper
                omega= 0.0D0
              ENDIF
          ENDIF

        af = ecc * cos(fr*omega + argp)
        ag = ecc * sin(fr*omega + argp)

        if (fr .gt. 0  ) THEN
            chi = dtan(incl*0.5D0) * sin(omega)
            psi = dtan(incl*0.5D0) * cos(omega)
          else
            chi = cot(incl*0.5D0) * sin(omega)
            psi = cot(incl*0.5D0) * cos(omega)
          ENDIF

c        IF (abs(ecc-1.0D0).le.small) THEN
c            n  = 2.0D0 * sqrt(mu/(p*p*p))
c          ELSE
            IF (a.gt.0.0D0) THEN
                n  = sqrt(mu/(a*a*a))
              ELSE
                n  = sqrt(-mu/(a*a*a))
              ENDIF
c          ENDIF

        meanlon = fr*omega + argp + m
        meanlon = dmod(meanlon, twopi)

      RETURN
      END
*
* ----------------------------------------------------------------------------
*
*                           function adbar2rv.m
*
*  this function transforms the adbarv elements (rtasc, decl, fpa, azimuth,
*    position and velocity magnitude) into eci position and velocity vectors.
*
*  author        : david vallado                  719-573-2600    9 jun 2002
*
*  revisions
*                -
*
*  inputs          description                    range / units
*    rmag        - eci position vector magnitude  km
*    vmag        - eci velocity vector magnitude  km/sec
*    rtasc       - right ascension of sateillite  rad
*    decl        - declination of satellite       rad
*    fpav        - sat flight path angle from vertrad
*    az          - sat flight path azimuth        rad
*
*  outputs       :
*    r           - eci position vector            km
*    v           - eci velocity vector            km/s
*
*  locals        :
*    none        -
*
*  coupling      :
*    none        -
*
*  references    :
*    vallado       2007, 117
*    chobotov            70
*
* ----------------------------------------------------------------------------

      SUBROUTINE adbar2rv    ( rmag, vmag, rtasc, decl, fpav, az, r, v )
        IMPLICIT NONE
        REAL*8 R(3), V(3), rmag, vmag, rtasc, decl, fpav, az

        ! --------------------  Implementation   ----------------------
        ! -------- form position vector
        r(1)= rmag*cos(decl)*cos(rtasc)
        r(2)= rmag*cos(decl)*sin(rtasc)
        r(3)= rmag*sin(decl)

        ! -------- form velocity vector
        v(1)= vmag*( cos(rtasc)*(-cos(az)*sin(fpav)*sin(decl) +
     &                cos(fpav)*cos(decl)) - sin(az)*sin(fpav)*
     &                sin(rtasc) )
        v(2)= vmag*( sin(rtasc)*(-cos(az)*sin(fpav)*sin(decl) +
     &                cos(fpav)*cos(decl)) + sin(az)*sin(fpav)*
     &                cos(rtasc) )
        v(3)= vmag*( cos(az)*cos(decl)*sin(fpav) + cos(fpav)*
     &               sin(decl) )

      RETURN
      END
*
* ----------------------------------------------------------------------------
*
*                           function rv2adbar.m
*
*  this function transforms a position and velocity vector into the adbarv
*    elements - rtasc, decl, fpav, azimuth, position and velocity magnitude.
*
*  author        : david vallado                  719-573-2600    9 jun 2002
*
*  revisions
*                -
*
*  inputs          description                    range / units
*    r           - eci position vector            km
*    v           - eci velocity vector            km/s
*
*  outputs       :
*    rmag        - eci position vector magnitude  km
*    vmag        - eci velocity vector magnitude  km/sec
*    rtasc       - right ascension of sateillite  rad
*    decl        - declination of satellite       rad
*    fpav        - sat flight path angle from vertrad
*    az          - sat flight path azimuth        rad
*
*  locals        :
*    none        -
*
*  coupling      :
*    none        -
*
*  references    :
*    vallado       2007, 117
*    chobotov            70
*
* ----------------------------------------------------------------------------

      SUBROUTINE rv2adbar( R, V, rmag, vmag, rtasc, decl, fpav, az)
        IMPLICIT NONE
        REAL*8 R(3), V(3), rmag, vmag, rtasc, decl, fpav, az


* -----------------------------  Locals  ------------------------------
        REAL*8 temp, temp1, h(3), hcrossr(3), rdotv, hmag

       
        ! --------------------  Implementation   ----------------------
        rmag = norm2(r)
        vmag = norm2(v)

        ! ---------------- calculate rtasc and decl -------------------
        temp= sqrt( r(1)*r(1) + r(2)*r(2) )
        if ( temp .lt. small ) THEN
            temp1= sqrt( v(1)*v(1) + v(2)*v(2) )
            if ( abs(temp1) .gt. small ) THEN
                rtasc= atan2( v(2) , v(1) )
              else
                rtasc= 0.0D0
              ENDIF
          else
            rtasc= atan2( r(2), r(1) )
          ENDIF
        decl= asin( r(3)/rmag )

        CALL cross(r, v, h)
        hmag = norm2(h)
        rdotv= dot_product(r, v)
        fpav = atan2(hmag, rdotv)

        CALL cross(h, r, hcrossr)
        az = atan2( r(1)*hcrossr(2) - r(2)*hcrossr(1),
     &               hcrossr(3)*rmag )

      RETURN
      END
*
* ------------------------------------------------------------------------------
*
*                           function rv2rsw
*
*  this function converts position and velocity vectors into radial, along-
*    track, and cross-track coordinates. note that sometimes the middle vector
*    is called in-track.
*
*  author        : david vallado                  719-573-2600    9 jun 2002
*
*  revisions
*                -
*
*  inputs          description                    range / units
*    r           - position vector                km
*    v           - velocity vector                km/s
*
*  outputs       :
*    rrsw        - position vector                km
*    vrsw        - velocity vector                km/s
*    transmat    - transformation matrix 
*
*  locals        :
*    tempvec     - temporary vector
*    rvec,svec,wvec - direction cosines
*
*  coupling      :
*
*
*  references    :
*    vallado       2007, 163
*
* ------------------------------------------------------------------------------

      SUBROUTINE rv2rsw      ( R, V, rrsw, vrsw, transmat )
        IMPLICIT NONE
        REAL*8 R(3), V(3), rrsw(3), vrsw(3), transmat(3,3)

* -----------------------------  Locals  ------------------------------
        REAL*8 rvec(3), svec(3), wvec(3), tempvec(3)

        ! --------------------  Implementation   ----------------------
        ! in order to work correctly each of the components must be
        ! unit vectors
        ! radial component
        CALL norm( r, rvec)

        ! cross-track component
        CALL cross(r, v, tempvec)
        CALL norm( tempvec,wvec )

        ! along-track component
        CALL cross(wvec, rvec, tempvec)
        CALL norm( tempvec, svec )

        ! assemble transformation matrix from to rsw frame (individual
        !  components arranged in row vectors)
        transmat(1,1) = rvec(1)
        transmat(1,2) = rvec(2)
        transmat(1,3) = rvec(3)
        transmat(2,1) = svec(1)
        transmat(2,2) = svec(2)
        transmat(2,3) = svec(3)
        transmat(3,1) = wvec(1)
        transmat(3,2) = wvec(2)
        transmat(3,3) = wvec(3)

        rrsw = matmul(transmat, r)
        vrsw = matmul(transmat, v)

*   alt approach
*       rrsw(1) = norm2(r)
*       rrsw(2) = 0.0D0
*       rrsw(3) = 0.0D0
*       vrsw(1) = dot_product(r, v)/rrsw(1)
*       vrsw(2) = sqrt(v(1)**2 + v(2)**2 + v(3)**2 - vrsw(1)**2)
*       vrsw(3) = 0.0D0

      RETURN
      END

* ------------------------------------------------------------------------------
*
*                           function rv2ntw
*
*  this function converts position and velocity vectors into in-radial,
*    velocity, and cross-track coordinates. note that sometimes the first
*    vector is called along-radial.
*
*  author        : david vallado                  719-573-2600    5 jul 2002
*
*  revisions
*                -
*
*  inputs          description                    range / units
*    r           - position vector                km
*    v           - velocity vector                km/s
*
*  outputs       :
*    rntw        - position vector                km
*    vntw        - velocity vector                km/s
*    transmat    - transformation matrix
*
*  locals        :
*    tempvec     - temporary vector
*    tvec,nvec,wvec - direction cosines
*
*  coupling      :
*
*
*  references    :
*    vallado       2007, 164
*
* ------------------------------------------------------------------------------

      SUBROUTINE rv2ntw      ( R, V, rntw, vntw, transmat )
        IMPLICIT NONE
        REAL*8 R(3), V(3), rntw(3), vntw(3), transmat(3,3)

* -----------------------------  Locals  ------------------------------
        REAL*8 nvec(3), tvec(3), wvec(3), tempvec(3)

        ! --------------------  Implementation   ----------------------
        ! in order to work correctly each of the components must be
        ! unit vectors
        ! in-velocity component
        CALL norm( v, nvec)

        ! cross-track component
        CALL cross(r, v, tempvec)
        CALL norm( tempvec,wvec )

        ! along-radial component
        CALL cross(nvec, wvec, tempvec)
        CALL norm( tempvec, tvec )

        ! assemble transformation matrix from to ntw frame (individual
        !  components arranged in row vectors)
        transmat(1,1) = tvec(1)
        transmat(1,2) = tvec(2)
        transmat(1,3) = tvec(3)
        transmat(2,1) = nvec(1)
        transmat(2,2) = nvec(2)
        transmat(2,3) = nvec(3)
        transmat(3,1) = wvec(1)
        transmat(3,2) = wvec(2)
        transmat(3,3) = wvec(3)

        rntw = matmul( transmat, r)
        vntw = matmul(transmat, v)

      RETURN
      END

* ------------------------------------------------------------------------------
*
*                           SUBROUTINE FINDC2C3
*
*  this subroutine calculates the C2 and C3 functions for use in the Universal
*    Variable calcuLation of z.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    ZNew        - Z variable                     rad2
*
*  Outputs       :
*    C2New       - C2 FUNCTION value
*    C3New       - C3 FUNCTION value
*
*  Locals        :
*    SqrtZ       - Square root of ZNew
*
*  Coupling      :
*    SINH        - Hyperbolic Sine
*    COSH        - Hyperbolic Cosine
*
*  References    :
*    Vallado       2007, 71, Alg 1
*
* ------------------------------------------------------------------------------

      SUBROUTINE FINDC2C3    ( ZNew, C2New, C3New )
        IMPLICIT NONE
        REAL*8 ZNew, C2New, C3New
* -----------------------------  Locals  ------------------------------
        REAL*8 SqrtZ

        

        ! --------------------  Implementation   ----------------------
        IF ( ZNew .gt. Small ) THEN
            SqrtZ = sqrt( ZNew )
            C2New = (1.0D0-cos( SqrtZ )) / ZNew
            C3New = (SqrtZ-sin( SqrtZ )) / ( SqrtZ**3 )
          ELSE
            IF ( ZNew .lt. -Small ) THEN
                SqrtZ = sqrt( -ZNew )
                C2New = (1.0D0-COSH( SqrtZ )) / ZNew 
                C3New = (SINH( SqrtZ ) - SqrtZ) / ( SqrtZ**3 )
              ELSE
                C2New = 0.5D0
                C3New = 1.0D0/6.0D0 
              ENDIF 
          ENDIF 
      RETURN
      END
*
* ------------------------------------------------------------------------------
*
*                           SUBROUTINE NEWTONE
*
*  this subroutine solves Keplers equation when the Eccentric, paraboic, .or.
*    Hyperbolic anomalies are known. The Mean anomaly and true anomaly are
*    calculated.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Ecc         - Eccentricity                   0.0D0 to
*    E0          - Eccentric Anomaly              -2Pi to 2Pi rad
*
*  Outputs       :
*    M           - Mean Anomaly                   0.0D0 to 2Pi rad
*    Nu          - True Anomaly                   0.0D0 to 2Pi rad
*
*  Locals        :
*    Sinv        - Sine of Nu
*    Cosv        - Cosine of Nu
*
*  Coupling      :
*    SINH        - Hyperbolic Sine
*    COSH        - Hyperbolic Cosine
*
*  References    :
*    Vallado       2007, 85, Alg 6
*
* ------------------------------------------------------------------------------

      SUBROUTINE NEWTONE     ( Ecc, E0, M, Nu )
        IMPLICIT NONE
        REAL*8 Ecc, E0, M, Nu
* -----------------------------  Locals  ------------------------------
        Real*8 Sinv, Cosv

        

        ! --------------------  Implementation   ----------------------
        ! ------------------------- Circular --------------------------
        IF ( abs( Ecc ) .lt. Small ) THEN
            M = E0
            Nu= E0 
          ELSE

            ! ----------------------- Elliptical ----------------------
            IF ( Ecc .lt. 0.999D0 ) THEN
                M= E0 - Ecc*sin(E0)
                Sinv= ( sqrt( 1.0D0-Ecc*Ecc ) * sin(E0) ) /
     &                ( 1.0D0-Ecc*cos(E0) )
                Cosv= ( cos(E0)-Ecc ) / ( 1.0D0 - Ecc*cos(E0) ) 
                Nu  = atan2( Sinv, Cosv )
              ELSE

                ! ---------------------- Hyperbolic  ------------------
                IF ( Ecc .gt. 1.0001D0 ) THEN
                    M= Ecc*SINH(E0) - E0
                    Sinv= ( sqrt( Ecc*Ecc-1.0D0 ) * SINH(E0) ) /
     &                    ( 1.0D0 - Ecc*COSH(E0) )
                    Cosv= ( COSH(E0)-Ecc ) / ( 1.0D0 - Ecc*COSH(E0) ) 
                    Nu  = atan2( Sinv, Cosv )
                  ELSE

                    ! -------------------- Parabolic ------------------
                    M= E0 + (1.0D0/3.0D0)*E0*E0*E0
                    Nu= 2.0D0*atan(E0) 
                  ENDIF 
              ENDIF
          ENDIF
      RETURN
      END
*
* ------------------------------------------------------------------------------
*
*                           SUBROUTINE NEWTONM
*
*  this subroutine performs the Newton Rhapson iteration to find the
*    Eccentric Anomaly given the Mean anomaly.  The True Anomaly is also
*    calculated.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Ecc         - Eccentricity                   0.0D0 to
*    M           - Mean Anomaly                   -2Pi to 2Pi rad
*
*  Outputs       :
*    E0          - Eccentric Anomaly              0.0D0 to 2Pi rad
*    Nu          - True Anomaly                   0.0D0 to 2Pi rad
*
*  Locals        :
*    E1          - Eccentric Anomaly, next value  rad
*    Sinv        - Sine of Nu
*    Cosv        - Cosine of Nu
*    Ktr         - Index
*    R1r         - CUBIC roots - 1 to 3
*    R1i         - imaginary component
*    R2r         -
*    R2i         -
*    R3r         -
*    R3i         -
*    S           - Variables for parabolic solution
*    W           - Variables for parabolic solution
*
*  Coupling      :
*    CUBIC       - Solves a CUBIC polynomial

*  References    :
*    Vallado       2001, 73, Alg 2, Ex 2-1
*
* ------------------------------------------------------------------------------

      SUBROUTINE NEWTONM     ( Ecc, M, E0, Nu )
        IMPLICIT NONE
        REAL*8 Ecc, M, E0, Nu

* -----------------------------  Locals  ------------------------------
        INTEGER Ktr, NumIter
        REAL*8 E1, Sinv, Cosv, R1r, R1i, R2r, R2i, R3r, R3i
        

        ! --------------------  Implementation   ----------------------
        NumIter =    50
        ! -------------------------- Hyperbolic  ----------------------
        IF ( (Ecc-1.0D0) .gt. Small ) THEN
           ! -------------------  Initial Guess -----------------------
            IF ( Ecc .lt. 1.6D0 ) THEN
                IF ( ((M.lt.0.0D0).and.(M.gt.-Pi)).or.(M.gt.Pi) ) THEN
                    E0= M - Ecc
                  ELSE
                    E0= M + Ecc
                  ENDIF
              ELSE
                IF ( (Ecc .lt. 3.6D0) .and. (abs(M) .gt. Pi) ) THEN
                    E0= M - DSIGN(1.0D0, M)*Ecc
                  ELSE
                    E0= M/(Ecc-1.0D0)
                  ENDIF
              ENDIF
            Ktr= 1
            E1 = E0 + ( (M-Ecc*SINH(E0)+E0) / (Ecc*COSH(E0) - 1.0D0) )
            DO WHILE ((abs(E1-E0).gt.Small ) .and. ( Ktr.le.NumIter ))
                E0= E1
                E1= E0 + ( ( M - Ecc*SINH(E0) + E0 ) /
     &                     ( Ecc*COSH(E0) - 1.0D0 ) )
                Ktr = Ktr + 1
              ENDDO
            ! ----------------  Find True Anomaly  --------------------
            Sinv= -( sqrt( Ecc*Ecc-1.0D0 ) * SINH(E1) ) /
     &             ( 1.0D0 - Ecc*COSH(E1) )
            Cosv= ( COSH(E1) - Ecc ) / ( 1.0D0 - Ecc*COSH(E1) )
            Nu  = atan2( Sinv, Cosv )
          ELSE
            ! --------------------- Parabolic -------------------------
            IF ( abs( Ecc-1.0D0 ) .lt. Small ) THEN
                CALL CUBIC( 1.0D0/3.0D0, 0.0D0, 1.0D0, -M, R1r, R1i,
     &                      R2r, R2i, R3r, R3i )
                E0= R1r
*                 S = 0.5D0 * (HalfPi - atan( 1.5D0*M ) )
*                 W = atan( DTAN( S )**(1.0D0/3.0D0) )
*                 E0= 2.0D0*COT(2.0D0*W)
                Ktr= 1
                Nu = 2.0D0 * atan(E0)
              ELSE
                ! -------------------- Elliptical ----------------------
                IF ( Ecc .gt. Small ) THEN
                    ! -----------  Initial Guess -------------
                    IF ( ((M .lt. 0.0D0) .and. (M .gt. -Pi)) .or.
     &                   (M .gt. Pi) ) THEN
                        E0= M - Ecc
                      ELSE
                        E0= M + Ecc
                      ENDIF
                    Ktr= 1
                    E1 = E0 + ( M - E0 + Ecc*sin(E0) ) /
     &                        ( 1.0D0 - Ecc*cos(E0) )
                    DO WHILE (( abs(E1-E0) .gt. Small ) .and.
     &                       ( Ktr .le. NumIter ))
                        Ktr = Ktr + 1
                        E0= E1
                        E1= E0 + ( M - E0 + Ecc*sin(E0) ) /
     &                           ( 1.0D0 - Ecc*cos(E0) )
                      ENDDO
                    ! -------------  Find True Anomaly  ---------------
                    Sinv= ( sqrt( 1.0D0-Ecc*Ecc ) * sin(E1) ) /
     &                    ( 1.0D0-Ecc*cos(E1) )
                    Cosv= ( cos(E1)-Ecc ) / ( 1.0D0 - Ecc*cos(E1) )
                    Nu  = atan2( Sinv, Cosv )
                  ELSE
                    ! -------------------- Circular -------------------
                    Ktr= 0
                    Nu= M
                    E0= M
                  ENDIF
              ENDIF
          ENDIF
      RETURN
      END
*
* ------------------------------------------------------------------------------
*
*                           SUBROUTINE NEWTONNU
*
*  this subroutine solves Keplers equation when the true anomaly is known.
*    The Mean and Eccentric, parabolic, or hyperbolic anomaly is also found.
*    The parabolic limit at 168 is arbitrary. The hyperbolic anomaly is also
*    limited. The hyperbolic sine is used because it's not double valued.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Ecc         - Eccentricity                   0.0D0 to
*    Nu          - True Anomaly                   -2Pi to 2Pi rad
*
*  Outputs       :
*    E0          - Eccentric Anomaly              0.0D0 to 2Pi rad       153.02D0
*    M           - Mean Anomaly                   0.0D0 to 2Pi rad       151.7425D0
*
*  Locals        :
*    E1          - Eccentric Anomaly, next value  rad
*    SinE        - Sine of E
*    CosE        - Cosine of E
*    Ktr         - Index
*
*  Coupling      :

*  References    :
*    Vallado       2007, 85, Alg 5
*
* ------------------------------------------------------------------------------

      SUBROUTINE NEWTONNU    ( Ecc, Nu, E0, M )
        IMPLICIT NONE
        REAL*8 Ecc, Nu, E0, M

* -----------------------------  Locals  ------------------------------
        REAL*8 SinE, CosE

        

        ! --------------------  Implementation   ----------------------
        E0= 999999.9D0
        M = 999999.9D0
        ! --------------------------- Circular ------------------------
        IF ( abs( Ecc ) .lt. 0.000001D0 ) THEN
            M = Nu
            E0= Nu 
          ELSE
            ! ---------------------- Elliptical -----------------------
            IF ( Ecc .lt. 0.999D0 ) THEN
                SinE= ( sqrt( 1.0D0-Ecc*Ecc ) * sin(Nu) ) /
     &                ( 1.0D0+Ecc*cos(Nu) )
                CosE= ( Ecc + cos(Nu) ) / ( 1.0D0 + Ecc*cos(Nu) )
                E0  = atan2( SinE, CosE )
                M   = E0 - Ecc*sin(E0) 
              ELSE
                ! -------------------- Hyperbolic  --------------------
                IF ( Ecc .gt. 1.0001D0 ) THEN
                    IF ( ((Ecc .gt. 1.0D0) .and. (abs(Nu)+0.00001D0
     &                     .lt. Pi-DACOS(1.0D0/Ecc)) ) ) THEN
                        SinE= ( sqrt( Ecc*Ecc-1.0D0 ) * sin(Nu) ) /
     &                        ( 1.0D0 + Ecc*cos(Nu) )
                        E0  = ASINH( SinE )
                        M   = Ecc*SINH(E0) - E0
                      ENDIF 
                  ELSE
                    ! ----------------- Parabolic ---------------------
                    IF ( abs(Nu) .lt. 168.0D0/57.29578D0 ) THEN
                        E0= DTAN( Nu*0.5D0 )
                        M = E0 + (E0*E0*E0)/3.0D0 
                      ENDIF
                  ENDIF
              ENDIF
          ENDIF

        IF ( Ecc .lt. 1.0D0 ) THEN
            M = DMOD( M, 2.0D0*Pi )
            IF ( M .lt. 0.0D0 ) THEN
                M= M + 2.0D0*Pi 
              ENDIF
            E0 = DMOD( E0, 2.0D0*Pi )
          ENDIF 
      RETURN
      END  ! end newtonnu
*
* ------------------------------------------------------------------------------
*
*                           SUBROUTINE KEPLER
*
*  this subroutine solves Keplers problem for orbit determination and returns a
*    future Geocentric Equatorial (IJK) position and velocity vector.  The
*    solution SUBROUTINE uses Universal variables.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*                                                               15 Jan 1993
*                     29 Nov 95 - needs repeat until - else it can't for 10.0D0TUs
*  Inputs          Description                    Range / Units
*    Ro          - IJK Position vector - initial  km
*    Vo          - IJK Velocity vector - initial  km / s
*    dtsec        - Length of time to propagate    s
*
*  OutPuts       :
*    R           - IJK Position vector            km
*    V           - IJK Velocity vector            km / s
*    Error       - Error flag                     'ok', ...
*
*  Locals        :
*    F           - f expression
*    G           - g expression
*    FDot        - f DOT expression
*    GDot        - g DOT expression
*    XOld        - Old Universal Variable X
*    XOldSqrd    - XOld squared
*    XNew        - New Universal Variable X
*    XNewSqrd    - XNew squared
*    ZNew        - New value of z
*    C2New       - C2(psi) FUNCTION
*    C3New       - C3(psi) FUNCTION
*    dtsec        - change in time                 s
*    TimeNew     - New time                       s
*    RDotV       - Result of Ro DOT Vo
*    A           - Semi .or. axis                 km
*    Alpha       - Reciprocol  1/a
*    SME         - Specific Mech Energy           km2 / s2
*    Period      - Time period for satellite      s
*    S           - Variable for parabolic case
*    W           - Variable for parabolic case
*    H           - Angular momentum vector
*    Temp        - Temporary REAL*8 value
*    i           - index
*
*  Coupling      :
*    MAG         - Magnitude of a vector
*    DOT         - DOT product of two vectors
*    COT         - Cotangent FUNCTION
*    FINDC2C3    - Find C2 and C3 functions
*    CROSS       - CROSS product of two vectors
*
*  References    :
*    Vallado       2007, 101, Alg 8, Ex 2-4
*
* ------------------------------------------------------------------------------

      SUBROUTINE KEPLER      ( Ro, Vo, dtsec, R, V, Error )
        IMPLICIT NONE
        REAL*8 Ro(3), Vo(3), dtsec, R(3), V(3), Dot
        CHARACTER*12 Error

* -----------------------------  Locals  ------------------------------
        INTEGER Ktr, i, NumIter
        REAL*8 H(3), F, G, FDot, GDot, Rval, XOld, XOldSqrd, XNew,
     &      XNewSqrd, ZNew, p, C2New, C3New, DtNew, RDotV, A,
     &      Alpha, SME, Period, S, W, Temp,
     &      magro, magvo, magh, magr

        
        

        ! --------------------  Implementation   ----------------------
        NumIter    =    35
        ! --------------------  Initialize values   -------------------
        ktr = 0
        xOld= 0.0D0
        ZNew= 0.0D0
        Error= 'ok' 

        IF ( abs( dtsec ) .gt. Small ) THEN
            magro = norm2( Ro )
            magvo = norm2( Vo )
            RDotV= dot_product( Ro, Vo )

            ! -------------  Find SME, Alpha, and A  ------------------
            SME= ( magvo*magvo*0.5D0 ) - ( mu/magro )
            Alpha= -SME*2.0D0 / mu

            IF ( abs( SME ) .gt. Small ) THEN
                A= -mu / ( 2.0D0*SME )
              ELSE
                A= Infinite
              ENDIF
            IF ( abs( Alpha ) .lt. Small ) THEN ! Parabola
                Alpha= 0.0D0
              ENDIF

            ! ------------   Setup initial guess for x  ---------------
            ! -----------------  Circle and Ellipse -------------------
            IF ( Alpha .ge. Small ) THEN
                Period= TwoPi * sqrt( abs(A)**3/mu )
                ! ------- Next IF needed for 2body multi-rev ----------
                IF ( abs( dtsec ) .gt. abs( Period ) ) THEN
                    dtsec= DMOD( dtsec, Period )
                  ENDIF
                IF ( abs(Alpha-1.0D0) .gt. Small ) THEN
                     XOld = sqrt(mu) * dtsec * Alpha
                  ELSE
                     ! 1st guess can't be too close. ie a circle, r=a
                     XOld= sqrt(mu) * dtsec*Alpha*0.97D0
                  ENDIF
              ELSE
                ! --------------------  Parabola  ---------------------
                IF ( abs( Alpha ) .lt. Small ) THEN
                    CALL CROSS( ro, vo, h )
                    magh = norm2(h)
                    p= magh*magh/mu
                    S= 0.5D0 * (HalfPi - atan( 3.0D0*sqrt( mu/
     &                         (p*p*p) )* dtsec ) )
                    W= atan( DTAN( S )**(1.0D0/3.0D0) )
                    XOld = sqrt(p) * ( 2.0D0*COT(2.0D0*W) )
                    Alpha= 0.0D0 
                  ELSE
                    ! ------------------  Hyperbola  ------------------
                    Temp= -2.0D0*mu*dtsec /
     &                   ( A*( RDotV + DSIGN(1.0D0, dtsec)*sqrt(-mu*A)*
     &                   (1.0D0-magro*Alpha) ) )
                    XOld= DSIGN(1.0D0, dtsec) * sqrt(-A) *DLOG(Temp)
                  ENDIF

              ENDIF
            Ktr= 1
            DtNew = -10.0D0
            DO WHILE ( (abs(DtNew-sqrt(mu)*dtsec).ge.Small).and.
     &                 (Ktr.lt.NumIter) )
                XOldSqrd = XOld*XOld 
                ZNew     = XOldSqrd * Alpha

                ! ------------- Find C2 and C3 functions --------------
                CALL FINDC2C3( ZNew, C2New, C3New )

                ! ------- Use a Newton iteration for New values -------
                DtNew= XOldSqrd*XOld*C3New + (RDotV/sqrt(mu))*XOldSqrd
     &                   *C2New + magro*XOld*( 1.0D0 - ZNew*C3New )
                Rval = XOldSqrd*C2New + (RDotV/sqrt(mu))*XOld*(1.0D0-
     &                   ZNew*C3New) + magro*( 1.0D0 - ZNew*C2New )

                ! ------------- Calculate New value for x -------------
                XNew = XOld + ( sqrt(mu)*dtsec - DtNew ) / Rval

* ------------------------------------------------------
*  Check If the orbit is an ellipse and xNew .gt. 2pi sqrt(a), the step
*  size must be changed.  This is accomplished by multiplying Rval
*  by 10.0D0.  NOTE !! 10.0D0 is arbitrary, but seems to produce good
*  results.  The idea is to keep XNew from increasing too rapidily.
* ------------------------------------------------------
*  including this doesn't work If you don't MOD the dtsec
*             IF ( ( A .gt. 0.0D0 ) .and. ( abs(XNew)>TwoPi*sqrt(A) ) .and. ( SME .lt. 0.0D0 ) ) THEN
*                 Dx= ( dtsec-DtNew ) / Rval  ! *7.0D0  * 10.0D0
*                 XNew = XOld + Dx / 7.0D0   ! /(1.0D0 + Dx)
*!   Alternate method to test various values of change
*XNew = XOld + ( dtsec-DtNew ) / ( RVal*10 chgamt  )
*               ENDIF 

                  Ktr = Ktr + 1
                XOld = XNew 
              ENDDO

            IF ( Ktr .ge. NumIter ) THEN
                Error= 'KNotConv'
c               Write(*,*) ' Not converged in ', NumIter:2,' iterations '
                DO i= 1 , 3
                    V(i)= 0.0D0
                    R(i)= V(i)
                  ENDDO
              ELSE
                ! --- Find position and velocity vectors at New time --
                XNewSqrd = XNew*XNew
                F = 1.0D0 - ( XNewSqrd*C2New / magro )
                G = dtsec - XNewSqrd*XNew*C3New/sqrt(mu)
                DO i= 1 , 3
                    R(i)= F*Ro(i) + G*Vo(i)
                  ENDDO
                magr = norm2( R )
                GDot = 1.0D0 - ( XNewSqrd*C2New / magr )
                FDot = ( sqrt(mu)*XNew / ( magro*magr ) ) *
     &                 ( ZNew*C3New-1.0D0 )
                DO i= 1 , 3
                    V(i)= FDot*Ro(i) + GDot*Vo(i)
                  ENDDO
                Temp= F*GDot - FDot*G 
                IF ( abs(Temp-1.0D0) .gt. 0.00001D0 ) THEN
                    Error= 'FandG'
                  ENDIF 
              ENDIF  ! IF (
          ELSE
            ! ----------- Set vectors to incoming since 0 time --------
            DO i=1, 3
                r(i)= ro(i)
                v(i)= vo(i)
              ENDDO
          ENDIF 

      RETURN
      END
*
* ------------------------------------------------------------------------------
*
*                           SUBROUTINE FINDTOF
*
*  this subroutine finds the time of flight given the initial position vectors,
*    Semi-parameter, and the sine and cosine values for the change in true
*    anomaly.  The result uses p-iteration theory to analytically find the result.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Ro          - Interceptor position vector    km
*    R           - TARGET position vector         km
*    p           - Semiparameter                  km
*
*  Outputs       :
*    Tof         - Time for transfer              s
*
*  Locals        :
*    SinDNu      - Sine of change in Nu           rad
*    CosDNu      - Cosine of change in Nu         rad
*    DeltaE      -
*    DeltaH      -
*    k           -
*    l           -
*    m           -
*    a           -
*    f           -
*    g           -
*    FDot        -
*    SinDeltaE   - Sine value
*    CosDeltaE   - Cosine value
*    RcrossR     - CROSS product of two positions
*
*  Coupling      :
*    CROSS       - CROSS product of two vectors

*  References    :
*    Vallado       2007, 134, Alg 11
*
* ------------------------------------------------------------------------------

      SUBROUTINE FINDTOF     ( Ro, R, p, Tof )
        IMPLICIT NONE
        REAL*8 Ro(3), R(3), p, Tof
* -----------------------------  Locals  ------------------------------
        REAL*8 RCrossR(3), CosDNu, SinDNu, c , s, alpha, DeltaE,
     &    DeltaH, DNu, k, l, m, a, f, g, FDot, SinDeltaE, CosDeltaE,
     &    magro, magr, magrcrossr

        

        ! --------------------  Implementation   ----------------------
        magro = norm2(ro)
        magr = norm2(r)
        CosDNu= dot_product(Ro, R)/(magro*magr)
        CALL CROSS( Ro, R, RCrossR )
        magrcrossr = norm2(rcrossr)
        SinDNu= magRCrossR/(magro*magr)

        k= magro * magr*( 1.0D0-CosDNu )
        l= magro + magr
        m= magro * magr*( 1.0D0+CosDNu )
        a= (m*k*p) / ((2.0D0*m-l*l)*p*p + 2.0D0*k*l*p - k*k) 

        ! ------  Use F and G series to find Velocity Vectors  --------
        F = 1.0D0 - ( magr/p )*(1.0D0-CosDNu)
        G = magro*magr*SinDNu/sqrt(mu*p)
        Alpha= 1.0D0/a 

        IF ( alpha .gt. Small ) THEN
            ! ------------------------ Elliptical ---------------------
            DNu  = atan2( SinDNu, CosDNu )
            FDot = sqrt(mu/p) * DTAN(DNu*0.5D0)*
     &              ( ((1.0D0-CosDNu)/p)-(1.0D0/magro)-(1.0D0/magr) )
            COSDeltaE= 1.0D0-(magro/a)*(1.0D0-f)
            SinDeltaE= (-magro*magr*FDot)/sqrt(mu*a)
            DeltaE   = atan2( SinDeltaE, CosDeltaE )
            Tof      = G + sqrt(a*a*a/mu)*(DeltaE-SinDeltaE)
          ELSE
            ! ------------------------ Hyperbolic ---------------------
            IF ( alpha .lt. -Small ) THEN
                DeltaH = ACOSH( 1.0D0-(magro/a)*(1.0D0-F) )
                Tof    = G + sqrt(-a*a*a/mu)*(SINH(DeltaH)-DeltaH)
              ELSE
                ! -------------------- Parabolic ----------------------
                DNu= atan2( SinDNu, CosDNu )
                c  = sqrt( magr*magr+magro*magro -
     &                     2.0D0*magr*magro*cos(DNu) )
                s  = (magro+magr+c ) * 0.5D0
                Tof= ( 2.0D0/3.0D0 ) * sqrt(s*s*s*0.5D0/mu) *
     &                     (1.0D0 -  ((s-c)/s)**1.5D0 )
              ENDIF 
          ENDIF

      RETURN
      END
*
* ------------------------------------------------------------------------------
*
*                           SUBROUTINE ijk2ll
*
*  These SUBROUTINEs convert a Geocentric Equatorial (IJK) position vector into
*    latitude and longitude.  Geodetic and Geocentric latitude are found.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    R           - IJK position vector            km
*    JD          - Julian Date                    days from 4713 BC
*
*  OutPuts       :
*    Latgc       - Geocentric Latitude            -Pi to Pi rad
*    Latgd       - Geodetic Latitude              -Pi to Pi rad
*    Lon         - Longitude (WEST -)             -2Pi to 2Pi rad
*    Hellp       - Height above the ellipsoid     km
*
*  Locals        :
*  Escobal:
*    Rc          - Range of SITE wrt earth center km
*    Height      - Height above earth wrt SITE    km
*    Alpha       - ANGLE from Iaxis to point, LST rad
*    OldDelta    - Previous value of DeltaLat     rad
*    DeltaLat    - Diff between Delta and
*                  Geocentric lat                 rad
*    Delta       - Declination ANGLE of R in IJK  rad
*    RSqrd       - Magnitude of r squared         km2
*    SinTemp     - Sine of Temp                   rad
*    c           -
*
*  Almanac:
*    Temp        - Diff between Geocentric/
*                  Geodetic lat                   rad
*    GST         - Greenwich SIDEREAL time        rad
*    SinTemp     - Sine of Temp                   rad
*    OldDelta    - Previous value of DeltaLat     rad
*    RtAsc       - Right ascension                rad
*    Decl        - Declination                    rad
*    i           - index
*
*  Borkowski:
*
*
*  Coupling      :
*    MAG         - Magnitude of a vector
*    GSTIME      - Greenwich SIDEREAL Time
*    gc2gd    - Converts between geocentric and geodetic latitude
*
*  References    :
*    Vallado       2007, 179, Alg 12 and Alg 13, Ex 3-3
*
* ------------------------------------------------------------------------------

      SUBROUTINE ijk2llA ( R, JD, Latgc, Latgd, Lon, Hellp )
        IMPLICIT NONE
        REAL*8 R(3), JD, Latgc, Latgd, Lon, Hellp

* -----------------------------  Locals  ------------------------------
        INTEGER i
        REAL*8 RtAsc, OldDelta, c, Decl, Temp, GST, SinTemp,  magr


        ! --------------------  Implementation   ----------------------
        magr = norm2( R )

        ! ----------------- Find Longitude value  ---------------------
        Temp = sqrt( R(1)*R(1) + R(2)*R(2) )
        IF ( abs( Temp ) .lt. Small ) THEN
            RtAsc= DSIGN(1.0D0, R(3))*Pi*0.5D0
          ELSE
            RtAsc= atan2( R(2) / Temp , R(1) / Temp )
          ENDIF
        GST  = GSTIME( JD ) 
        Lon  = RtAsc - GST 
        IF ( abs(Lon) .ge. Pi ) THEN ! Mod it ?
            IF ( Lon .lt. 0.0D0 ) THEN
                Lon= TwoPi + Lon
              ELSE
                Lon= Lon - TwoPi 
              ENDIF
          ENDIF
        Decl = asin( R(3) / magr )
        Latgd= Decl 

        ! ------------- Iterate to find Geodetic Latitude -------------
        i= 1 
        OldDelta = Latgd + 10.0D0

        DO WHILE ((abs(OldDelta-Latgd).ge.Small).and.(i.lt.10))
            OldDelta= Latgd 
            SinTemp = sin( Latgd ) 
            c       = rekm / (sqrt( 1.0D0-EESqrd*SinTemp*SinTemp ))
            Latgd= atan( (r(3)+c*EESqrd*SinTemp)/Temp )
            i = i + 1
          ENDDO

        Hellp   = (Temp/cos(Latgd)) - c

        CALL gc2gd( Latgc, 'FROM', Latgd )

      RETURN
      END



      SUBROUTINE ijk2llE ( R, JD, Latgc, Latgd, Lon, Hellp )
        IMPLICIT NONE
        REAL*8 R(3), JD, Latgc, Latgd, Lon, Hellp

* -----------------------------  Locals  ------------------------------
        INTEGER i
        Real*8 rsite, DeltaLat, RSqrd
        REAL*8 RtAsc, OldDelta, Decl,
     &         Temp, GST, SinTemp,  OneMinusE2,  magr
        CHARACTER Show

        
        

        ! --------------------  Implementation   ----------------------
        Show = 'N'

       ! -------------------  Initialize values   --------------------
        magr = norm2( R )
        OneMinuse2 = 1.0D0 - EeSqrd

       ! ---------------- Find Longitude value  ----------------------
        Temp = sqrt( R(1)*R(1) + R(2)*R(2) )
        IF ( abs( Temp ) .lt. Small ) THEN
            RtAsc= DSIGN(1.0D0, R(3))*Pi*0.5D0
          ELSE
            RtAsc= atan2( R(2) / Temp , R(1) / Temp )
          ENdif
        GST  = GSTIME( JD )
        Lon  = RtAsc - GST 

        IF ( abs(Lon) .ge. Pi ) THEN
            IF ( Lon .lt. 0.0D0 ) THEN
                Lon= TwoPi + Lon
              ELSE
                Lon= Lon - TwoPi 
              ENDIF
          ENDIF
       ! -------------- Set up initial latitude value  ---------------  
        Decl    = asin( R(3) / magr )
        Latgc= Decl 
        DeltaLat= 100.0D0 
        RSqrd   = magr**2

       ! ---- Iterate to find Geocentric .and. Geodetic Latitude  -----  
        i= 1 
        DO WHILE ( ( abs( OldDelta - DeltaLat ) .ge. Small ) .and.
     &             ( i .lt. 10 ))
            OldDelta = DeltaLat 
            rsite    = sqrt( OneMinuse2 / (1.0D0 -
     &                 EeSqrd*(cos(Latgc))**2 ) )
            Latgd = atan( DTAN(Latgc) / OneMinuse2 ) 
            Temp     = Latgd-Latgc 
            SinTemp  = sin( Temp ) 
            Hellp    = sqrt( RSqrd - rsite*rsite*SinTemp*SinTemp ) -
     &                 rsite*cos(Temp)
            DeltaLat = asin( Hellp*SinTemp / magr )
            Latgc = Decl - DeltaLat 
            i = i + 1
            IF ( Show .eq. 'Y' ) THEN
                write(*, *) 'E loops gc gd ', Latgc*57.29578D0,
     &                     Latgd*57.29578D0
              ENDIF
          ENDDO

        IF ( i .ge. 10 ) THEN
            Write(*, *) 'ijk2ll did NOT converge '
          ENDIF

      RETURN
      END

* ------------------------------- Borkowski method  --------------------------
      SUBROUTINE ijk2llB ( R , JD, Latgc, Latgd, Lon, Hellp )
        IMPLICIT NONE
        REAL*8 JD, R(3), Latgc, Latgd, Lon, Hellp

* -----------------------------  Locals  ------------------------------
        Real*8 a, b, RtAsc, sqrtp, third, e, f, p, q, d, nu, g, t,
     &         aTemp, Temp, GST

        

        ! --------------------  Implementation   ----------------------

        ! ---------------- Find Longitude value  ----------------------
        Temp = sqrt( R(1)*R(1) + R(2)*R(2) )
        IF ( abs( Temp ) .lt. Small ) THEN
            RtAsc= DSIGN(1.0D0, R(3))*Pi*0.5D0
          ELSE
            RtAsc= atan2( R(2) / Temp , R(1) / Temp )
          ENDIF
        GST  = GSTIME( JD )
        Lon  = RtAsc - GST 
        IF ( abs(Lon) .ge. Pi ) THEN
            IF ( Lon .lt. 0.0D0 ) THEN
                Lon= TwoPi + Lon
              ELSE
                Lon= Lon - TwoPi
              ENDIF
          ENDIF

        a= 1.0D0 
        b= DSIGN(1.0D0, r(3))*6356.75160056D0/6378.137D0
       ! -------------- Set up initial latitude value  ---------------  
        aTemp= 1.0D0/(a*Temp) 
        e= (b*r(3)-a*a+b*b)*atemp
        f= (b*r(3)+a*a-b*b)*atemp
        third= 1.0D0/3.0D0 
        p= 4.0D0*Third*(e*f + 1.0D0 ) 
        q= 2.0D0*(e*e - f*f) 
        d= p*p*p + q*q 

        IF ( d .gt. 0.0D0 ) THEN
            nu= (sqrt(d)-q)**third - (sqrt(d)+q)**third
          ELSE
            SqrtP= sqrt(-p)
            nu= 2.0D0*SqrtP*cos( third*DACOS(q/(p*SqrtP)) ) 
          ENDIF 
        g= 0.5D0*(sqrt(e*e + nu) + e) 
        t= sqrt(g*g + (f-nu*g)/(2.0D0*g-e)) - g 

        Latgd= atan(a*(1.0D0-t*t)/(2.0D0*b*t)) 
        hellp= (temp-a*t)*cos( Latgd) + (r(3)-b)*sin(Latgd)

        CALL gc2gd( Latgc, 'FROM', Latgd )
      RETURN
      END
*
* ------------------------------------------------------------------------------
*
*                           SUBROUTINE GC2GD
*
*  this subroutine converts from Geodetic to Geocentric Latitude for positions
*    on the surface of the Earth.  Notice that (1-f) squared = 1-eSqrd.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Latgd       - Geodetic Latitude              -Pi to Pi rad
*    Direction   - Which set of vars to output    FROM  TOO
*
*  Outputs       :
*    Latgc       - Geocentric Latitude            -Pi to Pi rad
*
*  Locals        :
*    None.
*
*  Coupling      :
*    None.
*
*  References    :
*    Vallado       2007, 148, Eq 3-11
*
* ------------------------------------------------------------------------------

      SUBROUTINE gc2gd       ( Latgc, Direction, Latgd )
        IMPLICIT NONE
        REAL*8 Latgc, Latgd
        CHARACTER*4 Direction

        

        ! --------------------  Implementation   ----------------------
        IF ( Direction .eq. 'FROM' ) THEN
            Latgc= atan( (1.0D0 - EESqrd)*DTAN(Latgd) )
          ELSE
            Latgd= atan( DTAN(Latgc)/(1.0D0 - EESqrd) )
          ENDIF
      RETURN
      END


* ------------------------------------------------------------------------------
*
*                           SUBROUTINE SIGHT
*
*  this subroutine takes the position vectors of two satellites and determines
*    If there is line-of-SIGHT between the two satellites.  An oblate Earth
*    with radius of 1 ER is assumed.  The process forms the equation of
*    a line between the two vectors.  Differentiating and setting to zero finds
*    the minimum value, and when plugged back into the original line equation,
*    gives the minimum distance.  The parameter TMin is allowed to range from
*    0.0D0 to 1.0D0.  Scale the K-component to account for oblate Earth because it's
*    the only qunatity that changes.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    R1          - Position vector of the 1st sat km
*    R2          - Position vector of the 2nd sat km
*    WhichKind   - Spherical .or. Ellipsoidal Earth 'S', 'E'*default
*
*  Outputs       :
*    LOS         - Line of SIGHT                  'YES', 'NO '
*
*  Locals        :
*    TR1         - Scaled R1 vector               km
*    TR2         - Scaled R2 vector               km
*    ADotB       - DOT product of a DOT b
*    TMin        - Minimum value of t from a to b
*    DistSqrd    - minute Distance squared to Earth  km
*    ASqrd       - Magnitude of A squared
*    BSqrd       - Magnitude of B squared
*
*  Coupling:
*    DOT         - DOT product of two vectors
*
*  References    :
*    Vallado       2007, 310, Alg 35, Ex 5-3
* ------------------------------------------------------------------------------

      SUBROUTINE SIGHT       ( R1, R2, WhichKind, LOS )

        REAL*8 R1(3), R2(3)
        INTEGER i
        CHARACTER WhichKind
        CHARACTER(3) LOS

* -----------------------------  Locals  ------------------------------
        REAL*8 TR1(3), TR2(3), ADotB, TMin, DistSqrd, ASqrd,
     &         BSqrd, Temp,  magtr1, magtr2

        

        ! --------------------  Implementation   ----------------------
        DO i=1 , 3
            TR1(i)= R1(i)
            TR2(i)= R2(i)
          ENDDO
        ! --------------------- Scale z component ---------------------
        IF ( WhichKind .eq. 'E' ) THEN
            Temp= 1.0D0/sqrt(1.0D0-EESqrd)
          ELSE
            Temp= 1.0D0
          ENDIF
        TR1(3)= TR1(3)*Temp
        TR2(3)= TR2(3)*Temp
        magtr1 = norm2(tr1)
        magtr2 = norm2(tr2)
        BSqrd= magTR2**2
        ASqrd= magTR1**2
        ADotB= dot_product( TR1, TR2 )
        ! ---------------------- Find TMin ----------------------------
        DistSqrd= 0.0D0
        IF ( abs(ASqrd + BSqrd - 2.0D0*ADotB) .lt. 0.0001D0 ) THEN
            TMin= 0.0D0
          ELSE
            TMin = ( ASqrd - ADotB ) / ( ASqrd + BSqrd - 2.0D0*ADotB )
          ENDIF
        ! ----------------------- Check LOS ---------------------------
        IF ( (TMin .lt. 0.0D0) .or. (TMin .gt. 1.0D0) ) THEN
            LOS= 'YES'
          ELSE
            DistSqrd= ( (1.0D0-TMin)*ASqrd + ADotB*TMin ) / rekm**2
            IF ( DistSqrd .gt. rekm ) THEN
                LOS= 'YES'
              ELSE
                LOS= 'NO '
              ENDIF
          ENDIF
      RETURN
      END

* ------------------------------------------------------------------------------
*
*                           SUBROUTINE SUN
*
*  this subroutine calculates the Geocentric Equatorial position vector
*    the SUN given the Julian Date.  This is the low precision formula .and.
*    is valid for years from 1950 to 2050.  Accuaracy of apparent coordinates
*    is 0.01D0 degrees.  Notice many of the calcuLations are performed in
*    degrees, and are not changed until Later.  This is due to the fact that
*    the Almanac uses degrees exclusively in their formuLations.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    JD          - Julian Date                    days from 4713 BC
*
*  Outputs       :
*    RSun        - IJK Position vector of the SUN AU
*    RtAsc       - Right Ascension                rad
*    Decl        - Declination                    rad
*
*  Locals        :
*    MeanLong    - Mean Longitude
*    MeanAnomaly - Mean anomaly
*    EclpLong    - Ecliptic Longitude
*    Obliquity   - Mean Obliquity of the Ecliptic
*    TUT1        - Julian Centuries of UT1 from
*                  Jan 1, 2000 12h
*    TTDB        - Julian Centuries of TDB from
*                  Jan 1, 2000 12h
*    Hr          - Hours                          0 .. 24              10
*    minute         - MiNutes                        0 .. 59              15
*    SEC         - Seconds                        0.0D0 .. 59.99D0         30.00D0
*    Temp        - Temporary variable
*    deg         - Degrees
*
*  Coupling      :
*    None.
*
*  References    :
*    Vallado       2007, 281, Alg 29, Ex 5-1
*
* ------------------------------------------------------------------------------

      SUBROUTINE SUN         ( JD, RSun, RtAsc, Decl )
        IMPLICIT NONE
        REAL*8 JD, RSun(3), RtAsc, Decl
* -----------------------------  Locals  ------------------------------
        REAL*8 MeanLong, MeanAnomaly,
     &        EclpLong, Obliquity, TUT1, TTDB, magrsun

        

        ! --------------------  Implementation   ----------------------
        ! -------------------  Initialize values   --------------------
        TUT1= ( JD - 2451545.0D0 )/ 36525.0D0 

        MeanLong= 280.4606184D0 + 36000.77005361D0*TUT1 
        MeanLong= DMOD( MeanLong, 360.0D0 )  !deg

        TTDB= TUT1 
        MeanAnomaly= 357.5277233D0 + 35999.05034D0*TTDB 
        MeanAnomaly= DMOD( MeanAnomaly*Deg2Rad, TwoPi )  !rad
        IF ( MeanAnomaly .lt. 0.0D0 ) THEN
            MeanAnomaly= TwoPi + MeanAnomaly
          ENDIF

        EclpLong= MeanLong + 1.914666471D0*sin(MeanAnomaly)
     &              + 0.019994643D0*sin(2.0D0*MeanAnomaly) !deg

        Obliquity= 23.439291D0 - 0.0130042D0*TTDB  !deg

        MeanLong = MeanLong*Deg2Rad 
        IF ( MeanLong .lt. 0.0D0 ) THEN
            MeanLong= TwoPi + MeanLong
          ENDIF
        EclpLong = EclpLong *Deg2Rad 
        Obliquity= Obliquity *Deg2Rad 

        ! ------- Find magnitude of SUN vector, ) THEN components -----
        magRSun= 1.000140612D0 - 0.016708617D0*cos( MeanAnomaly )
     &                         - 0.000139589D0*cos( 2.0D0*MeanAnomaly )    ! in AU's

        RSun(1)= magRSun*cos( EclpLong )
        RSun(2)= magRSun*cos(Obliquity)*sin(EclpLong)
        RSun(3)= magRSun*sin(Obliquity)*sin(EclpLong)

        RtAsc= atan( cos(Obliquity)*DTAN(EclpLong) )
        ! --- Check that RtAsc is in the same quadrant as EclpLong ----
        IF ( EclpLong .lt. 0.0D0 ) THEN
            EclpLong= EclpLong + TwoPi    ! make sure it's in 0 to 2pi range
          ENDIF
        IF ( abs( EclpLong-RtAsc ) .gt. Pi*0.5D0 ) THEN
            RtAsc= RtAsc + 0.5D0*Pi*DNINT( (EclpLong-RtAsc)/(0.5D0*Pi))
          ENDIF
        Decl = asin( sin(Obliquity)*sin(EclpLong) )

      RETURN
      END
*
*  References    :
*    Vallado       2007, 311, Eq 5-9

      SUBROUTINE SunIll   ( JD, Lat, Lon, SunIllum, SunAz, SunEl )
        IMPLICIT NONE
        REAL*8 SunIllum, Lat, Lon, JD, SunAz, SunEl

* -----------------------------  Locals  ------------------------------
        Real*8 RSun(3)
        Real*8 lst, gst, x, LHA, sinv, cosv
        Real*8 l0, l1, l2, l3, sRtAsc, sdecl

        

        ! --------------------  Implementation   ----------------------
        CALL SUN( JD, RSun, sRtAsc, sDecl ) ! AU's needed for Sun ill

        CALL LSTIME( Lon, JD, LST, GST )

        LHA = LST - sRtAsc

        SunEl  = asin( sin(sDecl)*sin(Lat) +
     &           cos(sDecl)*cos(Lat)*cos(LHA) )

        Sinv= -sin(LHA)*cos(sDecl)*cos(Lat)/(cos(SunEl)*cos(Lat))
        Cosv= ( sin(sDecl)-sin(SunEl)*sin(Lat) )/
     &        ( cos(SunEl)*cos(Lat) )
        SunAz  = atan2( Sinv, Cosv )

        SunEl= SunEl/Deg2Rad

        IF (SunEl .gt. -18.01D0) THEN
            x= SunEl/90.0D0

         IF (SunEl .ge. 20) THEN
             l0=  3.74
             l1=  3.97
             l2= -4.07
             l3=  1.47
           ELSEIF ((SunEl .ge. 5.0).and.(SunEl .lt. 20.0)) THEN
                 l0=   3.05
                 l1=  13.28
                 l2= -45.98
                 l3=  64.33
               ELSEIF ((SunEl .ge. -0.8).and.(SunEl .lt. 5.0)) THEN
                     l0=    2.88
                     l1=   22.26
                     l2= -207.64
                     l3= 1034.30
                   ELSEIF ((SunEl .ge. -5.0).and.(SunEl .lt. -0.8))
     &                THEN
                         l0=    2.88
                         l1=   21.81
                         l2= -258.11
                         l3= -858.36
                       ELSEIF ((SunEl .ge. -12.0).and.(SunEl .lt.
     &                          -5.0)) THEN
                             l0=    2.70
                             l1=   12.17
                             l2= -431.69
                             l3=-1899.83
                           ELSEIF ((SunEl .ge. -18.0).and.(SunEl .lt.
     &                              -12.0)) THEN
                                 l0=   13.84
                                 l1=  262.72
                                 l2= 1447.42
                                 l3= 2797.93
                               ELSE
                                 l0= 0.0
                                 l1= 0.0
                                 l2= 0.0
                                 l3= 0.0
                               ENDIF

         l1= l0 + l1*x + l2*x*x + l3*x*x*x
         SunIllum= 10.0** l1
         IF ((SunIllum .lt. -1D+36).or.(SunIllum .gt. 999.999D0)) THEN
             SunIllum= 0.0
           ENDIF
       ELSE
         SunIllum= 0.0D0
       ENDIF

      RETURN
      END
*
* ------------------------------------------------------------------------------
*
*                           SUBROUTINE MOON
*
*  this subroutine calculates the Geocentric Equatorial (IJK) position vector
*    for the MOON given the Julian Date.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    JD          - Julian Date                    days from 4713 BC
*
*  Outputs       :
*    RMoon       - IJK Position vector of MOON    ER
*    RtAsc       - Right Ascension                rad
*    Decl        - Declination                    rad
*
*  Locals        :
*    EclpLong    - Ecliptic Longitude
*    EclpLat     - Eclpitic Latitude
*    HzParal     - Horizontal Parallax
*    l           - Geocentric Direction Cosines
*    m           -             "     "
*    n           -             "     "
*    TTDB        - Julian Centuries of TDB from
*                  Jan 1, 2000 12h
*    Hr          - Hours                          0 .. 24
*    minute         - MiNutes                        0 .. 59
*    SEC         - Seconds                        0.0D0 .. 59.99D0
*    deg         - Degrees
*
*  Coupling      :
*    None.
*
*  References    :
*    Vallado       2007, 290, Alg 31, Ex 5-3
*
* ------------------------------------------------------------------------------
*
      SUBROUTINE MOON        ( JD, RMoon, RtAsc, Decl )
        IMPLICIT NONE
        REAL*8 JD, RMoon(3), RtAsc, Decl

* -----------------------------  Locals  ------------------------------
        REAL*8 TTDB, l, m, n, Obliquity, magrmoon, EclpLong, EclpLat,
     &         HzParal

        

        ! --------------------  Implementation   ----------------------
        TTDB = ( JD - 2451545.0D0 ) / 36525.0D0

        EclpLong= 218.32D0 + 481267.883D0*TTDB
     &              + 6.29D0*sin( (134.9D0+477198.85D0*TTDB)*Deg2Rad )
     &              - 1.27D0*sin( (259.2D0-413335.38D0*TTDB)*Deg2Rad )
     &              + 0.66D0*sin( (235.7D0+890534.23D0*TTDB)*Deg2Rad )
     &              + 0.21D0*sin( (269.9D0+954397.70D0*TTDB)*Deg2Rad )
     &              - 0.19D0*sin( (357.5D0+ 35999.05D0*TTDB)*Deg2Rad )
     &              - 0.11D0*sin( (186.6D0+966404.05D0*TTDB)*Deg2Rad )  ! Deg

        EclpLat =   5.13D0*sin( ( 93.3D0+483202.03D0*TTDB)*Deg2Rad )
     &              + 0.28D0*sin( (228.2D0+960400.87D0*TTDB)*Deg2Rad )
     &              - 0.28D0*sin( (318.3D0+  6003.18D0*TTDB)*Deg2Rad )
     &              - 0.17D0*sin( (217.6D0-407332.20D0*TTDB)*Deg2Rad )  ! Deg

        HzParal =  0.9508D0 + 0.0518D0*cos( (134.9D0+477198.85D0*TTDB)
     &              *Deg2Rad )
     &            + 0.0095D0*cos( (259.2D0-413335.38D0*TTDB)*Deg2Rad )
     &            + 0.0078D0*cos( (235.7D0+890534.23D0*TTDB)*Deg2Rad )
     &            + 0.0028D0*cos( (269.9D0+954397.70D0*TTDB)*Deg2Rad )  ! Deg

        EclpLong = DMOD( EclpLong*Deg2Rad, TwoPi )
        EclpLat  = DMOD( EclpLat*Deg2Rad, TwoPi )
        HzParal  = DMOD( HzParal*Deg2Rad, TwoPi )

        Obliquity= 23.439291D0 - 0.0130042D0*TTDB  !deg
        Obliquity= Obliquity *Deg2Rad

        ! ------------ Find the geocentric direction cosines ----------
        l= cos( EclpLat ) * cos( EclpLong )
        m= cos(Obliquity)*cos(EclpLat)*sin(EclpLong)
     &       - sin(Obliquity)*sin(EclpLat)
        n= sin(Obliquity)*cos(EclpLat)*sin(EclpLong)
     &       + cos(Obliquity)*sin(EclpLat)

        ! ------------- Calculate MOON position vector ----------------
        magRMoon= 1.0D0/sin( HzParal )
        RMoon(1)= magRMoon*l
        RMoon(2)= magRMoon*m
        RMoon(3)= magRMoon*n

        ! -------------- Find Rt Ascension and Declination ------------
        RtAsc= atan2( m, l )
        Decl = asin( n )

      RETURN
      END
*
*  References    :
*    Vallado       2007, 311, Eq 5-9

      SUBROUTINE MoonIll     ( MoonEl, f, MoonIllum )
        IMPLICIT NONE
        REAL*8 MoonEl, f, MoonIllum

* -----------------------------  Locals  ------------------------------
        REAL*8  x, l0, l1, l2, l3

        ! --------------------  Implementation   ----------------------
        x= MoonEl/90.0D0
c        g= 1.0

       IF (MoonEl .ge. 20) THEN
         l0= -1.95
         l1=  4.06
         l2= -4.24
         l3=  1.56
        ELSEIF ((MoonEl .ge. 5.0).and.(MoonEl .lt. 20.0)) THEN
             l0=  -2.58
             l1=  12.58
             l2= -42.58
             l3=  59.06
           ELSEIF ((MoonEl .gt. -0.8).and.(MoonEl .lt. 5.0)) THEN
                 l0=   -2.79
                 l1=   24.27
                 l2= -252.95
                 l3= 1321.29
               ELSE
                 l0= 0.0
                 l1= 0.0
                 l2= 0.0
                 l3= 0.0
                 f= 0.0
c                 g= 0.0
               ENDIF

       l1= l0 + l1*x + l2*x*x + l3*x*x*x
       l2= (-0.00868D0*f - 2.2D-9*f*f*f*f)

*       HzParal =   0.9508 + 0.0518*cos( (134.9+477198.85*TTDB)*Deg2Rad )
*                + 0.0095*cos( (259.2-413335.38*TTDB)*Deg2Rad )
*                + 0.0078*cos( (235.7+890534.23*TTDB)*Deg2Rad )
*                + 0.0028*cos( (269.9+954397.70*TTDB)*Deg2Rad )   { Deg }
*       HzParal  = REALMOD( HzParal*Deg2Rad, TwoPi )
*       l3= (2.0* POWER(10.0, (HzParal*rad / 0.951))*g ) { use g to eliminate neg el passes }

        MoonIllum= 10.0D0 ** ( l1 + l2 )
        IF ((MoonIllum .lt. -1d+36).or.(MoonIllum .gt. 0.999D0)) THEN
            MoonIllum= 0.0D0
          ENDIF
       RETURN
       END
*
* ------------------------------------------------------------------------------
*
*                           SUBROUTINE LIGHT
*
*  this subroutine determines If a spacecraft is sunlit .or. in the dark at a
*    particular time.  An oblate Earth and cylindrical shadow is assumed.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    R           - Position vector of sat         km
*    JD          - Julian Date at desired time    Days from 4713 BC
*    WhichKind   - Spherical .or. Ellipsoidal Earth 'S', 'E'*default
*
*  OutPuts       :
*    Vis         - Visibility Flag                'YES', 'NO '
*
*  Locals        :
*    RtAsc       - Suns Right ascension           rad
*    Decl        - Suns Declination               rad
*    RSun        - SUN vector                     AU
*    AUER        - Conversion from AU to ER
*
*  Coupling      :
*    SUN         - Position vector of SUN
*    SIGHT       - Does Line-of-SIGHT exist beteen vectors
*
*  References    :
*    Vallado       2007, 310, Alg 35, Ex 5-6
*
* ------------------------------------------------------------------------------

      SUBROUTINE LIGHT       ( R, JD, WhichKind, LIT )

        REAL*8 R(3), JD
        CHARACTER WhichKind
        character(3) Lit
* -----------------------------  Locals  ------------------------------
        REAL*8 RSun(3), RtAsc, Decl

        

        ! --------------------  Implementation   ----------------------

        CALL SUN( JD, RSun, RtAsc, Decl )
        RSun = AUER * RSun

        ! ------------ Is the satellite in the shadow? ----------------
        CALL SIGHT( RSun, R, WhichKind, Lit )

      RETURN
      END
*
* ------------------------------------------------------------------------------
*
*                           SUBROUTINE CHECKHITEARTH
*
*  this subroutine checks to see If the trajectory hits the earth during the
*    transfer.  The first check determines If the satellite is initially
*    heading towards perigee, and finally heading away from perigee.  IF (
*    this is the case, the radius of perigee is calculated.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    RInt        - Initial Position vector of Int km
*    V1t         - Initial Velocity vector of trnskm/TU
*    RTgt        - Initial Position vector of Tgt km
*    V2t         - Final Velocity vector of trns  km/TU
*
*  Outputs       :
*    HitEarth    - Is Earth was impacted          'Y' 'N'
*
*  Locals        :
*    SME         - Specific mechanical energy
*    rp          - Radius of Perigee              km
*    TransA      - Semi-.or. axis of transfer     km
*    TransE      - Eccentricity of transfer
*    TransP      - Semi-paramater of transfer     km
*    HBar        - Angular momentum vector of
*                  transfer orbit
*
*  Coupling      :
*    DOT         - DOT product of vectors
*    MAG         - Magnitude of a vector
*    CROSS       - CALL CROSS product of vectors
*
*  References    :
*    Vallado       2007, 500, Alg 60
*
* ------------------------------------------------------------------------------

      SUBROUTINE CHECKHITEARTH ( Rint, V1t, Rtgt, V2t, HitEarth )
        IMPLICIT NONE
        REAL*8 RInt(3), V1t(3), RTgt(3), V2t(3)
        CHARACTER HitEarth

* -----------------------------  Locals  ------------------------------
        REAL*8 HBar(3), SME, rp, TransP, TransA, TransE,
     &         magrint, magv1t, maghbar

        ! --------------------  Implementation   ----------------------
        HitEarth= 'N'

        ! ---------- Find If trajectory intersects Earth --------------
        IF ((dot_product(Rint,V1t).lt.0.0D0).and.
     &      (dot_product(RTgt,V2t).gt.0.0D0)) THEN

            ! ---------------  Find H N and E vectors   ---------------
            CALL CROSS( RInt, V1t, HBar )
            magrint = norm2( rint )
            magv1t  = norm2( v1t )
            maghbar = norm2( HBar )

            IF ( maghbar .gt. 0.00001D0 ) THEN
                ! ---------  Find a e and semi-Latus rectum   ---------
                SME    = magV1t**2*0.5D0 - ( 1.0D0/magRInt )
                TransP = maghbar*maghbar
                TransE = 1.0D0
                IF ( abs( SME ) .gt. 0.00001D0 ) THEN
                    TransA= -1.0D0 / (2.0D0*SME)
                    TransE= sqrt( (TransA - TransP)/TransA )
                    rp= TransA*(1.0D0-TransE) 
                  ELSE
                    rp= TransP*0.5D0   ! Parabola
                  ENDIF

                IF ( abs( rp ) .lt. 1.0D0 ) THEN
                    HitEarth= 'Y' 
                  ENDIF
              ELSE
                Write(*, *) 'The orbit does not exist '
              ENDIF
          ENDIF

      RETURN
      END
*
* ------------------------------------------------------------------------------
*
*                           SUBROUTINE SATFOV
*
*  this subroutine finds parameters reLating to a satellite's FOV.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Incl        - Inclination                    rad
*    Az          - Azimuth                        rad
*    SLatgd      - Geodetic Latitude of sat       rad
*    SLon        - Longitude of sat               rad
*    SAlt        - Altitudeof satellite           km
*    TFOV        - Total field of view            rad
*    EtaCtr      - Ctr where sensor looks         rad
*
*  Outputs       :
*    FovMax      - Maximum field of view          rad
*    TotalRng    -
*    RhoMax      -
*    RhoMin      -
*    TgtLat      -
*    TgtLon      -
*
*  Locals        :
*    r           -
*    etaHopriz   -
*    RhoHoriz    -
*    gamma       -
*    rho         -
*    FovMin      -
*    Lat         -
*    Lon         -
*    MaxLat      -
*    MinLKat     -
*    i           - Index
*
*  Coupling      :
*    PATH        - Finds tgt location given initial location, range, and az
*
*  References    :
*    Vallado       2007, 845, Eq 11-8 to Eq 11-13, Ex 11-1
*
* ------------------------------------------------------------------------------
*
      SUBROUTINE SATFOV      ( Incl, Az, SLatgd, SLon, SAlt, tFOV,
     &                         EtaCtr, FovMax, TotalRng, RhoMax, RhoMin,
     &                         TgtLat, TgtLon )
        IMPLICIT NONE
        REAL*8 Incl, Az, SLatgd, SLon, SAlt, tFOV, EtaCtr, FovMAx,
     &     TotalRng, RhoMax, RhoMin, TgtLat, TgtLon
* -----------------------------  Locals  ------------------------------
        INTEGER i
        REAL*8 r, EtaHoriz, rhoHoriz, gamma, rho, FovMin, Lat,
     &     Lon, maxLat, minLat

        

        ! --------------------  Implementation   ----------------------
        ! ------- Find satellite parameters and limiting cases --------
        r       = 1.0D0 + SAlt 
        EtaHoriz= asin(1.0D0/r) 
        RhoHoriz= r*cos(EtaHoriz) 

        ! ---------------- Find Ground range ANGLE --------------------
        FovMax= tFOV*0.5D0 + EtaCtr 
        Gamma = Pi - asin( r*sin(FovMax) )   ! must use larger ANGLE
        Rho   = cos( gamma ) + r*cos(FovMax) 
        RhoMax= asin( Rho*sin(FovMax) ) 

        ! -------- for minimum, If the sensor looks off axis ----------
        IF ( abs(EtaCtr) .gt. 0.00001D0 ) THEN
            FovMin  = EtaCtr - tFOV*0.5D0
            Gamma   = Pi - asin( r*sin(FovMin) )  ! use larger
            Rho     = cos( gamma ) + r*cos(FovMin) 
            RhoMin  = asin( Rho*sin(FovMin) ) 
            TotalRng= RhoMax - RhoMin 
          ELSE
            ! --------------------- Nadir pointing --------------------
            FovMin  = 0.0D0
            RhoMin  = 0.0D0 
            TotalRng= 2.0D0*RhoMax  ! equal sided
          ENDIF

        ! -------------- Find location of center of FOV ---------------
        IF ( abs(EtaCtr) .gt. 0.00001D0 ) THEN
            CALL PATH( SLatgd, SLon, RhoMin + TotalRng*0.5D0, Az,
     &                 Lat, Lon )
          ELSE
            Lat= SLatgd
            Lon= SLon 
          ENDIF 

        ! ----- Loop around the New circle with the sensor range ------
        DO i= 0 , 72
            Az= i*5.0D0/Rad2Deg
            CALL PATH( Lat, Lon, TotalRng*0.5D0, Az,  TgtLat, TgtLon )
            IF ( i .eq. 0 ) THEN
                MaxLat= TgtLat
              ENDIF 
            IF ( i .eq. 36 ) THEN
                MinLat= TgtLat
              ENDIF 
          ENDDO

      RETURN
      END
*
* ------------------------------------------------------------------------------
*
*                           SUBROUTINE RNGAZ
*
*  this subroutine calculates the Range and Azimuth between two specified
*    ground points on a spherical Earth.  Notice the range will ALWAYS be
*    within the range of values listed since you do not know the direction of
*    firing, Long .or. short.  The SUBROUTINE will calculate rotating Earth ranges
*    If the Tof is passed in other than 0.0D0. Range is calulated in rad .and.
*    converted to ER by s = rO, but the radius of the Earth = 1 ER, so it's
*    s = O.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    LLat        - Start Geocentric Latitude      -Pi/2 to  Pi/2 rad
*    LLon        - Start Longitude (WEST -)       0.0D0 to 2Pi rad
*    TLat        - ENDIF Geocentric Latitude        -Pi/2 to  Pi/2 rad
*    TLon        - ENDIF Longitude (WEST -)         0.0D0 to 2Pi rad
*    Tof         - Time of Flight If ICBM, .or. 0.0D0 TU
*
*  OutPuts       :
*    Range       - Range between points           ER
*    Az          - Azimuth                        0.0D0 to 2Pi rad
*
*  Locals        :
*    None.
*
*  Coupling      :
*    None.
*
*  References    :
*    Vallado       2007, 843, Eq 11-3, Eq 11-4, Eq 11-5
*
* ------------------------------------------------------------------------------

      SUBROUTINE RNGAZ       ( LLat, LLon, TLat, TLon, Tof, Range, Az )
        IMPLICIT NONE
        REAL*8 LLat, LLon, TLat, TLon, Tof, Range, Az

        
        

        ! --------------------  Implementation   ----------------------
        Range= DACOS( sin(LLat)*sin(TLat) +
     &         cos(LLat)*cos(TLat)*cos(TLon-LLon + OmegaEarth*Tof) )

        ! ------ Check If the Range is 0 .or. half the earth  ---------
        IF ( abs( sin(Range)*cos(LLat) ) .lt. Small ) THEN
            IF ( abs( Range - Pi ) .lt. Small ) THEN
                Az= Pi
              ELSE
                Az= 0.0D0
              ENDIF
          ELSE
            Az= DACOS( ( sin(TLat) - cos(Range) * sin(LLat)) /
     &                 ( sin(Range) * cos(LLat)) )
          ENDIF

        ! ------ Check If the Azimuth is grt than Pi ( 180deg ) -------
        IF ( sin( TLon - LLon + OmegaEarth*Tof ) .lt. 0.0D0 ) THEN
            Az= TwoPi - Az
          ENDIF
      RETURN
      END
*
* ------------------------------------------------------------------------------
*
*                           SUBROUTINE PATH
*
*  this subroutine determines the ENDIF position for a given range and azimuth
*    from a given point.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    LLat        - Start Geocentric Latitude      -Pi/2 to  Pi/2 rad
*    LLon        - Start Longitude (WEST -)       0.0D0 to 2Pi rad
*    Range       - Range between points           km
*    Az          - Azimuth                        0.0D0 to 2Pi rad
*
*  OutPuts       :
*    TLat        - ENDIF Geocentric Latitude        -Pi/2 to  Pi/2 rad
*    TLon        - ENDIF Longitude (WEST -)         0.0D0 to 2Pi rad
*
*  Locals        :
*    SinDeltaN   - Sine of Delta N                rad
*    CosDeltaN   - Cosine of Delta N              rad
*    DeltaN      - ANGLE between the two points   rad
*
*  Coupling      :
*    None.
*
*  References    :
*    Vallado       2007, 843, Eq 11-6, Eq 11-7
*
* ------------------------------------------------------------------------------

      SUBROUTINE PATH        ( LLat, LLon, Range, Az, TLat, TLon )
        IMPLICIT NONE
        REAL*8 LLat, LLon, Range, Az, TLat, TLon
* -----------------------------  Locals  ------------------------------
        REAL*8 SinDN, CosDN, DeltaN

        

        ! --------------------  Implementation   ----------------------
        Az= DMOD( Az, TwoPi )
        IF ( LLon .lt. 0.0D0 ) THEN
            LLon= TwoPi + LLon 
          ENDIF   
        IF ( Range .gt. TwoPi ) THEN
            Range= DMOD( Range, TwoPi )
          ENDIF

        ! ----------------- Find Geocentric Latitude  -----------------
        TLat = asin( sin(LLat)*cos(Range) +
     &         cos(LLat)*sin(Range)*cos(Az) )

        ! ---- Find Delta N, the ANGLE between the points -------------
        IF ( (abs(cos(TLat)) .gt. Small) .and.
     &        (abs(cos(LLat)) .gt. Small) ) THEN
            SinDN = sin(Az)*sin(Range) / cos(TLat)
            CosDN = ( cos(Range)-sin(TLat)*sin(LLat) ) /
     &                 ( cos(TLat)*cos(LLat) )
            DeltaN= atan2(SinDN, CosDN)
          ELSE
            ! ------ Case where launch is within 3nm of a Pole --------
            IF ( abs(cos(LLat)) .le. Small ) THEN
                IF ( (Range .gt. Pi) .and. (Range .lt. TwoPi) ) THEN
                    DeltaN= Az + Pi
                  ELSE
                    DeltaN= Az 
                  ENDIF
              ENDIF
            ! ----- Case where ENDIF point is within 3nm of a pole ----
            IF ( abs( cos(TLat) ) .le. Small ) THEN
                DeltaN= 0.0D0 
              ENDIF
          ENDIF 

        TLon= LLon + DeltaN
        IF ( abs(TLon) .gt. TwoPi ) THEN
            TLon= DMOD( TLon, TwoPi )
          ENDIF
        IF ( TLon .lt. 0.0D0 ) THEN
            TLon= TwoPi + TLon 
          ENDIF   
      RETURN
      END
      
      end module ast2body
