      module astmath
      use, intrinsic:: ieee_arithmetic
      use, intrinsic:: iso_fortran_env, only: wp=>real64
      use comm
      implicit none
      
      contains

*    this file contains miscellaneous math functions.
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
*              28 Jan 04  David Vallado
*                           Update headers                                    
*              14 Feb 03  David Vallado
*                           Change vectors to (3) from (4)
*              14 May 01  David Vallado
*                           2nd edition baseline
*              23 Nov 87  David Vallado
*                           Original baseline
*

*                           FUNCTION FACTORIAL
*
*  this function finds the value of a FACTORIAL.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    N           - Input value
*
*  Outputs       :
*    FUNCTION    - answer
*
*  Locals        :
*    i           - Index
*
*  Coupling      :
*    None.
*
* ------------------------------------------------------------------------------

      elemental REAL(wp) FUNCTION  FACTORIAL(N)
       INTEGER, intent(in) :: N
       INTEGER i

       Factorial = product([(i, i=1, N)])


      END FUNCTION  FACTORIAL

* ------------------------------------------------------------------------------
*
*                           FUNCTION BINOMIAL
*
*  this function finds the value of a BINOMIAL coefficient
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    i           -
*    j           -
*
*  Outputs       :
*    FUNCTION    - answer
*
*  Locals        :
*    None.
*
*  Coupling      :
*    FACTORIAL   - Finds the FACTORIAL of a number
*
* ------------------------------------------------------------------------------

      elemental REAL(wp) FUNCTION  BINOMIAL( i,j )

        INTEGER, intent(in) :: i,j

       BINOMIAL= FACTORIAL(j) / ( FACTORIAL(i)*FACTORIAL(j-i) )

      END FUNCTION  BINOMIAL 

* ------------------------------------------------------------------------------
*
*                           SUBROUTINE PLANE
*
* this subroutine calculates the equation of a PLANE given 3 points
*   pt1 - x1,y1,z1, pt2 - x2,y2,z2, pt3 - x3,y3,z3 , and outputs the
*   a b c d variables describing the PLANE. NOTE that the general equation
*   of a PLANE is defined here as:  ax + by + cz + d = 0  and the values
*   are obtained by solving the ordered determinant  x  y  z   1
*                                                    x1 y1 z1  1   =0
*                                                    x2 y2 z2  1
*                                                    x3 y3 z3  1
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    x1,y1,z1    - point # 1
*    x2,y2,z2    - point # 2
*    x3,y3,z3    - point # 3
*
*  OutPuts       :
*    a,b,c,d     - constants for the equation of the PLANE
*
*  Locals        :
*    z23
*
*  Coupling      :
*    None.
*
* ------------------------------------------------------------------------------  

      SUBROUTINE PLANE       ( x1,y1,z1,x2,y2,z2,x3,y3,z3, a,b,c,d )
        IMPLICIT NONE
        REAL*8 x1,y1,z1,x2,y2,z2,x3,y3,z3,a,b,c,d
* -----------------------------  Locals  ------------------------------
        REAL*8 z23,y23,x23,yz23,yz32,xz23,xz32,xy23,xy32

        z23= z2-z3 
        y23= y2-y3 
        x23= x2-x3 

        yz23= y2*z3 
        yz32= y3*z2 
        xz23= x2*z3
        xz32= x3*z2 
        xy23= x2*y3 
        xy32= x3*y2 

        a= y1*z23 - z1*y23 + yz23 - yz32 
        b= -(x1*z23 - z1*x23 + xz23 - xz32)
        c= x1*y23 - y1*x23 + xy23 - xy32 
        d= -(x1*(yz23 - yz32) - y1*(xz23-xz32) + z1*(xy23 -xy32))
      RETURN
      END
*
* ------------------------------------------------------------------------------
*
*                           FUNCTION COT
*
*  this function finds the Cotangent of an ANGLE in radians.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    XVal        - ANGLE to take Cotangent of                   rad
*
*  OutPuts       :
*    COT         - Result
*
*  Locals        :
*    Temp        - Temporary Real variable
*
* ------------------------------------------------------------------------------

      Elemental REAL(wp) FUNCTION COT(XVal)
        IMPLICIT NONE
        REAL(wp), intent(in) :: XVal

        COT = 1/TAN(XVal)

      END FUNCTION COT




* ------------------------------------------------------------------------------
*
*                           SUBROUTINE CROSS
*
*  this subroutine crosses two vectors.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Vec1        - Vector number 1
*    Vec2        - Vector number 2
*
*  OutPuts       :
*    OutVec      - Vector result of A x B
*
*  Locals        :
*    None.
*
*  Coupling      :
*
* ------------------------------------------------------------------------------  

      SUBROUTINE CROSS       ( Vec1,Vec2, OutVec )
        IMPLICIT NONE
        REAL*8 Vec1(3), Vec2(3), OutVec(3)

        ! --------------------  Implementation   ----------------------
        OutVec(1)= Vec1(2)*Vec2(3) - Vec1(3)*Vec2(2)
        OutVec(2)= Vec1(3)*Vec2(1) - Vec1(1)*Vec2(3)
        OutVec(3)= Vec1(1)*Vec2(2) - Vec1(2)*Vec2(1)

      RETURN
      END


* ------------------------------------------------------------------------------
*
*                           SUBROUTINE NORM
*
*  this subroutine calculates a unit vector given the original vector.  If  a
*    zero vector is input, the vector is set to zero.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Vec         - Vector
*
*  OutPuts       :
*    OutVec      - Unit Vector
*
*  Locals        :
*    MagVec      - Magnitude of vector
*    i           - Index
*
*  Coupling      :
*    MAG           Magnitude of a vector
*
* ------------------------------------------------------------------------------  

      SUBROUTINE NORM        ( Vec, OutVec )
        IMPLICIT NONE
        REAL*8 Vec(3), OutVec(3)

* -----------------------------  Locals  ------------------------------
        REAL*8 MagVec
        INTEGER i

        ! --------------------  Implementation   ----------------------
        MagVec = norm2( Vec )
        IF ( MagVec .gt. Small ) THEN
            DO i= 1, 3
                OutVec(i)= Vec(i)/MagVec
              ENDDO
          ELSE
            DO i= 1, 3
                OutVec(i)= 0.0D0
              ENDDO
          ENDIF

      RETURN
      END

* ------------------------------------------------------------------------------
*
*                           SUBROUTINE ROTi
*
*  this subroutine performs a rotation about the ith axis. i is specified
*    for each operation.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Vec         - Input vector
*    XVal        - ANGLE of rotation              rad
*
*  OutPuts       :
*    OutVec      - Vector result
*
*  Locals        :
*    c           - Cosine of the ANGLE XVal
*    s           - Sine of the ANGLE XVal
*    Temp        - Temporary REAL*8 value
*
*  Coupling      :
*    None.
*
* ------------------------------------------------------------------------------

      SUBROUTINE ROT1        ( Vec, XVal, OutVec )
        IMPLICIT NONE
        REAL*8 Vec(3), XVal, OutVec(3)
* -----------------------------  Locals  ------------------------------
        REAL*8 c, s, Temp

        Temp= Vec(3)
        c= DCOS( XVal )
        s= DSIN( XVal )

        OutVec(3)= c*Vec(3) - s*Vec(2)
        OutVec(2)= c*Vec(2) + s*Temp
        OutVec(1)= Vec(1)
      RETURN
      END


      pure SUBROUTINE ROT2        ( Vec, XVal, OutVec )

        REAL(wp), intent(in) :: Vec(3), XVal
        real(wp), intent(out) :: OutVec(3)
* -----------------------------  Locals  ------------------------------
        REAL*8 c, s, Temp

        ! --------------------  Implementation   ----------------------
        Temp= Vec(3)
        c= COS( XVal )
        s= SIN( XVal )

        OutVec(3)= c*Vec(3) + s*Vec(1)
        OutVec(1)= c*Vec(1) - s*Temp
        OutVec(2)= Vec(2)

      END


      SUBROUTINE ROT3        ( Vec, XVal, OutVec )
        IMPLICIT NONE
        REAL*8 Vec(3), XVal, OutVec(3)
* -----------------------------  Locals  ------------------------------
        REAL*8 c, s, Temp

        ! --------------------  Implementation   ----------------------
        Temp= Vec(2)
        c= DCOS( XVal )
        s= DSIN( XVal )

        OutVec(2)= c*Vec(2) - s*Vec(1)
        OutVec(1)= c*Vec(1) + s*Temp
        OutVec(3)= Vec(3)
      RETURN
      END
*
* ------------------------------------------------------------------------------
*
*                                  rot1mat
*
*  this function sets up a rotation matrix for an input angle about the first
*    axis.
*
*  author        : david vallado                  719-573-2600   10 jan 2003
*
*  revisions
*                -
*
*  inputs          description                    range / units
*    xval        - angle of rotation              rad
*
*  outputs       :
*    outmat      - matrix result
*
*  locals        :
*    c           - cosine of the angle xval
*    s           - sine of the angle xval
*
*  coupling      :
*    none.
*
* ----------------------------------------------------------------------------- }

      SUBROUTINE ROT1MAT     ( XVal, OutMat )
        IMPLICIT NONE
        REAL*8 XVal, OutMat(3,3)

* -----------------------------  Locals  ------------------------------
        REAL*8 c, s

        c= Dcos( xval )
        s= Dsin( xval )

        outmat(1,1)= 1.0D0
        outmat(1,2)= 0.0D0
        outmat(1,3)= 0.0D0

        outmat(2,1)= 0.0D0
        outmat(2,2)= c
        outmat(2,3)= s

        outmat(3,1)= 0.0D0
        outmat(3,2)=  c
        outmat(3,3)= -s

      RETURN
      END

      SUBROUTINE ROT2MAT     ( XVal, OutMat )
        IMPLICIT NONE
        REAL*8 XVal, OutMat(3,3)

* -----------------------------  Locals  ------------------------------
        REAL*8 c, s

        c= Dcos( xval )
        s= Dsin( xval )

        outmat(1,1)= c
        outmat(1,2)= 0.0D0
        outmat(1,3)= -s

        outmat(2,1)= 0.0D0
        outmat(2,2)= 1.0D0
        outmat(2,3)= 0.0D0

        outmat(3,1)= c
        outmat(3,2)= 0.0D0
        outmat(3,3)= s

      RETURN
      END

      SUBROUTINE ROT3MAT     ( XVal, OutMat )
        IMPLICIT NONE
        REAL*8 XVal, OutMat(3,3)

* -----------------------------  Locals  ------------------------------
        REAL*8 c, s

        c= Dcos( xval )
        s= Dsin( xval )

        outmat(1,1)= c
        outmat(1,2)= s
        outmat(1,3)= 0.0D0

        outmat(2,1)= c
        outmat(2,2)= -s
        outmat(2,3)= 0.0D0

        outmat(3,1)= 0.0D0
        outmat(3,2)= 0.0D0
        outmat(3,3)= 1.0D0

      RETURN
      END



* ---- this subroutine takes an arbitrary number of points .and.
*  fits a polynomial of degree <degree>. The datapoints are stored in a
*  matrix having rows <NumPts> by 2 <one for x, and y>.
*  See Matthews pg 230-232.
*
*    Assumes the form of
*       ax2 + bx + c = y
*    coeff indices are
*       1 = a, 2 = b, 3 = c
*    MIN is at
*       2ax + b = 0
*    y is found by plugging a back into the polynomial
*
*  dav  3 Apr 97
* ---------------------------------------------------------------------

      SUBROUTINE POLYFIT     ( Degree,NumPts,DataPoints,Coeff,
     &                         MaxRow,MaxDeg,MinX,MinY )
        IMPLICIT NONE
        INTEGER Degree, NumPts,MaxRow,MaxDeg
        REAL*8 DataPoints(MaxRow,2),Coeff(MaxDeg,1), MinX, MinY
* -----------------------------  Locals  ------------------------------
        INTEGER j,k,r,c,MaxR
        PARAMETER (MaxR = 20)
        Real*8 p, ai(MaxR,MaxR), a(MaxR,MaxR),b(MaxR,1), parr(MaxR,1)

        ! --------------------  Implementation   ----------------------
        ! ----------- Find the sum of the product terms (x*y) ---------
        DO k= 1, NumPts
            p= 1
            DO r= 1, Degree+1
                b(r,1)=  b(r,1) + DataPoints(k,2)*p
                p= p*DataPoints(k,1)
              ENDDO
          ENDDO

        ! ------------ Find the sum of powers for x (x**) -------------
        DO k= 1, NumPts
            p= DataPoints(k,1)
            DO j= 1, 2*Degree
                  parr(j,1)= parr(j,1) + p
                p= p*DataPoints(k,1)
              ENDDO
          ENDDO

        ! ---------- Find the matrix entries for the equations --------
        DO r= 1, Degree+1
            DO c= 1, Degree+1
                IF ( r+c .ne. 2 ) THEN
                      a(r,c)= parr(r+c-2,1)
                  ELSE
                      a(r,c)= NumPts
                    ENDIF
              ENDDO
          ENDDO

        ! --------- Solve linear equations for coeffeicients ----------
        CALL MATINVERSE( A      , 3,Maxr         , AI )
        coeff = MATMUL( AI,b)

        ! ----------------- Find minimum values -----------------------
        MinX= -coeff(2,1) / ( 2*coeff(3,1) )
        MinY= coeff(3,1)*MinX*MinX + coeff(2,1)*MinX + coeff(1,1)

      RETURN
      END
*
* ------------------------------------------------------------------------------
*
*                           SUBROUTINE ANGLE
*
*  this subroutine calculates the ANGLE between two vectors.  The output is
*    set to 999999.1D0 to indicate an undefined value.  Be SURE to check  
*    this at the output phase.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Vec1        - Vector number 1
*    Vec2        - Vector number 2
*
*  OutPuts       :
*    Theta       - ANGLE between the two vectors  -Pi to Pi
*
*  Locals        :
*    Temp        - Temporary REAL variable
*
*  Coupling      :
*    DACOS        Arc Cosine FUNCTION
*
* ------------------------------------------------------------------------------  

      SUBROUTINE ANGLE( Vec1,Vec2, Theta )
        IMPLICIT NONE
        REAL*8 Vec1(3), Vec2(3), Theta, magvec1, magvec2
* -----------------------------  Locals  ------------------------------
        REAL*8 Temp
        

        ! --------------------  Implementation   ----------------------
        magvec1 = norm2(vec1)
        magvec2 = norm2(vec2)
        IF ( magVec1*magVec2 .gt. Small**2 ) THEN
            Temp= DOT_product(Vec1,Vec2) / (magVec1*magVec2)
            IF ( DABS( Temp ) .gt. 1.0D0 ) THEN
                Temp= DSIGN(1.0D0, Temp)
              ENDIF
            Theta= DACOS( Temp ) 
          ELSE
            Theta= Undefined
          ENDIF

      END subroutine angle


      SUBROUTINE DMulRSub    ( ALPR,ALPI,BETR,BETI )
        IMPLICIT NONE
        REAL*8 Alpr(4), Alpi(4), Betr(4), Beti(4)
* -----------------------------  Locals  ------------------------------
        REAL*8 Te1,Te2,Te3,Te4,Te5,Te6,Te7,Te8,Te9,Te10,Te11,Te12,
     &         Te13,Te14,Te15,Te16,TeM,DE15,DE16,TemTe7,TemTe8

        ! --------------------  Implementation   ----------------------
         Te1= ALPR(1)-ALPR(3)
         Te2= ALPI(1)-ALPI(3)
         Te5= ALPR(3)-ALPR(2)
         Te6= ALPI(3)-ALPI(2)
         TeM= Te5*Te5+Te6*Te6 

         ! --------------- Check for zero values of tem ---------------
         IF ( DABS( Tem ) .gt. 1.0D-20 ) THEN
             Te3= (Te1*Te5+Te2*Te6)/TeM
             Te4= (Te2*Te5-Te1*Te6)/TeM 
           ELSE
             Te3= 0.0D0
             Te4= 0.0D0 
           ENDIF 

         Te7 = Te3+1.0D0 
         Te9 = Te3*Te3-Te4*Te4 
         Te10= 2.0D0*Te3*Te4 
         DE15= Te7*BETR(3)-Te4*BETI(3)
         DE16= Te7*BETI(3)+Te4*BETR(3)
         Te11= Te3*BETR(2)-Te4*BETI(2)+BETR(1)-DE15
         Te12= Te3*BETI(2)+Te4*BETR(2)+BETI(1)-DE16

         Te7 = Te9-1.0D0 
         Te1 = Te9*BETR(2)-Te10*BETI(2)
         Te2 = Te9*BETI(2)+Te10*BETR(2)
         Te13= Te1-BETR(1)-Te7*BETR(3)+Te10*BETI(3)
         Te14= Te2-BETI(1)-Te7*BETI(3)-Te10*BETR(3)
         Te15= DE15*Te3-DE16*Te4 
         Te16= DE15*Te4+DE16*Te3 

         Te1= Te13*Te13-Te14*Te14-4.0D0*(Te11*Te15-Te12*Te16)
         Te2= 2.0D0*Te13*Te14-4.0D0*(Te12*Te15+Te11*Te16) 

* ---------------------------------------------------------------------
*   Sometimes, for stiff systems (the roots vary widely in order
*   of magnitude), Te1 and Te2 get large enough to have their
*   squares overflow the floating point range.  To prevent this,
*   when either one is large, they are scaled by 10**10 for the
*   purpose of finding TeM.  The scale is restored when the
*   magnitude computation is completed.  This should not affect
*   the accuracy of the computations, since the mantissa is not
*   affected, only the exponent.
* ---------------------------------------------------------------------

         IF ( ( Te1 .gt. 1.0D15 ) .or. ( Te2 .gt. 1.0D15 ) ) THEN
             Te1= Te1*1.0D-10
             Te2= Te2*1.0D-10
             TeM= 1.0D10*DSQRT(Te1*Te1+Te2*Te2)
           ELSE
             TeM= DSQRT(Te1*Te1+Te2*Te2)
           ENDIF

         IF ( Te1 .gt. 0.0D0 ) THEN
             Te3= DSQRT(0.5D0*(TeM+Te1))
             IF ( Te2 .lt. 0.0D0 ) THEN
                 Te3= -Te3
               ENDIF
             ! ---------------- check for zero values of te3 ----------
             IF ( DABS( Te3 ) .lt. 1.0D-15 ) THEN
                 Te4= 0.0D0
               ELSE
                 Te4= 0.5D0*Te2/Te3 
               ENDIF
           ELSE
             Te4= DSQRT(0.5D0*(TeM-Te1))
             ! -------------------- Check for underflows --------------
             IF ( DABS( Te4 ) .lt. 1.0D-15 ) THEN
                 Te3= 0.0D0
               ELSE
                 Te3= 0.5D0*Te2/Te4 
               ENDIF
           ENDIF 

         Te7 = Te13+Te3
         Te8 = Te14+Te4 
         Te9 = Te13-Te3
         Te10= Te14-Te4 
         Te1 = 2.0D0*Te15 
         Te2 = 2.0D0*Te16 

         IF ( (Te7*Te7+Te8*Te8-Te9*Te9-Te10*Te10) .le. 0.0D0 ) THEN
             Te7= Te9
             Te8= Te10 
           ENDIF 
         TeM= Te7*Te7+Te8*Te8 

         TemTe7= Tem*Te7 
         TemTe8= Tem*Te8 

         ! ------------ Check for values of almost zero ---------------
         IF ( DABS( TemTe7 ) .lt. 1.0D-20 ) THEN
             Te3= 0.0D0
             Te4= 0.0D0 
           ELSE
             Te3= Te1/TeM*Te7
             Te4= Te2/TeM*Te7 
           ENDIF 
         IF ( DABS( TemTe8 ) .gt. 1.0D-20 ) THEN
             Te3= Te3 + Te2/TeM*Te8
             Te4= Te4 - Te1/TeM*Te8 
           ENDIF

         ALPR(4)= ALPR(3)+Te3*Te5-Te4*Te6
         ALPI(4)= ALPI(3)+Te3*Te6+Te4*Te5
      RETURN
      END

* ------------------------------------------------------------------------------
*
*                           SUBROUTINE FACTOR
*
*  this subroutine is a root finding algorithm.  It takes the polynomial .and.
*    returns the roots (real and imaginary) in the Roots array.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Poly        - Array of 16 coefficients
*                    representing the polynomial
*                    (1) is x^8th, (2) is x^7th, ...
*                    others are zero
*    NRoots      - Number of roots
*
*  OutPuts       :
*    Roots       - Array containing roots (real,imag)
*
*  Locals        :
*                -
*                -
*                -
*
*  Coupling      :
*    DMulRSub    -
*
*  References    :
*    This is Bairstows method?
*
* ------------------------------------------------------------------------------  

      SUBROUTINE FACTOR      ( Poly, NRootS, RootS )
        IMPLICIT NONE
        REAL*8 Poly(16), Roots(15,2)
        INTEGER NRootS
* -----------------------------  Locals  ------------------------------
        REAL*8 Alpr(4), Alpi(4), Betr(4), Beti(4), Small, DPoly(16)
        REAL*8 Temp1,    Temp2,    AXR, AXi, PMax, Tem1, Tem2,
     &   TempReal, TempImag, Temp7
        LOGICAL Skip
        INTEGER Mode,LoopCnt,kk,i,j,l,RootCnt

        ! --------------------  Implementation   ----------------------
         Small = 0.000001D0
         PMax= 0.0D0 
         DO KK= 1, NRootS+1
             IF ( DABS(Poly(KK)) .gt. PMax ) THEN
                 PMax= Poly(KK)
               ENDIF
           ENDDO

         IF ( DABS(PMax) .lt. Small ) THEN
             PMax= 1.0D0
           ENDIF
      
         DO KK= 1, NRootS+1
             DPoly(kk)= Poly(kk)/PMax
           ENDDO

         IF ( NRootS .gt. 0 ) THEN
             RootCnt= 0
             i= NRootS+1 

             DO WHILE ((DABS(DPoly(i)) .lt. Small) .and.
     &                 (RootCnt .ne. NRoots))
                 RootCnt= RootCnt+1
                 RootS(RootCnt,1)= 0.0D0
                 RootS(RootCnt,2)= 0.0D0
                 i= i-1
               ENDDO

             IF ( RootCnt .ne. NRootS ) THEN
                 AXR= 0.8D0
                 AXi= 0.0D0
                 L  = 1 
                 LOOPCNT= 1 
                 AlpR(1)= AXR
                 Alpi(1)= AXi
                 MODE= 1 

                 DO WHILE (RootCnt .lt. NRoots)
                     BetR(4)= DPoly(1)
                     Beti(4)= 0.0D0
                     DO i= 1, NRootS
                         J= i+1
                         Temp1= BetR(4)*AXR - Beti(4)*AXi
                         Beti(4)= Beti(4)*AXR + BetR(4)*AXi
                         BetR(4)= Temp1 + DPoly(J)
                        ENDDO

                     TempReal= BetR(4)
                     Tempimag= Beti(4)

                     IF ( RootCnt .ne. 0 ) THEN
                         DO i= 1, RootCnt
                             Tem1   = AXR - RootS(i,1)
                             Tem2   = AXi - RootS(i,2)
                             Temp1  = Tem1*Tem1 + Tem2*Tem2 
                             Temp2  = (BetR(4)*Tem1+Beti(4)*Tem2)/Temp1
                             Beti(4)= (Beti(4)*Tem1-BetR(4)*Tem2)/Temp1
                             BetR(4)= Temp2
                           ENDDO
                       ENDIF 

                     IF (Mode .eq. 1 ) THEN
                             BetR(1)= BetR(4)
                             Beti(1)= Beti(4)
                             AXR    = 0.85D0 
                             AlpR(2)= AXR
                             Alpi(2)= AXi
                             MODE   = 2 
                           ENDIF
                     IF (Mode .eq. 2 ) THEN
                             BetR(2)= BetR(4)
                             Beti(2)= Beti(4)
                             AXR    = 0.9D0 
                             AlpR(3)= AXR
                             Alpi(3)= AXi
                             MODE   = 3
                           ENDIF
                     IF (Mode .eq. 3 ) THEN
                             BetR(3)= BetR(4)
                             Beti(3)= Beti(4)
                             CALL DMULRSUB( AlpR,Alpi,BetR,Beti )
                             AXR    = AlpR(4)
                             AXi    = Alpi(4)
                             MODE   = 4 
                           ENDIF
                     IF (Mode .eq. 5 ) THEN
                             BetR(1)=  BetR(4)
                             Beti(1)=  Beti(4)
                             AXR    =  AlpR(2)
                             AXi    = -Alpi(2)
                             Alpi(2)= -Alpi(2)
                             MODE   =  6 
                           ENDIF
                     IF (Mode .eq. 6 ) THEN
                             BetR(2)=  BetR(4)
                             Beti(2)=  Beti(4)
                             AXR    =  AlpR(3)
                             AXi    = -Alpi(3)
                             Alpi(3)= -Alpi(3)
                             L      =  2 
                             MODE   =  3 
                           ENDIF 

                     ! --------------- the convergence mode -----------
                     IF (Mode .eq. 4 ) THEN
                         Skip= .FALSE.
                         IF ( DABS(TempReal)+DABS(Tempimag) .gt.
     &                          1.0D-20 ) THEN
                             Temp7= DABS(AlpR(3)-AXR)+DABS(Alpi(3)-AXi)
                         IF ( Temp7/(DABS(AXR)+DABS(AXi)) .gt.
     &                          1.0D-7 ) THEN
                             LoopCnt = LoopCnt + 1
                             DO i= 1, 3
                                 AlpR(i)= AlpR(i+1)
                                 Alpi(i)= Alpi(i+1)
                                 BetR(i)= BetR(i+1)
                                 Beti(i)= Beti(i+1)
                               ENDDO
                             IF ( LOOPCNT .lt. 100 ) THEN
                                 CALL DMULRSUB( AlpR,Alpi,BetR,Beti )
                                 AXR = AlpR(4)
                                 AXi = Alpi(4)
                                 MODE= 4
                                 Skip= .TRUE.
                               ENDIF
                           ENDIF

                         IF ( Skip .eqv. .FALSE. ) THEN
                             RootCnt= RootCnt+1
                             RootS(RootCnt,1)= AlpR(4)
                             RootS(RootCnt,2)= Alpi(4)
                             LOOPCNT= 0

                             IF ( RootCnt .lt. NRootS ) THEN
                                 IF ( DABS(RootS(RootCnt,2)) .gt.
     &                                 Small ) THEN
                                     IF ( L .eq. 1 ) THEN
                                         AXR    =  AlpR(1)
                                         AXi    = -Alpi(1)
                                         Alpi(1)= -Alpi(1)
                                         MODE   =  5
                                       ELSE
                                         AXR= 0.8D0
                                         AXi= 0.0D0
                                         L  = 1
                                         LOOPCNT= 1
                                         AlpR(1)= AXR
                                         Alpi(1)= AXi
                                         MODE   = 1
                                       ENDIF
                                   ELSE
                                     AXR= 0.8D0
                                     AXi= 0.0D0
                                     L  = 1
                                     LOOPCNT= 1
                                     AlpR(1)= AXR
                                     Alpi(1)= AXi
                                     MODE   = 1
                                   ENDIF
                               ENDIF
                           ENDIF   ! IF ( Skip
                           ENDIF 
                     ENDIF   ! Case

                   ENDDO   ! DO WHILE

             ENDIF   ! If  RootCnt .ne. NRoots

           ENDIF   ! If  nRoots .gt. 0

      RETURN
      END

* ------------------------------------------------------------------------------
*
*                           SUBROUTINE QUADRATIC
*
*  this subroutine solves for the two Roots of a QUADRATIC equation.  There are
*    no restrictions on the coefficients, and imaginary results are passed
*    out as separate values.  The general form is y = ax2 + bx + c.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    a           - Coefficient of x squared term
*    b           - Coefficient of x term
*    c           - Constant
*
*  OutPuts       :
*    R1r         - Real portion of Root 1
*    R1i         - Imaginary portion of Root 1
*    R2r         - Real portion of Root 2
*    R2i         - Imaginary portion of Root 2
*
*  Locals        :
*    Discrim     - Discriminate b2 - 4ac
*
*  Coupling      :
*    None.
*
*  References    :
*    Vallado       2007, pg 974
*
* ------------------------------------------------------------------------------  

      SUBROUTINE QUADRATIC   ( a,b,c, R1r,R1i,R2r,R2i )
        IMPLICIT NONE
        REAL*8 R1r,R1i,R2r,R2i,a,b,c
* -----------------------------  Locals  ------------------------------
        REAL*8 Discrim

        ! --------------------  Implementation   ----------------------
        R1r= 0.0D0 
        R1i= 0.0D0 
        R2r= 0.0D0 
        R2i= 0.0D0 

        Discrim= b*b - 4.0D0*a*c 

        ! ---------------------  Real Roots  --------------------------
        IF ( Discrim .gt. 0.0D0 ) THEN
            R1r= ( -b + DSQRT(Discrim) ) / ( 2.0D0*a )
            R2r= ( -b - DSQRT(Discrim) ) / ( 2.0D0*a ) 
          ELSE
            ! -------------------- Complex Roots ----------------------
            R1r= -b / ( 2.0D0*a )
            R2r= R1r
            R1i=  DSQRT(-Discrim) / ( 2.0D0*a )
            R2i= -DSQRT(-Discrim) / ( 2.0D0*a )
          ENDIF
      RETURN
      END
*
* ------------------------------------------------------------------------------
*
*                           SUBROUTINE CUBIC
*
*  this subroutine solves for the three Roots of a CUBIC equation.  There are
*    no restrictions on the coefficients, and imaginary results are passed
*    out as separate values.  The general form is y = ax3 + bx2 + cx + d.  Note
*    that R1i will ALWAYS be ZERO since there is ALWAYS at least one REAL Root.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    a           - Coefficient of x cubed term
*    b           - Coefficient of x squared term
*    c           - Coefficient of x term
*    d           - Constant
*
*  OutPuts       :
*    R1r         - Real portion of Root 1
*    R1i         - Imaginary portion of Root 1
*    R2r         - Real portion of Root 2
*    R2i         - Imaginary portion of Root 2
*    R3r         - Real portion of Root 3
*    R3i         - Imaginary portion of Root 3
*
*  Locals        :
*    Temp1       - Temporary value
*    Temp2       - Temporary value
*    Root1       - Temporary value of the Root
*    Root2       - Temporary value of the Root
*    Root3       - Temporary value of the Root
*    P           - Coefficient of x squared term where x cubed term is 1.0D0
*    Q           - Coefficient of x term where x cubed term is 1.0D0
*    R           - Coefficient of constant term where x cubed term is 1.0D0
*    Delta       - Discriminator for use with Cardans formula
*    E0          - ANGLE holder for trigonometric solution
*    Phi         - ANGLE used in trigonometric solution
*    CosPhi      - Cosine of Phi
*    SinPhi      - Sine of Phi
*
*  Coupling      :
*    None.
*
*  References    :
*    Vallado       2007, pg 975
*
* ------------------------------------------------------------------------------
*
      SUBROUTINE CUBIC       ( a,b,c,d, R1r,R1i,R2r,R2i,R3r,R3i )
        IMPLICIT NONE
        REAL*8 R1r,R1i,R2r,R2i,R3r,R3i,a,b,c,d
* -----------------------------  Locals  ------------------------------
        REAL*8 OneThird,
     &   Temp1, Temp2, Root1, Root2, Root3, P, Q, R, Delta,
     &   E0, CosPhi, SinPhi, Phi
        

        ! --------------------  Implementation   ----------------------
        OneThird  = 1.0D0/3.0D0
        R1r  = 0.0D0
        R1i  = 0.0D0
        R2r  = 0.0D0
        R2i  = 0.0D0
        R3r  = 0.0D0
        R3i  = 0.0D0
        Root1= 0.0D0 
        Root2= 0.0D0 
        Root3= 0.0D0 

        ! ----------- Force coefficients into std form ----------------
        P= B/A 
        Q= C/A 
        R= D/A 

        a= OneThird*( 3.0D0*Q - P*P ) 
        b= (1.0D0/27.0D0)*( 2.0D0*P*P*P - 9.0D0*P*Q + 27.0D0*R )

        Delta= (a*a*a/27.0D0) + (b*b*0.25D0) 

        ! ------------------ Use Cardans formula ----------------------
        IF ( Delta .gt. Small ) THEN
            Temp1= (-b*0.5D0)+DSQRT(Delta)
            Temp2= (-b*0.5D0)-DSQRT(Delta) 
            IF ( DABS(Temp1) .gt. Small ) THEN
                Temp1= Temp1**OneThird
              ENDIF
            IF ( DABS(Temp2) .gt. Small ) THEN
                Temp2= Temp2**OneThird
              ENDIF
            Root1= Temp1 + Temp2 
            Root2= -0.5D0*(Temp1 + Temp2) 
            Root3= -0.5D0*(Temp1 + Temp2) 
            R2i= -0.5D0*DSQRT( 3.0D0 )*(Temp1 - Temp2) 
            R3i= -R2i 
          ELSE
            ! --------------- Evaluate zero point ---------------------
            IF ( DABS( Delta ) .lt. Small  ) THEN
                IF ( DABS(b) .gt. Small ) THEN
                    Root1= -2.0D0 * (b*0.5D0)**OneThird
                    Root2=  (b*0.5D0)**OneThird
                    Root3= Root2 
                  ENDIF 
              ELSE
                ! ------------ Use trigonometric identities -----------
                E0     = 2.0D0*DSQRT(-a*OneThird)
                CosPhi = (-b/(2.0D0*DSQRT(-a*a*a/27.0D0)) ) 
                SinPhi = DSQRT( 1.0D0-CosPhi*CosPhi ) 
                Phi    = DATAN2( SinPhi,CosPhi ) 
                Root1= E0*DCOS( Phi*OneThird ) 
                Root2= E0*DCOS( Phi*OneThird + 120.0D0*deg2Rad )
                Root3= E0*DCOS( Phi*OneThird + 240.0D0*deg2Rad )
              ENDIF 
          ENDIF 

        R1r= Root1 - P*OneThird 
        R2r= Root2 - P*OneThird 
        R3r= Root3 - P*OneThird 
      RETURN
      END
*
* ------------------------------------------------------------------------------
*
*                           SUBROUTINE QUARTIC
*
*  this subroutine solves for the four Roots of a QUARTIC equation.  There are
*    no restrictions on the coefficients, and imaginary results are passed
*    out as separate values.  The general form is y = ax4 + bx3 + cx2 + dx + e.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    a           - Coeficient of x fourth term
*    b           - Coefficient of x cubed term
*    c           - Coefficient of x squared term
*    d           - Coefficient of x term
*    e           - Constant
*
*  OutPuts       :
*    R1r         - Real portion of Root 1
*    R1i         - Imaginary portion of Root 1
*    R2r         - Real portion of Root 2
*    R2i         - Imaginary portion of Root 2
*    R3r         - Real portion of Root 3
*    R3i         - Imaginary portion of Root 3
*    R4r         - Real portion of Root 4
*    R4i         - Imaginary portion of Root 4
*
*  Locals        :
*    Temp1       - Temporary value
*    Temp2       - Temporary value
*    Root1       - Temporary value of the Root
*    Root2       - Temporary value of the Root
*    Root3       - Temporary value of the Root
*    s           - alternate variable
*    h           - Temporary value
*    hSqr        - h squared
*    hCube       - h Cubed
*    P           - Term in auxillary equation
*    Q           - Term in auxillary equation
*    R           - Term in auxillary equation
*    Delta       - Discriminator for use with Cardans formula
*    E0          - ANGLE holder for trigonometric solution
*    Phi         - ANGLE used in trigonometric solution
*    CosPhi      - Cosine of Phi
*    SinPhi      - Sine of Phi
*    RPrime      - Values of Roots before final work
*    Temp        - Temporary variable in finding MAX RPrime
*    Eta         - Constant coefficient in QUADRATIC solutions
*    Beta        - Constant coefficient in QUADRATIC solutions
*
*  Coupling      :
*    QUADRATIC     Find roots of a QUADRATIC polynomial
*
*  References    :
*    Vallado       2007, pg 976
*
* ------------------------------------------------------------------------------  

      SUBROUTINE QUARTIC     ( a,b,c,d,e, R1r,R1i,R2r,R2i,
     &                         R3r,R3i,R4r,R4i )
        IMPLICIT NONE
        REAL*8 R1r,R1i,R2r,R2i,R3r,R3i,R4r,R4i,a,b,c,d,e
* -----------------------------  Locals  ------------------------------
        REAL*8 OneThird,
     &   Temp1, Temp2, Root1, Root2, Root3, s, h, P, Q, R, Delta, E0,
     &   CosPhi, SinPhi, Phi, RPrime, hSqr, HCube, Eta, Beta
        

        ! --------------------  Implementation   ----------------------
        OneThird  = 1.0D0/3.0D0
        R1r  = 0.0D0 
        R1i  = 0.0D0 
        R2r  = 0.0D0 
        R2i  = 0.0D0 
        R3r  = 0.0D0
        R3i  = 0.0D0 
        R4r  = 0.0D0 
        R4i  = 0.0D0
        Root1= 0.0D0 
        Root2= 0.0D0 
        Root3= 0.0D0 
        ! ----------- Force coefficients into std form ----------------
        b= B/A 
        c= C/A 
        d= D/A 
        e= E/A 

        H= -b/4 
        HSqr= H**2
        HCube= HSqr * H 

        P=                      6.0D0*HSqr   + 3.0D0*b*h + c 
        Q=          4.0D0*HCube + 3.0D0*b*HSqr + 2.0D0*c*h + d 
        R= h*HCube +  b*HCube +   c*HSqr   +    d*h  + e 

        a=  OneThird*( -P*P-12.0D0*R )
        b=  (1.0D0/27.0D0)*( -2.0D0*P*P*P+72.0D0*P*R-27.0D0*Q*Q )
        s= -2.0D0*OneThird*P 

        Delta= (a*a*a/27.0D0) + (b*b*0.25D0) 

        IF ( DABS(Q) .gt. Small ) THEN
            ! ------------------ Use Cardans formula ------------------
            IF ( Delta .gt. Small ) THEN
                Temp1= (-b*0.5D0)+DSQRT(Delta)
                Temp2= (-b*0.5D0)-DSQRT(Delta) 
                IF ( DABS(Temp1) .gt. Small ) THEN
                    Temp1= Temp1**OneThird
                  ENDIF
                IF ( DABS(Temp2) .gt. Small ) THEN
                    Temp2= Temp2**OneThird
                  ENDIF   
                Root1= Temp1 + Temp2 
                Root2= -0.5D0*(Temp1 + Temp2) 
                Root3= -0.5D0*(Temp1 + Temp2) 
                R2i= -0.5D0*DSQRT( 3.0D0 )*(Temp1 - Temp2) 
                R3i= -R2i 
              ELSE
                ! --------------- Evaluate zero point -----------------
                IF ( DABS( Delta ) .lt. Small  ) THEN
                    IF ( DABS(b) .gt. Small ) THEN
                        Root1= -2.0D0* (b*0.25D0)**OneThird
                        Root2=  (b*0.25D0)**OneThird
                        Root3= Root2 
                      ENDIF 
                  ELSE
                    ! ------ Use trigonometric identities -------------
                    E0    = 2.0D0*DSQRT(-a*OneThird)
                    CosPhi= (-b/(2.0D0*DSQRT(-a*a*a/27.0D0)) ) 
                    SinPhi= DSQRT( 1.0D0-CosPhi*CosPhi ) 
                    Phi   = DATAN2( SinPhi,CosPhi ) 
                    Root1= E0*DCOS( Phi*OneThird ) 
                    Root2= E0*DCOS( Phi*OneThird + 120.0D0*deg2Rad )
                    Root3= E0*DCOS( Phi*OneThird + 240.0D0*deg2Rad )
                  ENDIF 
              ENDIF 

             ! --------------- Find largest value of Root -------------
             RPrime= Root1+s 
             IF ( (RPrime .lt. Root2+s) .and.
     &            (DABS(R2i).lt.0.0001D0) ) THEN
                 RPrime= Root2+s
              ENDIF
             IF ( (RPrime .lt. Root3+s) .and.
     &            (DABS(R3i).lt.0.0001D0) ) THEN
                 RPrime= Root3+s 
              ENDIF

             ! -- Evaluate coefficients of two resulting Quadratics ---
             IF ( RPrime .gt. Small ) THEN
                 Eta = 0.5D0*( P + RPrime - Q/DSQRT(RPrime) )
                 Beta= 0.5D0*( P + RPrime + Q/DSQRT(RPrime) ) 
               ELSE
                 Eta = 0.5D0*P
                 Beta= 0.5D0*P 
                 IF ( RPrime .lt. 0.0D0 ) THEN
                     RPrime= -RPrime 
                   ENDIF
               ENDIF 

             CALL QUADRATIC( 1.0D0, DSQRT(RPrime),Eta,
     &                       R1r,R1i,R2r,R2i )
             CALL QUADRATIC( 1.0D0,-DSQRT(RPrime),Beta,
     &                       R3r,R3i,R4r,R4i )

           ELSE
             ! ------ Case where solution reduces to a QUADRATIC -------
             CALL QUADRATIC( 1.0D0,P,R,   R1r,R1i,R3r,R3i )
             R  = DSQRT( R1r*R1r + R1i*R1i ) 
             IF ( DABS(R) .gt. Small ) THEN
                 Phi= DATAN2( R1i/R,R1r/R )
               ELSE
                 Phi= 0.0D0
               ENDIF
             R1r= DSQRT(R) * DCOS(Phi*0.5D0)
             R1i= DSQRT(R) * DSIN(Phi*0.5D0) 
             IF ( DABS(R1i) .gt. Small ) THEN
                 R2r= R1r
               ELSE
                 R2r= -R1r
               ENDIF
             R2i= -R1i 

             R  = DSQRT( R3r*R3r + R3i*R3i ) 
             IF ( DABS(R) .gt. Small ) THEN
                 Phi= DATAN2( R3i/R,R3r/R )
               ELSE
                 Phi= 0.0D0
               ENDIF
             R3r= DSQRT(R) * DCOS(Phi*0.5D0) 
             R3i= DSQRT(R) * DSIN(Phi*0.5D0) 
             IF ( DABS(R3i) .gt. Small ) THEN
                 R4r= R3r
               ELSE
                 R4r= -R3r
               ENDIF
             R4i= -R3i 
           ENDIF 

        R1r= R1r + h 
        R2r= R2r + h
        R3r= R3r + h 
        R4r= R4r + h 
      RETURN
      END




* ------------------------------------------------------------------------------
*
*                           SUBROUTINE MAKEMAT
*
*  this subroutine forms a rotation matrix for a given axis of rotation.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Angl        - ANGLE through which to rotate
*    Numbr       - Number of axis for rotation
*
*  OutPuts       :
*    Matr        - Matrix of rotation
*
*  Locals        :
*    None.
*
*  Coupling      :
*
* ------------------------------------------------------------------------------  

      SUBROUTINE MAKEMAT     ( Angl, Numbr, Matr )
        IMPLICIT NONE
        REAL*8 Angl, Matr(3,3)
        INTEGER Numbr

        ! --------------------  Implementation   ----------------------
        IF (Numbr .eq.1 ) THEN
            Matr(1,1)=  1.0D0
            Matr(2,2)=  DCOS(Angl)
            Matr(2,3)=  DSIN(Angl)
            Matr(3,2)= -DSIN(Angl)
            Matr(3,3)=  DCOS(Angl)
          ENDIF
        IF (Numbr .eq.2 ) THEN
            Matr(2,2)=  1.0D0
            Matr(1,1)=  DCOS(Angl)
            Matr(1,3)= -DSIN(Angl)
            Matr(3,1)=  DSIN(Angl)
            Matr(3,3)=  DCOS(Angl)
          ENDIF
        IF (Numbr .eq.3 ) THEN
            Matr(1,1)=  DCOS(Angl)
            Matr(1,2)=  DSIN(Angl)
            Matr(2,1)= -DSIN(Angl)
            Matr(2,2)=  DCOS(Angl)
            Matr(3,3)=  1.0D0
          ENDIF
      RETURN
      END
*
* ------------------------------------------------------------------------------
*
*                           SUBROUTINE LUDECOMP
*
*  this subroutine decomposes a matrix into an LU form.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*                  Dave Vallado                     719-554-3638    1 Aug 1989
*
*  Inputs          Description                    Range / Units
*    Order       - Order of matrix
*
*  OutPuts       :
*    LU          - LU decomposition matrix
*    Index       - Index vector for pivoting
*
*  Locals        :
*    i           - Index
*    j           - Index
*    k           - Index
*    IMax        - Pivot row pointer
*    Scale       - Scale FACTOR vector
*    Sum         - Temporary Variables
*    AMax        - Temporary Variables
*    Dum         - Temporary Variables
*
*  Coupling      :
*
*  References    :
*    Numerical Recipes - Flannery
*
* ------------------------------------------------------------------------------
*
      SUBROUTINE LUDeComp    ( LU, Index, Order, MaxRow )
        IMPLICIT NONE
        INTEGER Order, Index(Order), MaxRow
        REAL*8 LU(MaxRow,MaxRow)
* -----------------------------  Locals  ------------------------------
        REAL*8 Small, Scale(30), Sum, AMax, Dum ! Hardwire Scale to 30
        INTEGER I, J, K, IMax

        ! --------------------  Implementation   ----------------------
        Small= 0.000001D0
        IMax = 0 
        DO I = 1, Order
            AMax = 0.0D0
            DO J = 1, Order
                IF ( (DABS(  LU(i,j)) .gt. AMax) ) THEN
                    AMax = DABS( LU(i,j))
                  ENDIF
              ENDDO
            IF ( DABS(AMax) .lt. Small ) THEN
                Write(*,*) ' Singular Matrix AMax '
              ENDIF 
            Scale(i)= 1.0D0 / AMax
          ENDDO

        DO j = 1, Order
            DO i = 1, j - 1
                Sum =   LU(i,j)
                DO k = 1, i - 1
                    Sum = Sum -   LU(i,k)*  LU(k,j)
                  ENDDO
                LU(i,j)= Sum
              ENDDO
            AMax = 0.0D0 
            DO i = j, Order
                Sum =   LU(i,j)
                DO k = 1, j - 1
                    Sum = Sum -   LU(i,k)*  LU(k,j)
                  ENDDO
                LU(i,j)=  Sum
                Dum =   Scale(i )*DABS(Sum)
                IF ( (Dum .ge. AMax) ) THEN
                    IMax = i
                    AMax = Dum 
                  ENDIF 
              ENDDO
            IF ( (j .ne. iMax) ) THEN
                DO k = 1, Order
                    Dum =   LU(iMax,K)
                    LU(iMax,k)=   LU(j,k )
                    LU(j,k)=   Dum
                  ENDDO
                Scale(iMax)=   Scale(j)
              ENDIF 
            Index(j) = IMax
            IF ( DABS(  LU(j,j)) .lt. Small ) THEN
                Write(*,*) ' Matrix is Singular LU '
              ENDIF 
            IF ( (j .ne. Order) ) THEN
                Dum = 1.0D0 /   LU(j,j)
                DO i = j + 1, Order
                    LU(i,j)=  Dum*  LU(i,j)
                  ENDDO
              ENDIF 
          ENDDO
      RETURN
      END
*
* ------------------------------------------------------------------------------
*
*                           SUBROUTINE LUBkSUB
*
*  this subroutine finds the inverse of a matrix using LU decomposition.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Order       - Order of matrix
*    LU          - LU decomposition matrix
*    Index       - Index vector for pivoting
*
*  OutPuts       :
*    B           - Solution Vector
*
*  Locals        :
*    i           - Index
*    j           - Index
*    I0          - Pointer to non-zero element
*    IPtr        - Pivot Rwo Pointer
*    Sum         - Temporary Variables
*
*  Coupling      :
*
*  References    :
*    Numerical Recipes - Flannery
*
* ------------------------------------------------------------------------------  

      SUBROUTINE LUBkSub     ( LU, Index, Order, B, MaxRow )
        IMPLICIT NONE
        INTEGER Order, MaxRow, Index(Order)
        REAL*8 LU(MaxRow,MaxRow), B(MaxRow)
* -----------------------------  Locals  ------------------------------
        INTEGER I, J, IPtr, I0
        REAL*8 Sum

        ! --------------------  Implementation   ----------------------
        I0 = 0 
        DO i = 1, Order
            IPtr = Index(i)
            Sum =   B(IPtr)
            B(Iptr)=   B(i)
            IF ( (I0 .ne. 0) ) THEN
                DO j = I0, i - 1
                    Sum = Sum -   LU(i,j)*  B(j)
                  ENDDO
              ELSE
                IF ( (Sum .ne. 0.0D0) ) THEN
                   I0 = I
                 ENDIF
              ENDIF
            B(i)=  Sum
          ENDDO

        B(Order)=   B(Order)/  LU(Order,Order)

        DO i = (Order - 1), 1, -1
            Sum =   B(i)
            DO j = i + 1, Order
                Sum = Sum -   LU(i,j)*  B(j)
              ENDDO
            B(i)= Sum /   LU(i,i)
          ENDDO
      RETURN
      END
*
* ------------------------------------------------------------------------------
*
*                           SUBROUTINE MATINVERSE
*
*  this subroutine finds the inverse of a matrix using LU decomposition.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Mat         - Matrix to invert
*    Order       - Order of matrix
*
*  OutPuts       :
*    MatInv      - Inverted matrix
*
*  Locals        :
*    i           - Index
*    j           - Index
*    Index       - Index vector for pivoting
*    LU          - LU decomposition matrix
*    B           - Operational vector to form MatInv
*
*  Coupling      :
*    LUDeComp    -
*    LUBkSub     -
*
*  References    :
*    Numerical Recipes - Flannery
*
* ------------------------------------------------------------------------------  

      SUBROUTINE MATINVERSE  ( Mat, Order,MaxRow, MatInv )
        IMPLICIT NONE
        INTEGER Order,MaxRow
        REAL*8 Mat(MaxRow,MaxRow), MatInv(MaxRow,MaxRow)
* -----------------------------  Locals  ------------------------------
        INTEGER MaxR
        PARAMETER (MaxR = 10)
        INTEGER I, J, Index( MaxR )
        REAL*8 LU(MaxR,MaxR),B(MaxR)

        ! --------------------  Implementation   ----------------------
         DO i = 1, Order
             Index(i) = i
             DO j = 1, Order
                 LU(i,j)=   Mat(i,j)
               ENDDO
           ENDDO 

         CALL LUDeComp(LU, Index, Order, MaxR )

         DO j = 1, Order
             DO i = 1, Order
                IF ( (i .eq. j) ) THEN
                    B(i)= 1.0D0
                  ELSE
                    B(i)= 0.0D0
                  ENDIF

               ENDDO
             CALL LUBkSub(LU, Index, Order, B, MaxR)

             DO i = 1, Order
                 MatInv(i,j)=    B(i)
               ENDDO
           ENDDO

      RETURN
      END
*
* ------------------------------------------------------------------------------
*
*                           SUBROUTINE WRITEMAT
*
*  this subroutine WRITEs a matrix.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Mat1        - Matrix to WRITE out
*    DecNum      - Number of decimals to write
*
*  OutPuts       :
*    None.
*
*  Locals        :
*    Row         - Row Index
*    Col         - Column Index
*
*  Coupling      :
*    GETVAL      - Gets a value from the matrix
*
* ------------------------------------------------------------------------------

      SUBROUTINE WRITEMAT    ( Mat1, Mat1r,Mat1c,Max1r,Max1c,Title )
        IMPLICIT NONE
        INTEGER Mat1r,Mat1c,Max1r,Max1c
        REAL*8 Mat1( Max1r,Max1c )
        CHARACTER*64 Title
* -----------------------------  Locals  ------------------------------
        INTEGER  Row,Col

        ! --------------------  Implementation   ----------------------
        Write(*,'(A)') Title
        DO Row= 1, Mat1r
            Write(*,20) (Mat1(Row,Col),Col= 1, Mat1c)
          ENDDO
   20   FORMAT( 20(F12.7,1X) )
      RETURN
      END

* ------------------------------------------------------------------------------
*
*                           SUBROUTINE READMAT
*
*  this subroutine reads a matrix.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Mat1        - Matrix to WRITE out
*    Mat1r       - Matrix number 1 rows
*    Mat1c       - Matrix number 1 columns
*
*  OutPuts       :
*    None.
*
*  Locals        :
*    Row         - Row Index
*    Col         - Column Index
*
*  Coupling      :
*    None.
*
* ------------------------------------------------------------------------------  

      SUBROUTINE READMAT    ( Mat1, Mat1r,Mat1c,Max1r,Max1c )
        IMPLICIT NONE
        INTEGER Mat1r,Mat1c,Max1r,Max1c
        REAL*8 Mat1( Max1r,Max1c )
* -----------------------------  Locals  ------------------------------
        INTEGER  Row,Col

        ! --------------------  Implementation   ----------------------
        DO Row= 1, Mat1r
            Read(10,*) (Mat1(Row,Col),Col= 1, Mat1c)
          ENDDO
      RETURN
      END
*
      SUBROUTINE FILEWRITEMAT ( Mat1, Mat1r,Mat1c,Max1r,Max1c,Title )
        IMPLICIT NONE
        INTEGER Mat1r,Mat1c,Max1r,Max1c
        REAL*8 Mat1( Max1r,Max1c )
        CHARACTER*64 Title
* -----------------------------  Locals  ------------------------------
        INTEGER row,col

        ! --------------------  Implementation   ----------------------
        Write(20,'(A)') Title
        DO Row= 1, Mat1r
            Write(20,*) (Mat1(Row,Col),Col= 1, Mat1c)
          ENDDO
      RETURN
      END


* ----- WRITE a matrix using exponential formats -------  

      SUBROUTINE FILEEXPWRITEMAT ( Mat1,Mat1r,Mat1c,Max1r,Max1c,Title )
        IMPLICIT NONE
        INTEGER Mat1r, Mat1c,Max1r,Max1c
        REAL*8 Mat1(Max1r,Max1c)
        CHARACTER*64 Title
* -----------------------------  Locals  ------------------------------
        INTEGER  row,col

        ! --------------------  Implementation   ----------------------
        Write(20,'(A)') Title
        DO Row= 1, Mat1r
            Write(20,*) (' ',Mat1(Row,Col),Col= 1,Mat1c)
          ENDDO
      RETURN
      END
*
* ------------------------------------------------------------------------------
*
*                           FUNCTION DETERMINANT
*
*  this function calculates the determinant value using L-U decompisition.
*  The formula must have a NON-ZERO number in the 1,1 position.  If  the
*  FUNCTION senses a NON-ZERO number in row 1, it exchanges row1 for a row
*  WITH a NON-ZERO number.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*                                                                16 Jun 1993
*  Inputs          Description                    Range / Units
*    Order       - Order of determinaent (# of rows)
*    Mat1        - Matrix to find determinant of
*
*  OutPuts       :
*    Determinant - Result
*
*  Locals        :
*    i           - Index
*    j           - Index
*    k           - Index
*    n           - Index
*    Temp        -
*    D           -
*    Sum         -
*    L           -
*    U           -
*    Small       - Tolarance for comparing to 0.0D0
*
*  Coupling      :
*
*  References    :
*    Marion        pg. 168-172, 126-127
*
* ------------------------------------------------------------------------------  

      REAL(wp) FUNCTION DETERMINANT ( Mat1, Order )

        INTEGER, intent(in) :: Order
        REAL(wp), intent(inout) :: Mat1(:,:)
* -----------------------------  Locals  ------------------------------
        REAL(wp) Small, Temp, D, Sum
        INTEGER i,j,k,MaxR
        PARAMETER (MaxR = 10)
        REAL(wp) L(MaxR,MaxR), U(MaxR,MaxR)

        ! --------------------  Implementation   ----------------------
        Small = 0.00000001D0
        Sum   = 0.0D0
        ! ---------- Switch a non zero row to the first row -----------
        IF ( ABS(   Mat1(1,1 ) ) .lt. Small ) THEN
            j= 1
            DO WHILE (j .le. Order)
                IF ( ABS(   Mat1(j,1 ) ) .gt. Small ) THEN
                    DO k= 1, Order
                        Temp=    Mat1(1,k )
                          Mat1(1,k)=   Mat1(j,k )
                          Mat1(j,k)= Temp
                      ENDDO
                    j= Order + 1 
                  ENDIF 
                j= j+1 
              ENDDO
          ENDIF  ! If  DABS(Mat1(1,1)) .lt. Small

        DO i= 1, Order
              L(i,1)=  Mat1(i,1 )
          ENDDO
        DO j= 1, Order
              U(1,j)=  Mat1(1,j ) /   L(1,1 )
          ENDDO
        DO j= 2, Order
            DO i= j, Order
                Sum= 0.0D0
                DO k= 1, j-1
                    Sum= Sum+   L(i,k )*   U(k,j )
                  ENDDO
                L(i,j)=  Mat1(i,j )- Sum
              ENDDO   ! DO i
              U(j,j)=  1.0D0
            IF ( j .ne. Order ) THEN
                DO i= j+1, Order
                    Sum= 0.0D0
                    DO k= 1, j-1
                        Sum= Sum +    L(j,k)*   U(k,i )
                      ENDDO
                      U(j,i)=  (Mat1(j,i)-Sum)/   L(j,j )
                  ENDDO   ! DO i
              ENDIF   ! If  j
          ENDDO   ! DO j
        D= 1.0D0 
        DO i= 1, Order
            D= D*   L(i,i)
          ENDDO
        Determinant= D 


      END
      
      end module astmath
