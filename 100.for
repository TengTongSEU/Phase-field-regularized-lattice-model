C     *************************************************************************************************
C
C     Common variable module, which defines some common variables 
C     and allocates memory
C
C     *************************************************************************************************
C      
      module vars_module
          parameter (NodeNum = 111, NumEle= 10000, NPT=2)
          real*8,save :: allD(NumEle),allH(NumEle)
          real*8,save :: allElem(NumEle)
      end module

      
C     *************************************************************************************************
C
C     VUEL for the displacement field  
C     Activated DOFs: 1 and 2
C
C     *************************************************************************************************
C         
      subroutine vumat(
C Read only -
     1  nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     2  stepTime, totalTime, dt, cmname, coordMp, charLength,
     3  props, density, strainInc, relSpinInc,
     4  tempOld, stretchOld, defgradOld, fieldOld,
     3  stressOld, stateOld, enerInternOld, enerInelasOld,
     6  tempNew, stretchNew, defgradNew, fieldNew,
C Write only -
     5  stressNew, stateNew, enerInternNew, enerInelasNew )
C
      use vars_module
      include 'vaba_param.inc'
C
C All arrays dimensioned by (*) are not used in this algorithm
      dimension props(nprops), density(nblock), coordMp(nblock,*),
     1  charLength(nblock), strainInc(nblock,ndir+nshr),
     2  relSpinInc(nblock,nshr), tempOld(nblock),
     3  stretchOld(nblock,ndir+nshr),
     4  defgradOld(nblock,ndir+nshr+nshr),
     5  fieldOld(nblock,nfieldv), stressOld(nblock,ndir+nshr),
     6  stateOld(nblock,nstatev), enerInternOld(nblock),
     7  enerInelasOld(nblock), tempNew(nblock),
     8  stretchNew(nblock,ndir+nshr),
     8  defgradNew(nblock,ndir+nshr+nshr),
     9  fieldNew(nblock,nfieldv),
     1  stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     2  enerInternNew(nblock), enerInelasNew(nblock)
C
      character*80 cmname
C
      parameter( zero = 0., one = 1., two = 2., three = 3.,
     1  third = one/three, half = .5, twoThirds = two/three,
     2  threeHalfs = 1.5 )
C      
      real*8 p, a1, a2, a3
      real*8 lb, phase
      real*8 omega,domega
      real*8 area
      integer jElemUid
C
*     ================================================================================================     *
*     Statev(1)    total strain                                                                            *
*     Statev(2)    total strain                                                                            *
*     Statev(15)   History field, H                                                                        *
*     Statev(16)   phase field, d                                                                          *
*     ================================================================================================     *
*     Parameter initilization
*
      ea   =  30.0d3      ! props(1) -- Young's modulus
      ft   =  4.0d0       ! props(3) -- failure strength
      Gf   =  0.1d0      ! props(4) -- fracture energy
      lb   =  2.0d0       ! props(5) -- length scale            
      thk  =  1.0d0       ! props(6) -- thickness
*
      c0   =  3.1415926535897932384626433832d0
*
      do 100 i = 1,nblock
*     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     *
*       strain update
*          
        strain_old    =  stateOld(i,1) 
        strain_new    =  strain_old + strainInc(i,1)
*        
        jElemUid      =  stateOld(i,8)
*     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     *
        if ((jElemUid .le. 75) .and.(jElemUid .ge. 35)) then
            ft = 2.0d0
        else
            ft = 4.0d0       ! props(3) -- failure strength
        end if
*     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     *
        area = 1.0d0
*                        
        a1   =  4.d0/(c0*lb)*ea*Gf/(ft*ft)
*       concrete softening
        p  =  2.0d0
        a2 =  1.3868d0
        a3 =  0.6567d0
*        
        hist_threshold = 0.5d0*(ft)**2.0d0/ea
*        
        stateNew(i,9) = hist_threshold
*     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     *
*       Trial stress
*         
        sig1   =  ea * strain_new
*
        if (sig1 .gt. 0.0) then
          stressPower = 0.5d0 * sig1 * sig1 / ea
        else
          stressPower = 0.0d0
        end if
*       stressNew(i,1) = sig1
*
        phase = allD(jElemUid)
*       
        call energeticFunc(omega,domega,phase,a1,a2,a3,p)
*     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     *
*        
        if (sig1 .gt. 0.0) then
            sig_pos = sig1 
            sig_neg = 0.0d0
        else
            sig_pos = 0.0d0 
            sig_neg = sig1
        end if
*        
        stressNew(i,1) = (sig_pos * omega + sig_neg) * area    
*
*     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     *
*       update state variables
*
*       total strain, E         
        stateNew(i,1)  =  strain_new
*        
        stateNew(i,2)  =  (sig_pos * omega + sig_neg)
*       stress power      
        stateNew(i,11) =  stressPower
*       histroy variable, H  
        stateNew(i,14) =  max(stateOld(i,14), stateNew(i,11))
*       phase field, d          
        stateNew(i,15) =  phase   
        
        stateNew(i,13) =  ft
        
        allH(jElemUid) = stateNew(i,14)

*     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     *
*       enerInternNew(i) - ALLIE (internal energy)
*
*        enerInternNew(i) = (0.5*omega*sig1*strain_new) / density(i)
*
*       enerInelasNew(i) - ALLPD (plastic dessipation)
*
*        enerInelasNew(i) = (1.0d0 -omega) * (0.5*sig1*strain_new) /
*     +                      density(i)
*
*     ================================================================================================     *
  100 continue
*
      return
      end      
C     *************************************************************************************************
C
C     VUEL for the phase field  
C     Activated DOFs: 11
C
C     *************************************************************************************************
C         
      subroutine vuel(
     *     nblock,
c          to be defined
     *     rhs,amass,dtimeStable,
     *     svars,nsvars,
     *     energy,
c          
     *     nnode,ndofel,
     *     props,nprops,
     *     jprops,njprops,
     *     coords,ncrd,
     *     u,du,v,a,
     *     jtype,jelem,
     *     time,period,dtimeCur,dtimePrev,kstep,kinc,lflags,
     *     dMassScaleFactor,
     *     predef,npredef,
     *     jdltyp,adlmag)
C     
      use vars_module
      include 'vaba_param.inc'
C
C     operation code
      parameter ( jMassCalc            = 1,
     *            jIntForceAndDtStable = 2,
     *            jExternForce         = 3)
C
c     flags
      parameter (iProcedure = 1,
     *           iNlgeom    = 2,
     *           iOpCode    = 3,
     *           nFlags     = 3)
C
c     time
      parameter (iStepTime  = 1,
     *           iTotalTime = 2,
     *           nTime      = 2)
C
c     procedure flags
      parameter ( jDynExplicit = 17 )
C
c     energies 
      parameter ( iElPd = 1,
     *            iElCd = 2,
     *            iElIe = 3,
     *            iElTs = 4,
     *            iElDd = 5,
     *            iElBv = 6,
     *            iElDe = 7,
     *            iElHe = 8,
     *            iElKe = 9,
     *            iElTh = 10,
     *            iElDmd = 11,
     *            iElDc = 12,
     *            nElEnergy = 12)
C
c     predefined variables
      parameter ( iPredValueNew = 1,
     *            iPredValueOld = 2,
     *            nPred         = 2)    
C
c     indexing in a 3-long vector
C
      parameter (factorStable = 0.99d0)
      parameter ( zero = 0.d0, half = 0.5d0, one = 1.d0, two=2.d0 )
C
      dimension rhs(nblock,ndofel), amass(nblock,ndofel,ndofel),
     *     dtimeStable(nblock),
     *     svars(nblock,nsvars), energy(nblock,nElEnergy),
     *     props(nprops), jprops(njprops),
     *     jelem(nblock), time(nTime), lflags(nFlags),
     *     coords(nblock,nnode,ncrd), u(nblock,ndofel),
     *     du(nblock,ndofel), v(nblock,ndofel), a(nblock, ndofel),
     *     dMassScaleFactor(nblock),
     *     predef(nblock, nnode, npredef, nPred), adlmag(nblock)
C      
      dimension UD(2),DUD(2)
C
      integer i,j,k
C
      real*8 Awt(NPT),xii(NPT,1),Bp(NPT), N(2)
C
      real*8  lc,lb, Gf, ea, ft 
      real*8  p, a1, a2, a3
      real*8  phase, dalpha, dphase
      real*8  omega,domega
      real*8  hist, hist_threshold
      real*8  c0
*     ================================================================================================     *
*     Statev(1)    total strain                                                                            *
*     Statev(2)    total strain                                                                            *
*     Statev(15)   History field, H                                                                        *
*     Statev(16)   phase field, d                                                                          *
*     ================================================================================================     *
*     Parameter initilization
*
      eta = 1.0d-4
*
      ea   =  30.0d3      ! props(1) -- Young's modulus
      ft   =  4.0d0       ! props(3) -- failure strength
      Gf   =  0.1d0      ! props(4) -- fracture energy
      lb   =  2.0d0       ! props(5) -- length scale            
      thk  =  1.0d0       ! props(6) -- thickness
*
      c0   =  3.1415926535897932384626433832d0
*     exponential softening
*
      area = 1.0 
*
      if ( lflags(iOpCode).eq.jMassCalc ) then
*     ================================================================================================     *
          do kblock = 1, nblock
*         
              alenX = (coords(kblock,2,1) - coords(kblock,1,1))    !coords(nblock,nnode,ncrd)
              alenY = (coords(kblock,2,2) - coords(kblock,1,2))
              alen  = sqrt(alenX*alenX + alenY*alenY)
*
              amass(kblock,1,1) = eta * alen / 2.0d0
              amass(kblock,2,2) = eta * alen / 2.0d0
*              
          end do         
*     ================================================================================================     *
      else if ( lflags(iOpCode) .eq. jIntForceAndDtStable) then
*     ================================================================================================     *
          WF = Gf/2.0/lb      
          WO = zero    
*              
          do kblock = 1, nblock             
*     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     *
*
             if ((jElemUid .le. 10075) .and.(jElemUid .ge. 10035)) then
                  ft = 2.0d0
             else
                  ft = 4.0d0       ! props(3) -- failure strength
             end if
*         
             a1   =  4.d0/(c0*lb)*ea*Gf/(ft*ft)
*            concrete softening
             p  =  2.0d0
             a2 =  1.3868d0
             a3 =  0.6567d0       
*        
             alenX = (coords(kblock,2,1) - coords(kblock,1,1))    !coords(nblock,nnode,ncrd)
             alenY = (coords(kblock,2,2) - coords(kblock,1,2))
             alen  = sqrt(alenX*alenX + alenY*alenY)
*             
*     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     *
             hist = allH(jelem(kblock) - NumEle)
             if(hist .lt. zero)  hist = zero 

             hist_threshold = 0.5d0*(ft)**2.0d0/ea
*     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     *
             do i = 1, nnode
                UD(i)  = U(kblock,i)
                DUD(i) = DU(kblock,i) 
             end do
*
             phase   = (UD(1) + UD(2)) / 2.0d0
             dphase  = (UD(2) - UD(1)) / alen
             
             if (phase .le. 0.0) then
                 phase = 0.0d0
             end if
*             
             if (phase .ge. 1.0d0) then
                 phase = 1.0d0
                 dphase = 0.0d0
             end if
             
             allD(jelem(kblock) - NumEle) = phase            
             
             call geometricFunc(dalpha,0.0d0) 
             call energeticFunc(omega,domega,0.0d0,a1,a2,a3,p)  
             hist_0 = - Gf * dalpha / domega / c0 / lb
*            
             hist = max(hist, hist_0)             
*
             call geometricFunc(dalpha,phase)                         ! geometric function
             call energeticFunc(omega,domega,phase,a1,a2,a3,p)        ! energetic function
*
             phi_source  = domega *  hist + Gf/(c0*lb)*dalpha
*
*     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     *
             rhs(kblock,1)=   2.d0*lb*Gf/c0*dphase/alen - phi_source
             rhs(kblock,2)= - 2.d0*lb*Gf/c0*dphase/alen - phi_source
*             
             rhs(kblock,1) =  - rhs(kblock,1) * alen
             rhs(kblock,2) =  - rhs(kblock,2) * alen    
*             
*     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     *
*             if(phase .le. 0.3d0) then 
                 dtimeStable(kblock) = factorStable*eta/(500.0*(WF-WO))
*             else
*                 dtimeStable(kblock) = factorStable*eta/(5000.0*(WF-WO))
*             end if
*     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     *
          end do
*     ================================================================================================     *
      end if
c
      return
      end         


C     *************************************************************************************************
C
C     VUSDFLD for converting the state variables in VUMAT into global variables  
C     
C     *************************************************************************************************
C       
      subroutine vusdfld(
c Read only -
     *   nblock, nstatev, nfieldv, nprops, ndir, nshr, 
     *   jElemUid, kIntPt, kLayer, kSecPt, 
     *   stepTime, totalTime, dt, cmname, 
     *   coordMp, direct, T, charLength, props, 
     *   stateOld, 
c Write only -
     *   stateNew, field )
C
      use vars_module
      include 'vaba_param.inc'
*
      dimension props(nprops),
     *          jElemUid(nblock), coordMp(nblock, *), 
     *          direct(nblock, 3, 3), T(nblock,3,3), 
     *          charLength(nblock),
     *          stateOld(nblock, nstatev), 
     *          stateNew(nblock, nstatev),
     *          field(nblock, nfieldv)
      character*80 cmname
*
      character*3 cData(maxblk)
      dimension jData(maxblk)
      dimension eqps(maxblk)
*
      parameter ( zero = 0.d0 )
*
      do k = 1, nblock
            stateOld(k,8)   =  jElemUid(k)
            stateNew(k,8)   =  jElemUid(k)
      end do
*      
      return
      end subroutine
      
       
C     *************************************************************************************************
C
C     VUFIELD for converting the PHASE FIELD to the field adopted in VUMAT  
C     
C     *************************************************************************************************
C       
      subroutine VUField( 
C Write only - 
     *     rUserField, 
C Read only - 
     *     nBlock, nField, kField, nComp,
     *     kStep, kInc, jNodeUid, time, 
     *     coords, U, V, A )
*
      include 'vaba_param.inc'

      dimension rUserField(nBlock,nComp,nField)
      dimension jNodeUid(nBlock), time(4), coords(3,nBlock)
      dimension U(8,nBlock), V(8,nBlock), A(8,nBlock)
*
      parameter ( i_ufld_Current   = 1, 
     *            i_ufld_Increment = 2, 
     *            i_ufld_Period    = 3, 
     *            i_ufld_Total     = 4 )
*
      parameter ( i_ufld_CoordX = 1,
     *            i_ufld_CoordY = 2,
     *            i_ufld_CoordZ = 3 )
*
      parameter ( i_ufld_SpaDisplX = 1,
     *            i_ufld_SpaDisplY = 2,
     *            i_ufld_SpaDisplZ = 3,
     *            i_ufld_RotDisplX = 4,
     *            i_ufld_RotDisplY = 5,
     *            i_ufld_RotDisplZ = 6, 
     *            i_ufld_AcoPress  = 7,
     *            i_ufld_Temp      = 8 )
*
      parameter ( i_ufld_SpaVelX   = 1,
     *            i_ufld_SpaVelY   = 2,
     *            i_ufld_SpaVelZ   = 3,
     *            i_ufld_RotVelX   = 4,
     *            i_ufld_RotVelY   = 5,
     *            i_ufld_RotVelZ   = 6,
     *            i_ufld_DAcoPress = 7,
     *            i_ufld_DTemp     = 8 )
*
      parameter ( i_ufld_SpaAccelX  = 1,
     *            i_ufld_SpaAccelY  = 2,
     *            i_ufld_SpaAccelZ  = 3,
     *            i_ufld_RotAccelX  = 4,
     *            i_ufld_RotAccelY  = 5,
     *            i_ufld_RotAccelZ  = 6, 
     *            i_ufld_DDAcoPress = 7,
     *            i_ufld_DDTemp     = 8 )

      parameter (oneHundred = 100.d0, twoHundred = 200.d0)
*
      if (kField .eq. 1) then
*
         do kComp = 1, nComp
            do kNod = 1, nBlock
               rUserField(kNod,kComp,1) = U(i_ufld_Temp,kNod) 
            end do
         end do
         
      end if
*
      return
      end subroutine
      
      
C     *************************************************************************************************
C
C     energetic degradation function  omega  
C     
C     *************************************************************************************************
C        
      subroutine energeticFunc(omega,domega,phi,a1,a2,a3,p)
*
      include 'vaba_param.inc'
      
      real*8:: omega, domega, phi
      real*8:: fac1, dfac1, ddfac1, fac2, dfac2, ddfac2
      real*8:: p, a1, a2, a3
*     ================================================================================================     *
*      
      fac1    =  (1.d0 - phi)**p
      dfac1   = - p*(1.d0 - phi)**(p - 1.d0); 
      ddfac1  =  p*(p - 1.d0)*(1.d0 - phi)**(p - 2.d0)
*        
      fac2   =  fac1   + a1*phi + a1*a2*phi**2.d0 + a1*a2*a3*phi**3.d0
      dfac2  =  dfac1  + a1 + 2.d0*a1*a2*phi + 3.d0*a1*a2*a3*phi**2.d0
      ddfac2  =  ddfac1 + 2.d0*a1*a2 + 6.d0*a1*a2*a3*phi
*        
      omega   =  fac1/fac2        
      domega  =  (dfac1*fac2  - fac1*dfac2)/(fac2**2.d0)
*      
*      if (phi .gt. 0.9999999d0) then
*          omega  = 0.0
*          domega = 0.0
*      end if    
*     ================================================================================================     *
      return
      end subroutine energeticFunc
      
!**********************************************************************************************************
!
      subroutine geometricFunc(dalpha,phase)
!
!**********************************************************************************************************
      include 'vaba_param.inc'
      real*8  dalpha, phase
*     ================================================================================================     *       
      dalpha  = 2.d0 - 2.d0*phase
*     ================================================================================================     *        
      return 
      end subroutine geometricFunc  

