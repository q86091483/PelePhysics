c****************************************************************************
c ... DO NOT EDIT.  AUTOMATICALLY GENERATED FILE.
c****************************************************************************
c ...
c ... CKWYP-compatible GETRATES: return molar production rates for all
c ... species given pressure (P), temperature (T), and mass fractions (Y).
c ...
c ... Generated by cgetrates 1.0b
c ...
c****************************************************************************

      subroutine getrates(P,T,Y,ICKWRK,RCKWRK,WSPL)

      use thermchem_m, only : gibbsEnrg_all_dimT    !HK
      use chemkin_m, only : molwt_c                 !HK

      implicit none

      integer    NSMax               ! maximum number of species
      parameter (NSMax = 1000)
      integer    NRMax               ! maximum number of reactions
      parameter (NRMax = 5000)

      real PA            ! atm. pressure; dyne/cm^2
      parameter (PA  = 1.013250d+06)
      real R0            ! gas constant; erg/(mol K)
      parameter (R0  = 8.314510d+07)
      real R0c           ! gas constant; cal/(mol K)
      parameter (R0c = 1.9872155832)
      real DLn10         ! log(10)
      parameter (DLn10 = 2.3025850929940459d0)

      integer          ICKWRK(*)     ! ...
      real RCKWRK(*)     ! ...
      real P             ! pressure
      real T             ! temperature (K)
      real Y(NSMax)      ! ...
      real WSPL(NSMax)   ! ...

      real FCent         ! ...
      real FDenom        ! ...
      real FLogPR        ! ...
      real FQuan         ! ...
      real CTOT          ! ...
      real OPRT          ! ...
      real ORTC          ! ...
      real OTC           ! ...
      real PR            ! ...
      real PRT           ! ...
      real RC            ! ...
      real RR_F          ! forward reaction rate
      real RR_R          ! reverse reaction rate
      real RR_K0         ! ...
      real RR_KInf       ! ...
      real TC            ! ...
      real THBCTEMP      ! ...
      real VLNPRT        ! ...
      real VLNTEMP       ! ...
      real XIK           ! ...

      real CGSPL(NSMax)  ! ...
      real CSPL(NSMax)   ! ...
      real ROPL(NRMax)   ! ...

      real or0tc         ! HK

      integer          KSP           ! ...
      integer          I             ! loop index

      intrinsic exp, log
      external  CKRHOY, CKYTCP

      dimension thbctemp(4)

      integer NC_H2, NC_O2, NC_O, NC_OH, NC_H2O, NC_H, NC_HO2, NC_H2O2
      integer NC_CO, NC_CO2, NC_HCO, NC_N2

      parameter (NC_H2                =   1)
      parameter (NC_O2                =   2)
      parameter (NC_O                 =   3)
      parameter (NC_OH                =   4)
      parameter (NC_H2O               =   5)
      parameter (NC_H                 =   6)
      parameter (NC_HO2               =   7)
      parameter (NC_H2O2              =   8)
      parameter (NC_CO                =   9)
      parameter (NC_CO2               =  10)
      parameter (NC_HCO               =  11)
      parameter (NC_N2                =  12)

      integer    NS                 ! number of species
      parameter (NS      =   12)

c----------------------------------------------------------------------------
c ... loop over all cells

c$doacross local(i,j,tc,rc,otc,vlntemp,cgspl,ksp,wspl,cspl,rr_f,rr_r,
c$&              xik,ctot,ropl,pr,fcent,thbctemp,rr_k0,rr_kinf,flogpr,
c$&              fdenom,fquan,prt,vlnprt),
c$&        share(m,n,temp,rho,temp1tab,odtab,ns,cgsptab,ys,ospwt,spwt,
c$&              temptab,wchem,tempref)

cmic$  do all
cmic$* private(xik,ctot,ropl,pr,fcent,thbctemp,rr_k0,rr_kinf,flogpr)
cmic$* private(i,j,tc,rc,otc,vlntemp,cgspl,ksp,wspl,cspl,rr_f,rr_r)
cmic$* private(fdenom,fquan,prt,vlnprt)
cmic$* shared(m,n,temp,rho,temp1tab,odtab,ns,cgsptab,ys,ospwt,spwt)
cmic$* shared(temptab,wchem,tempref)

c----------------------------------------------------------------------------

      TC = T

      call CKRHOY(P,T,Y,ICKWRK,RCKWRK,RC)

      otc     = 1.0 / tc
      ortc    = 1.0 / (tc * R0c)
      vlntemp = log(TC)
      or0tc   = 1.0/(R0 * TC)
      prt     = PA *or0tc
!      prt     = PA / (R0 * TC)
      oprt    = 1.0 / prt
      vlnprt  = log(prt)

c ... get gibbs free energies

c inlined
c     call getcgsp(i,j,cgspl)

c      call CKGML(T,ICKWRK,RCKWRK,CGSPL)
c      do I=1,NS
c         CGSPL(I) = CGSPL(I) / R0
c      end do
c      evatt - this is faster
      call gibbsEnrg_all_dimT(T, cgspl(:))

c ... initialize

c     do ksp=1,ns
c        cspl(ksp)  = rc*ys(i,j,ksp)*ospwt(ksp)           ! concentrations
c     end do

c      call CKYTCP(P,T,Y,ICKWRK,RCKWRK,CSPL)
c      evatt - this is faster
      cspl(1:ns) = y(1:ns)*molwt_c(1:ns)
      rc = p*or0tc/sum(cspl(1:ns))
      cspl(1:ns) = rc*cspl(1:ns) 

c      prt = pa / (r0 * tc)
c      vlnprt = log(prt)
c----------------------------------------
c ... calculate third-body concentrations

      ctot = 0.0
      do ksp = 1, ns
        ctot = ctot + cspl(ksp)
      end do
      thbctemp(2) = ctot + 1.5d+000*cspl(NC_H2) + 1.1d+001*cspl(NC_H2O)
     &   + 8.999999999999999d-001*cspl(NC_CO) + 2.8d+000*cspl(NC_CO2)
      thbctemp(3) = ctot + cspl(NC_H2) - 2.2d-001*cspl(NC_O2) + 1.0d+001
     &  *cspl(NC_H2O) + 8.999999999999999d-001*cspl(NC_CO) + 2.8d+000*
     &  cspl(NC_CO2)
      thbctemp(4) = ctot + 1.5d+000*cspl(NC_H2) + 5.0d+000*cspl(NC_H2O)
     &   + 8.999999999999999d-001*cspl(NC_CO) + 2.8d+000*cspl(NC_CO2)
      
c----------------------------------------------------
c ... loop over all rxns to evaluate rate-of-progress

c   1)  O2 + H <=> OH + O
      rr_f = 3.547d+015 * exp(-4.06d-001*vlntemp - 1.6599d+004*ortc)
      xik = -cgspl(NC_O2) - cgspl(NC_H) + cgspl(NC_OH) + cgspl(NC_O)
      rr_r = rr_f * exp(xik*otc)
      ropl(1) = rr_f*cspl(NC_O2)*cspl(NC_H) - rr_r*cspl(NC_OH)
     &  *cspl(NC_O)
      
c   2)  H2 + O <=> OH + H
      rr_f = 5.08d+004 * exp(2.67d+000*vlntemp - 6.29d+003*ortc)
      xik = -cgspl(NC_H2) - cgspl(NC_O) + cgspl(NC_OH) + cgspl(NC_H)
      rr_r = rr_f * exp(xik*otc)
      ropl(2) = rr_f*cspl(NC_H2)*cspl(NC_O) - rr_r*cspl(NC_OH)
     &  *cspl(NC_H)
      
c   3)  OH + H2 <=> H + H2O
      rr_f = 2.16d+008 * exp(1.51d+000*vlntemp - 3.43d+003*ortc)
      xik = -cgspl(NC_OH) - cgspl(NC_H2) + cgspl(NC_H) + cgspl(NC_H2O)
      rr_r = rr_f * exp(xik*otc)
      ropl(3) = rr_f*cspl(NC_OH)*cspl(NC_H2) - rr_r*cspl(NC_H)
     &  *cspl(NC_H2O)
      
c   4)  H2O + O <=> 2 OH
      rr_f = 2.97d+006 * exp(2.02d+000*vlntemp - 1.34d+004*ortc)
      xik = -cgspl(NC_H2O) - cgspl(NC_O) + 2.0 * cgspl(NC_OH)
      rr_r = rr_f * exp(xik*otc)
      ropl(4) = rr_f*cspl(NC_H2O)*cspl(NC_O) - rr_r*cspl(NC_OH)
     &  *cspl(NC_OH)
      
c   5)  H2 + M <=> 2 H + M 
      rr_f = 4.577d+019 * exp(-1.4d+000*vlntemp - 1.0438d+005*ortc)
      xik = -cgspl(NC_H2) + 2.0 * cgspl(NC_H)
      rr_r = rr_f * exp(xik*otc) * oprt
      ropl(5) = rr_f*cspl(NC_H2) - rr_r*cspl(NC_H)*cspl(NC_H)
      ropl(5) = ropl(5) * thbctemp(2)
      
c   6)  2 O + M <=> O2 + M 
      rr_f = 6.165d+015 * exp(-5.0d-001*vlntemp)
      xik = -2.0 * cgspl(NC_O) + cgspl(NC_O2)
      rr_r = rr_f * exp(xik*otc) * prt
      ropl(6) = rr_f*cspl(NC_O)*cspl(NC_O) - rr_r*cspl(NC_O2)
      ropl(6) = ropl(6) * thbctemp(2)
      
c   7)  H + O + M <=> OH + M 
      rr_f = 4.714d+018 * otc
      xik = -cgspl(NC_H) - cgspl(NC_O) + cgspl(NC_OH)
      rr_r = rr_f * exp(xik*otc) * prt
      ropl(7) = rr_f*cspl(NC_H)*cspl(NC_O) - rr_r*cspl(NC_OH)
      ropl(7) = ropl(7) * thbctemp(2)
      
c   8)  OH + H + M <=> H2O + M 
      rr_f = 3.8d+022 * otc * otc
      xik = -cgspl(NC_OH) - cgspl(NC_H) + cgspl(NC_H2O)
      rr_r = rr_f * exp(xik*otc) * prt
      ropl(8) = rr_f*cspl(NC_OH)*cspl(NC_H) - rr_r*cspl(NC_H2O)
      ropl(8) = ropl(8) * thbctemp(2)
      
c   9)  O2 + H (+M) <=> HO2 (+M) 
      rr_k0 = 6.366d+020 * exp(-1.72d+000*vlntemp - 5.248d+002*ortc)
      rr_kinf = 1.475d+012 * exp(6.0d-001*vlntemp)
      pr = rr_k0 / rr_kinf * thbctemp(3)
c     evatt - fcent always log10(0.8) 
c      fcent = log10(2.0d-001 * exp(-9.999999999999999d+029 * tc) + 
c     &  8.0d-001 * exp(-9.999999999999999d-031 * tc))
      fcent = -9.691001300e-02
      flogpr = log10(pr) - 0.4 - 0.67 * fcent
      fdenom = 0.75 - 1.27 * fcent - 0.14 * flogpr
      fquan = flogpr / fdenom
      fquan = fcent / (1.0 + fquan * fquan)
      rr_f = rr_kinf * pr/(1.0 + pr) * exp(fquan*dln10)
      xik = -cgspl(NC_O2) - cgspl(NC_H) + cgspl(NC_HO2)
      rr_r = rr_f * exp(xik*otc) * prt
      ropl(9) = rr_f*cspl(NC_O2)*cspl(NC_H) - rr_r*cspl(NC_HO2)
      
c  10)  H + HO2 <=> O2 + H2
      rr_f = 1.66d+013 * exp(-8.23d+002*ortc)
      xik = -cgspl(NC_H) - cgspl(NC_HO2) + cgspl(NC_O2) + cgspl(NC_H2)
      rr_r = rr_f * exp(xik*otc)
      ropl(10) = rr_f*cspl(NC_H)*cspl(NC_HO2) - rr_r*cspl(NC_O2)
     &  *cspl(NC_H2)
      
c  11)  H + HO2 <=> 2 OH
      rr_f = 7.079d+013 * exp(-2.95d+002*ortc)
      xik = -cgspl(NC_H) - cgspl(NC_HO2) + 2.0 * cgspl(NC_OH)
      rr_r = rr_f * exp(xik*otc)
      ropl(11) = rr_f*cspl(NC_H)*cspl(NC_HO2) - rr_r*cspl(NC_OH)
     &  *cspl(NC_OH)
      
c  12)  O + HO2 <=> OH + O2
      rr_f = 3.25d+013
      xik = -cgspl(NC_O) - cgspl(NC_HO2) + cgspl(NC_OH) + cgspl(NC_O2)
      rr_r = rr_f * exp(xik*otc)
      ropl(12) = rr_f*cspl(NC_O)*cspl(NC_HO2) - rr_r*cspl(NC_OH)
     &  *cspl(NC_O2)
      
c  13)  OH + HO2 <=> O2 + H2O
      rr_f = 2.89d+013 * exp(4.97d+002*ortc)
      xik = -cgspl(NC_OH) - cgspl(NC_HO2) + cgspl(NC_O2) + cgspl(NC_H2O)
     &  
      rr_r = rr_f * exp(xik*otc)
      ropl(13) = rr_f*cspl(NC_OH)*cspl(NC_HO2) - rr_r*cspl(NC_O2)
     &  *cspl(NC_H2O)
      
c  14, 15)  2 HO2 <=> O2 + H2O2
      rr_f = 4.2d+014 * exp(-1.1982d+004*ortc)
      rr_f = rr_f + 1.3d+011 * exp(1.6293d+003*ortc)
      xik = -2.0 * cgspl(NC_HO2) + cgspl(NC_O2) + cgspl(NC_H2O2)
      rr_r = rr_f * exp(xik*otc)
      ropl(14) = rr_f*cspl(NC_HO2)*cspl(NC_HO2) - rr_r*cspl(NC_O2)
     &  *cspl(NC_H2O2)
      
c  16)  H2O2 (+M) <=> 2 OH (+M) 
      rr_k0 = 1.202d+017 * exp(-4.55d+004*ortc)
      rr_kinf = 2.951d+014 * exp(-4.843d+004*ortc)
      pr = rr_k0 / rr_kinf * thbctemp(2)
c     evatt - this is always equal to log10(0.5)
c      fcent = log10(5.0d-001 * exp(-9.999999999999999d+029 * tc) + 
c     &  5.0d-001 * exp(-9.999999999999999d-031 * tc))
      fcent = -3.0102999566398e-01
      flogpr = log10(pr) - 0.4 - 0.67 * fcent
      fdenom = 0.75 - 1.27 * fcent - 0.14 * flogpr
      fquan = flogpr / fdenom
      fquan = fcent / (1.0 + fquan * fquan)
      rr_f = rr_kinf * pr/(1.0 + pr) * exp(fquan*dln10)
      xik = -cgspl(NC_H2O2) + 2.0 * cgspl(NC_OH)
      rr_r = rr_f * exp(xik*otc) * oprt
      ropl(16) = rr_f*cspl(NC_H2O2) - rr_r*cspl(NC_OH)*cspl(NC_OH)
      
c  17)  H + H2O2 <=> OH + H2O
      rr_f = 2.41d+013 * exp(-3.97d+003*ortc)
      xik = -cgspl(NC_H) - cgspl(NC_H2O2) + cgspl(NC_OH) + cgspl(NC_H2O)
     &  
      rr_r = rr_f * exp(xik*otc)
      ropl(17) = rr_f*cspl(NC_H)*cspl(NC_H2O2) - rr_r*cspl(NC_OH)
     &  *cspl(NC_H2O)
      
c  18)  H + H2O2 <=> H2 + HO2
      rr_f = 4.82d+013 * exp(-7.95d+003*ortc)
      xik = -cgspl(NC_H) - cgspl(NC_H2O2) + cgspl(NC_H2) + cgspl(NC_HO2)
     &  
      rr_r = rr_f * exp(xik*otc)
      ropl(18) = rr_f*cspl(NC_H)*cspl(NC_H2O2) - rr_r*cspl(NC_H2)
     &  *cspl(NC_HO2)
      
c  19)  O + H2O2 <=> HO2 + OH
      rr_f = 9.55d+006 * tc * tc * exp(-3.97d+003*ortc)
      xik = -cgspl(NC_O) - cgspl(NC_H2O2) + cgspl(NC_HO2) + cgspl(NC_OH)
     &  
      rr_r = rr_f * exp(xik*otc)
      ropl(19) = rr_f*cspl(NC_O)*cspl(NC_H2O2) - rr_r*cspl(NC_HO2)
     &  *cspl(NC_OH)
      
c  20, 21)  OH + H2O2 <=> H2O + HO2
      rr_f = 1.0d+012
      rr_f = rr_f + 5.8d+014 * exp(-9.557d+003*ortc)
      xik = -cgspl(NC_OH) - cgspl(NC_H2O2) + cgspl(NC_H2O) + 
     &  cgspl(NC_HO2)
      rr_r = rr_f * exp(xik*otc)
      ropl(20) = rr_f*cspl(NC_OH)*cspl(NC_H2O2) - rr_r*cspl(NC_H2O)
     &  *cspl(NC_HO2)
      
c  22)  O + CO (+M) <=> CO2 (+M) 
      rr_k0 = 1.55d+024 * exp(-2.79d+000*vlntemp - 4.191d+003*ortc)
      rr_kinf = 1.8d+010 * exp(-2.384d+003*ortc)
      pr = rr_k0 / rr_kinf * thbctemp(2)
      rr_f = rr_kinf * pr/(1.0 + pr)
      xik = -cgspl(NC_O) - cgspl(NC_CO) + cgspl(NC_CO2)
      rr_r = rr_f * exp(xik*otc) * prt
      ropl(22) = rr_f*cspl(NC_O)*cspl(NC_CO) - rr_r*cspl(NC_CO2)
      
c  23)  O2 + CO <=> O + CO2
      rr_f = 2.53d+012 * exp(-4.77d+004*ortc)
      xik = -cgspl(NC_O2) - cgspl(NC_CO) + cgspl(NC_O) + cgspl(NC_CO2)
      rr_r = rr_f * exp(xik*otc)
      ropl(23) = rr_f*cspl(NC_O2)*cspl(NC_CO) - rr_r*cspl(NC_O)
     &  *cspl(NC_CO2)
      
c  24)  HO2 + CO <=> OH + CO2
      rr_f = 3.01d+013 * exp(-2.3d+004*ortc)
      xik = -cgspl(NC_HO2) - cgspl(NC_CO) + cgspl(NC_OH) + cgspl(NC_CO2)
     &  
      rr_r = rr_f * exp(xik*otc)
      ropl(24) = rr_f*cspl(NC_HO2)*cspl(NC_CO) - rr_r*cspl(NC_OH)
     &  *cspl(NC_CO2)
      
c  25)  OH + CO <=> H + CO2
      rr_f = 2.229d+005 * exp(1.89d+000*vlntemp + 1.1587d+003*ortc)
      xik = -cgspl(NC_OH) - cgspl(NC_CO) + cgspl(NC_H) + cgspl(NC_CO2)
      rr_r = rr_f * exp(xik*otc)
      ropl(25) = rr_f*cspl(NC_OH)*cspl(NC_CO) - rr_r*cspl(NC_H)
     &  *cspl(NC_CO2)
      
c  26)  HCO + M <=> CO + H + M 
      rr_f = 4.7485d+011 * exp(6.59d-001*vlntemp - 1.4874d+004*ortc)
      xik = -cgspl(NC_HCO) + cgspl(NC_CO) + cgspl(NC_H)
      rr_r = rr_f * exp(xik*otc) * oprt
      ropl(26) = rr_f*cspl(NC_HCO) - rr_r*cspl(NC_CO)*cspl(NC_H)
      ropl(26) = ropl(26) * thbctemp(4)
      
c  27)  O2 + HCO <=> HO2 + CO
      rr_f = 7.58d+012 * exp(-4.1d+002*ortc)
      xik = -cgspl(NC_O2) - cgspl(NC_HCO) + cgspl(NC_HO2) + cgspl(NC_CO)
     &  
      rr_r = rr_f * exp(xik*otc)
      ropl(27) = rr_f*cspl(NC_O2)*cspl(NC_HCO) - rr_r*cspl(NC_HO2)
     &  *cspl(NC_CO)
      
c  28)  H + HCO <=> H2 + CO
      rr_f = 7.23d+013
      xik = -cgspl(NC_H) - cgspl(NC_HCO) + cgspl(NC_H2) + cgspl(NC_CO)
      rr_r = rr_f * exp(xik*otc)
      ropl(28) = rr_f*cspl(NC_H)*cspl(NC_HCO) - rr_r*cspl(NC_H2)
     &  *cspl(NC_CO)
      
c  29)  O + HCO <=> H + CO2
      rr_f = 3.0d+013
      xik = -cgspl(NC_O) - cgspl(NC_HCO) + cgspl(NC_H) + cgspl(NC_CO2)
      rr_r = rr_f * exp(xik*otc)
      ropl(29) = rr_f*cspl(NC_O)*cspl(NC_HCO) - rr_r*cspl(NC_H)
     &  *cspl(NC_CO2)
      
c-----------------------------------------------------
c ... evaluate contributions to reactants and products

c  1. H2
      wspl(NC_H2) = -ropl(2)-ropl(3)-ropl(5)+ropl(10)+ropl(18)+ropl(28)
c  2. O2
      wspl(NC_O2) = -ropl(1)+ropl(6)-ropl(9)+ropl(10)+ropl(12)+ropl(13)
     &  +ropl(14)-ropl(23)-ropl(27)
c  3. O
      wspl(NC_O) = ropl(1)-ropl(2)-ropl(4)-2.0*ropl(6)-ropl(7)-ropl(12)
     &  -ropl(19)-ropl(22)+ropl(23)-ropl(29)
c  4. OH
      wspl(NC_OH) = ropl(1)+ropl(2)-ropl(3)+2.0*ropl(4)+ropl(7)-ropl(8)
     &  +2.0*ropl(11)+ropl(12)-ropl(13)+2.0*ropl(16)+ropl(17)+ropl(19)
     &  -ropl(20)+ropl(24)-ropl(25)
c  5. H2O
      wspl(NC_H2O) = ropl(3)-ropl(4)+ropl(8)+ropl(13)+ropl(17)+ropl(20)
c  6. H
      wspl(NC_H) = -ropl(1)+ropl(2)+ropl(3)+2.0*ropl(5)-ropl(7)-ropl(8)
     &  -ropl(9)-ropl(10)-ropl(11)-ropl(17)-ropl(18)+ropl(25)+ropl(26)
     &  -ropl(28)+ropl(29)
c  7. HO2
      wspl(NC_HO2) = ropl(9)-ropl(10)-ropl(11)-ropl(12)-ropl(13)
     &  -2.0*ropl(14)+ropl(18)+ropl(19)+ropl(20)-ropl(24)+ropl(27)
c  8. H2O2
      wspl(NC_H2O2) = ropl(14)-ropl(16)-ropl(17)-ropl(18)-ropl(19)
     &  -ropl(20)
c  9. CO
      wspl(NC_CO) = -ropl(22)-ropl(23)-ropl(24)-ropl(25)+ropl(26)
     &  +ropl(27)+ropl(28)
c 10. CO2
      wspl(NC_CO2) = ropl(22)+ropl(23)+ropl(24)+ropl(25)+ropl(29)
c 11. HCO
      wspl(NC_HCO) = -ropl(26)-ropl(27)-ropl(28)-ropl(29)
c 12. N2
      wspl(NC_N2) = 0.0
      
      return
      end

