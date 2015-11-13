#include "helixHelpers.hh"

#include "aidaTT-Units.hh"

namespace aidaTT
{



    void calculateStartHelix(const Vector3D& x1, const Vector3D& x2,   const Vector3D& x3 , trackParameters& tp , bool backward)
    {

      //-------------------------------------------------------------------------------------------------
      // modified version of original code copied from KalTest::THelicalTrack (2003/10/03  K.Fujii  )
      //------------------------------------------------------------------------------------------------
      
      static const Vector3D ez(0., 0., 1.);
      
      Vector3D x12 = x2 - x1 ;
      Vector3D x13 = x3 - x1 ;
      Vector3D x23 = x3 - x2 ;
      
      x12.z() = 0.;
      x13.z() = 0.;
      x23.z() = 0.;
      
      double x12mag = x12.r();
      x12 = x12.unit();
      double x13mag = x13.r();
      x13 = x13.unit();
      double x23mag = x23.r();
      x23 = x23.unit();
      
      double sinHalfPhi23 = x12.cross(x13).z();
      
      double cosHalfPhi23 = 0.5 * (x13mag / x12mag + (1. - x23mag / x12mag) * (x12mag + x23mag) / x13mag);
      double halfPhi23 = std::atan2(sinHalfPhi23, cosHalfPhi23);
      
      double r  = -0.5 * x23mag / sinHalfPhi23;
      
      Vector3D xc = 0.5 * (x2 + x3) + r * cosHalfPhi23 * x23.cross(ez);
      
      if(backward) r = -r;
      
      tp.parameters()(OMEGA)  =  1. / r ;
      tp.parameters()(TANL)  = (x2.z() - x3.z()) / (r * 2 * halfPhi23) ;
      tp.parameters()(PHI0)  =  std::atan2(r * (xc.y() - x1.y()),  r * (xc.x() - x1.x())) + M_PI / 2. ;
      tp.parameters()(D0)  =  0. ;
      tp.parameters()(Z0)  =  0. ;
      
      tp.referencePoint() = x1 ;
      
    }
  
  
  
  double moveHelixTo(trackParameters& tp,  const Vector3D& refNew,  bool CovMat)
  {
    
    //-------------------------------------------------------------------------------------------------
    // modified version of original code copied from KalTest::THelicalTrack (2003/10/03  K.Fujii  )
    //------------------------------------------------------------------------------------------------
    
    // ---------------------------------------------------
    // (0) Preparation
    // ---------------------------------------------------
    //   Define some numerical constants.
    //
    
    static const double kPi    = M_PI ;
    static const double kTwoPi = 2.0 * kPi;

    //   Copy helix parmeters to local variables
    
    double dr    = - tp(D0) ;
    double fi0   =   tp(PHI0) - M_PI / 2. ;
    
    while(fi0 < 0.)      fi0 += kTwoPi;
    while(fi0 > kTwoPi)  fi0 -= kTwoPi;
    
    double cpa   =  tp(OMEGA)  ;    //fixme: ???
    double dz    =  tp(Z0) ;
    double tnl   =  tp(TANL) ;
    
    double x0    = tp.referencePoint().x();
    double y0    = tp.referencePoint().y();
    double z0    = tp.referencePoint().z();
    
    double xv    = refNew.x();
    double yv    = refNew.y();
    double zv    = refNew.z();
    
    // ---------------------------------------------------
    // (1) Calculate a' = f_k-1(a_k-1)
    // ---------------------------------------------------
    //        a' = (dr', fi0', cpa', dz', tnl')
    //        a  = (dr , fi0 , cpa , dz , tnl )
    //
    double r     = 1. / cpa;
    double rdr   = r + dr;
    double csf0  = std::cos(fi0);
    double snf0  = std::sqrt(std::max(0.0, (1.0 - csf0) * (1.0 + csf0)));
    if(fi0 > kPi) snf0 = -snf0;
    
    double xc    = x0 + rdr * csf0;
    double yc    = y0 + rdr * snf0;
    
    // std::cout << " ***++ xc = " << xc << ", yc = " << yc << " XCenter " << calculateXCenter( tp ) << " YCenter " << calculateYCenter( tp ) << std::endl ;
    
    double fi0p  = 0.;
    
    if(cpa > 0.) fi0p = std::atan2((yc - yv), (xc - xv));
    if(cpa < 0.) fi0p = std::atan2((yv - yc), (xv - xc));
    
    while(fi0p < 0.)      fi0p += kTwoPi;
    while(fi0p > kTwoPi)  fi0p -= kTwoPi;
    
    double csf   = std::cos(fi0p);
    double snf   = std::sqrt(std::max(0.0, (1.0 - csf) * (1.0 + csf)));
    if(fi0p > kPi) snf = -snf;
    
    double anrm  = 1.0 / std::sqrt(csf * csf + snf * snf);
    csf  *= anrm;
    snf  *= anrm;
    double csfd  = csf * csf0 + snf * snf0;
    double snfd  = snf * csf0 - csf * snf0;
    
    double fid   = fi0p - fi0;
    while(fid < 0)      fid += kTwoPi;
    while(fid > kTwoPi) fid -= kTwoPi;
    if(fid > kPi)    fid -= kTwoPi;
    
    double drp   = (xc - xv) * csf + (yc - yv) * snf - r;
    double dzp   = z0 - zv + dz - r * tnl * fid;
    
    //fg:
    fi0p += M_PI / 2. ;
    while(fi0p < -kPi)  fi0p += kTwoPi;
    while(fi0p >  kPi)  fi0p -= kTwoPi;
    
    tp(D0)    = - drp ;
    tp(PHI0)  =   fi0p ;
    tp(OMEGA) =   cpa ;
    tp(Z0)    =   dzp ;
    tp(TANL)  =   tnl;
    
    tp.setReferencePoint(refNew) ;
    
    
    if ( CovMat ){

      fiveByFiveMatrix F;
      
      // ---------------------------------------------------
      // (2) Calculate @a'/@a = @a'/a = F_k-1
      // ---------------------------------------------------
      //        a' = (dr', fi0', cpa', dz', tnl')
      //        a  = (dr , fi0 , cpa , dz , tnl )
      //
      
      double rdrpr = 1.0/(r+drp);
      double rcpar = r/cpa;
      
      //YV: mapping the different order of parameters from LCIO to aida 
      
      // @drho'/@a
      F(3,3) = csfd;
      F(3,2) = rdr*snfd;
      F(3,0) = rcpar*(1.0-csfd);
      F(3,4) = 0;
      F(3,1) = 0;
      
      // @phi0'/@a
      F(2,3) = -rdrpr*snfd;
      F(2,2) =  rdr*rdrpr*csfd;
      F(2,0) =  rcpar*rdrpr*snfd;
      F(2,4) =  0;
      F(2,1) =  0;
      
      // @kappa'/@a
      F(0,3) = 0;
      F(0,2) = 0;
      F(0,0) = 1;
      F(0,4) = 0;
      F(0,1) = 0;
      
      // @dz'/@a
      F(4,3) =  r*rdrpr*tnl*snfd;
      F(4,2) =  r*tnl*(1.0-rdr*rdrpr*csfd);
      F(4,0) =  rcpar*tnl*(fid-r*rdrpr*snfd);
      F(4,4) =  1;
      F(4,1) = -r*fid;
      
      // @tanl'/@a
      F(1,3) = 0;
      F(1,2) = 0;
      F(1,0) = 0;
      F(1,4) = 0;
      F(1,1) = 1;
      
      fiveByFiveMatrix C = tp.covarianceMatrix();
      
      fiveByFiveMatrix FT(F);
      FT.Transpose();
      
      fiveByFiveMatrix Cp = F * C * FT ;
      
      tp.setCovarianceMatrix( Cp ) ;
    }
    
    return fid ;
    
    
    
    
    }








}
