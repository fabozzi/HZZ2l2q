#include "HiggsAnalysis/Higgs2l2b/interface/Helicity.h"

#include "TLorentzVector.h"
#include <iostream>


void Helicity::calculateAngles(TLorentzVector thep4H, TLorentzVector thep4Z1, TLorentzVector thep4M11, TLorentzVector thep4M12, TLorentzVector thep4Z2, TLorentzVector thep4M21, TLorentzVector thep4M22, double& costheta1, double& costheta2, double& phi, double& costhetastar, double& phistar1, double& phistar2, double& phistar12, double& phi1, double& phi2, bool &swappedZ1) {

   //std::cout << "In calculate angles..." << std::endl;

   double norm;
   TVector3 boostX = -(thep4H.BoostVector());
   TLorentzVector thep4Z1inXFrame( thep4Z1 );
   TLorentzVector thep4Z2inXFrame( thep4Z2 );       thep4Z1inXFrame.Boost( boostX );
   thep4Z2inXFrame.Boost( boostX );
   TVector3 theZ1X_p3 = TVector3( thep4Z1inXFrame.X(), thep4Z1inXFrame.Y(), thep4Z1inXFrame.Z() );
   TVector3 theZ2X_p3 = TVector3( thep4Z2inXFrame.X(), thep4Z2inXFrame.Y(), thep4Z2inXFrame.Z() );
   // calculate phi1, phi2, costhetastar
   phi1 = theZ1X_p3.Phi();
   phi2 = theZ2X_p3.Phi();
      ///////////////////////////////////////////////
   // check for z1/z2 convention, redefine all 4 vectors with convention
   TLorentzVector p4H, p4Z1, p4M11, p4M12, p4Z2, p4M21, p4M22;
   p4H = thep4H;
   if ((phi1 < 0)&&(phi1 >= -TMath::Pi())){
       // std::cout<<">>> Swapping the Z in calculateAngles (z1/z2 convention)  (phi1="<<phi1<<", phi2="<<phi2<<")"<<std::endl;
       p4Z1 = thep4Z2; p4M11 = thep4M21; p4M12 = thep4M22;
       p4Z2 = thep4Z1; p4M21 = thep4M11; p4M22 = thep4M12;               
       costhetastar = theZ2X_p3.CosTheta();
       swappedZ1=true;
   }
   else{
       p4Z1 = thep4Z1; p4M11 = thep4M11; p4M12 = thep4M12;
       p4Z2 = thep4Z2; p4M21 = thep4M21; p4M22 = thep4M22;
       costhetastar = theZ1X_p3.CosTheta();
       swappedZ1=false;
   }

   //std::cout << "phi1: " << phi1 << ", phi2: " << phi2 << std::endl;
   // now helicity angles................................
   // ...................................................
   TVector3 boostZ1 = -(p4Z1.BoostVector());
   TLorentzVector p4Z2Z1(p4Z2);
   p4Z2Z1.Boost(boostZ1);
   //find the decay axis
   /////TVector3 unitx_1 = -Hep3Vector(p4Z2Z1);
   TVector3 unitx_1( -p4Z2Z1.X(), -p4Z2Z1.Y(), -p4Z2Z1.Z() );
   norm = 1/(unitx_1.Mag());
   unitx_1*=norm;
   //boost daughters of z2
   TLorentzVector p4M21Z1(p4M21);
   TLorentzVector p4M22Z1(p4M22);
   p4M21Z1.Boost(boostZ1);
   p4M22Z1.Boost(boostZ1);
   //create z and y axes
   /////TVector3 unitz_1 = Hep3Vector(p4M21Z1).cross(Hep3Vector(p4M22Z1));
   TVector3 p4M21Z1_p3( p4M21Z1.X(), p4M21Z1.Y(), p4M21Z1.Z() );
   TVector3 p4M22Z1_p3( p4M22Z1.X(), p4M22Z1.Y(), p4M22Z1.Z() );
   TVector3 unitz_1 = p4M21Z1_p3.Cross( p4M22Z1_p3 );
   norm = 1/(unitz_1.Mag());
   unitz_1 *= norm;
   TVector3 unity_1 = unitz_1.Cross(unitx_1);
   //caculate theta1
   TLorentzVector p4M11Z1(p4M11);
   p4M11Z1.Boost(boostZ1);
   TVector3 p3M11( p4M11Z1.X(), p4M11Z1.Y(), p4M11Z1.Z() );
   TVector3 unitM11 = p3M11.Unit();
   double x_m11 = unitM11.Dot(unitx_1); 
   double y_m11 = unitM11.Dot(unity_1); 
   double z_m11 = unitM11.Dot(unitz_1);
   TVector3 M11_Z1frame(y_m11, z_m11, x_m11);
   costheta1 = M11_Z1frame.CosTheta();
   //std::cout << "theta1: " << M11_Z1frame.Theta() << std::endl;
   //////-----------------------old way of calculating phi---------------/////////
   phi = M11_Z1frame.Phi();
   //set axes for other system
   TVector3 boostZ2 = -(p4Z2.BoostVector());
   TLorentzVector p4Z1Z2(p4Z1);
   p4Z1Z2.Boost(boostZ2);
   TVector3 unitx_2( -p4Z1Z2.X(), -p4Z1Z2.Y(), -p4Z1Z2.Z() );
   norm = 1/(unitx_2.Mag());
   unitx_2*=norm;
   //boost daughters of z2
   TLorentzVector p4M11Z2(p4M11);
   TLorentzVector p4M12Z2(p4M12);
   p4M11Z2.Boost(boostZ2);
   p4M12Z2.Boost(boostZ2);
   TVector3 p4M11Z2_p3( p4M11Z2.X(), p4M11Z2.Y(), p4M11Z2.Z() );
   TVector3 p4M12Z2_p3( p4M12Z2.X(), p4M12Z2.Y(), p4M12Z2.Z() );
   TVector3 unitz_2 = p4M11Z2_p3.Cross( p4M12Z2_p3 );
   norm = 1/(unitz_2.Mag());
   unitz_2*=norm;
   TVector3 unity_2 = unitz_2.Cross(unitx_2);
   //calcuate theta2
   TLorentzVector p4M21Z2(p4M21);
   p4M21Z2.Boost(boostZ2);
   TVector3 p3M21( p4M21Z2.X(), p4M21Z2.Y(), p4M21Z2.Z() );
   TVector3 unitM21 = p3M21.Unit();
   double x_m21 = unitM21.Dot(unitx_2); 
   double y_m21 = unitM21.Dot(unity_2); 
   double z_m21 = unitM21.Dot(unitz_2);
   TVector3 M21_Z2frame(y_m21, z_m21, x_m21);
   costheta2 = M21_Z2frame.CosTheta();
      // calculate phi
   //calculating phi_n
   TLorentzVector n_p4Z1inXFrame( p4Z1 );
   TLorentzVector n_p4M11inXFrame( p4M11 );
   n_p4Z1inXFrame.Boost( boostX );
   n_p4M11inXFrame.Boost( boostX );           
   TVector3 n_p4Z1inXFrame_unit = n_p4Z1inXFrame.Vect().Unit();
   TVector3 n_p4M11inXFrame_unit = n_p4M11inXFrame.Vect().Unit();     
   //// y-axis is defined by neg lepton cross z-axis
   //// the subtle part is here...
   TVector3 n_unitz_1( n_p4Z1inXFrame_unit );//z' -> Z1 fly direction 
   TVector3 n_unity_1 = n_p4M11inXFrame_unit.Cross( n_unitz_1 );//vector orthog to neg lept from Z1 and z' (i.e., Z1)-> vect normal to Z1 decay plane
   TVector3 n_unitx_1 = n_unity_1.Cross( n_unitz_1 );//in the Z1 decay plane, points opposite to neg lepton  
      ///////-----------------new way of calculating phi-----------------///////
   TVector3 n_p4PartoninXFrame_unit( 0.0, 0.0, 1.0 );
   TVector3 n_p4PartoninXFrame_unitprime( n_p4PartoninXFrame_unit.Dot(n_unitx_1), n_p4PartoninXFrame_unit.Dot(n_unity_1), n_p4PartoninXFrame_unit.Dot(n_unitz_1) );
   //transform the beam axis to the coordinate system defined by unit_x, unit_y e unit_z as above
   //take the phi of it:phistar1
   phistar1 = n_p4PartoninXFrame_unitprime.Phi();
   // and the calculate phistar2
   TLorentzVector n_p4M21inXFrame( p4M21 );
   n_p4M21inXFrame.Boost( boostX );
   TVector3 n_p4M21inXFrame_unit = n_p4M21inXFrame.Vect().Unit();//3-vector normalized parallel to pos charge lepton from Z1
   //   rotate into other plane
   TVector3 n_p4M21inXFrame_unitprime( n_p4M21inXFrame_unit.Dot(n_unitx_1), n_p4M21inXFrame_unit.Dot(n_unity_1), n_p4M21inXFrame_unit.Dot(n_unitz_1) );
   // std::cout<<"n_p4M21inXFrame_unitprime = ("<< n_p4M21inXFrame_unitprime.X()<<", "<<n_p4M21inXFrame_unitprime.Y() <<", "<<n_p4M21inXFrame_unitprime.Z() <<" )"<<std::endl; 
   TLorentzVector n_p4Z2inXFrame( p4Z2 );
   n_p4Z2inXFrame.Boost( boostX );
   TVector3 n_p4Z2inXFrame_unit = n_p4Z2inXFrame.Vect().Unit();
   //// y-axis is defined by neg lepton cross z-axis
   //// the subtle part is here...
   TVector3 n_unitz_2( n_p4Z2inXFrame_unit );
   TVector3 n_unity_2 = n_p4M21inXFrame_unit.Cross( n_unitz_2 );//lepton #1 from Z2
   TVector3 n_unitx_2 = n_unity_2.Cross( n_unitz_2 );
   TVector3 n_p4PartoninZ2PlaneFrame_unitprime( n_p4PartoninXFrame_unit.Dot(n_unitx_2), n_p4PartoninXFrame_unit.Dot(n_unity_2), n_p4PartoninXFrame_unit.Dot(n_unitz_2) );
   phistar2 = n_p4PartoninZ2PlaneFrame_unitprime.Phi();
   double phistar12_0 = phistar1 + phistar2;
   if (phistar12_0 > TMath::Pi()) phistar12 = phistar12_0 - 2*TMath::Pi();
   else if (phistar12_0 < (-1.)*TMath::Pi()) phistar12 = phistar12_0 + 2*TMath::Pi();
   else phistar12 = phistar12_0;

} //end calculateangles



