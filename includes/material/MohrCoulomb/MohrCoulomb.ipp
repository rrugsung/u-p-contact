mpm::material::MohrCoulomb::MohrCoulomb() {
    setProperty("density", density_);
    setProperty("porosity", porosity_);
    setProperty("permeability", permeability_);
    setProperty("youngModulus", youngs_modulus_);
    setProperty("poissonRatio", poisson_ratio_);
    setProperty("frictionAngle", phi_);
    setProperty("cohesion", coh_);
    setProperty("dilationAngle", psi_);
    setProperty("sigt", sigt_);

    // initialise the plastic deviatoric strain vector
    PDS_ = Eigen::Matrix<double,6,1>::Zero();
    epds_ = 0.;
    phi_ = phi_ * (pi / 180.);
    psi_ = psi_ * (pi / 180.);

    // Bulk modulus
    bulk_modulus_ = youngs_modulus_ / (3.0 * (1. - 2. * poisson_ratio_));
    // Shear modulus
    shear_modulus_ = youngs_modulus_ / (2.0 * (1 + poisson_ratio_));
    this->compute_elastic_stiffness_matrix();
    std::cout << "de: " << de_ ;
    std::cout << "sigt" << sigt_ ;
}

void mpm::material::MohrCoulomb::compute_elastic_stiffness_matrix() {
  // Shear modulus
  const double G = shear_modulus_;
  const double a1 = bulk_modulus_ + (4.0 / 3.0) * G;
  const double a2 = bulk_modulus_ - (2.0 / 3.0) * G;
  // compute elastic stiffness matrix
  // clang-format off
  de_(0,0)=a1;    de_(0,1)=a2;    de_(0,2)=a2;    de_(0,3)=0;    de_(0,4)=0;    de_(0,5)=0;
  de_(1,0)=a2;    de_(1,1)=a1;    de_(1,2)=a2;    de_(1,3)=0;    de_(1,4)=0;    de_(1,5)=0;
  de_(2,0)=a2;    de_(2,1)=a2;    de_(2,2)=a1;    de_(2,3)=0;    de_(2,4)=0;    de_(2,5)=0;
  de_(3,0)= 0;    de_(3,1)= 0;    de_(3,2)= 0;    de_(3,3)=G;    de_(3,4)=0;    de_(3,5)=0;
  de_(4,0)= 0;    de_(4,1)= 0;    de_(4,2)= 0;    de_(4,3)=0;    de_(4,4)=G;    de_(4,5)=0;
  de_(5,0)= 0;    de_(5,1)= 0;    de_(5,2)= 0;    de_(5,3)=0;    de_(5,4)=0;    de_(5,5)=G;

}

//! plane strain
//! Mohr Coulomb model is implemented based on the work by Dilan ROBERT
// param[in] dstrain: strain increment
// param[in] stress: stress
// param[in] epds: equivalent plastic deviatoric strain
// param[out] stress: stress is updated
// param[out] epds: equivalent plastic deviatoric strain is updated
void mpm::material::MohrCoulomb::compute_stress(const Eigen::Matrix<double,1,dof>& dstrain, Eigen::Matrix<double,1,6>& stress, double& epds, double& pc, double& void_ratio) {
  // Get equivalent plastic deviatoric strain
  const double pdstrain = epds;
  VectorD6x1 F = Eigen::Matrix<double,6,1>::Zero();
  VectorD6x1 dS = Eigen::Matrix<double,6,1>::Zero();  
  VectorD6x1 S = stress.transpose();
  F(0) = dstrain(0);
  F(1) = dstrain(1);
  F(3) = dstrain(2);

  double fx, fy, fz, fxy;
  double mohr_c, mohr_r, f1, f2;
  double sigma1, sigma2, sigma3; 
  unsigned mohr_flag; 
  double nphi, npsi, f_s, f_t, alphap, heichi;
  double lamda_s, lamda_t;
  double cs2,si2,dc2,dss; 
  double K, G;
  double a1, a2;
  double phi = phi_;
  double psi = psi_;
  int yield = 0;  

  //-------------------------------------------------------------------------
  // Elastic-predictor stage: compute the trial stress
  S = S + (this->de_ * F);
  fx = S(0);
  fy = S(1);
  fz = S(2);
  fxy = S(3);
   
  mohr_c = 0.5 * (fx + fy); 
  mohr_r = 0.5 * sqrt((fx - fy) * (fx - fy) + 4.0 * fxy * fxy); 
 
  f1 = mohr_c - mohr_r; 
  f2 = mohr_c + mohr_r; 

  if( f1 > fz ){ 
	  sigma1 = fz; 
	  sigma2 = f1; 
	  sigma3 = f2; 
	  mohr_flag = 2; 
  } else if( f2 < fz ){ 
	  sigma1 = f1; 
	  sigma2 = f2; 
	  sigma3 = fz; 
	  mohr_flag = 3; 
  } else { 
	  sigma1 = f1; 
	  sigma2 = fz; 
	  sigma3 = f2; 
	  mohr_flag = 1; 
  }  
    
  //Mohr-Coulomb failure criteria 
  nphi = ( 1.0 + sin( phi ) ) / ( 1.0 - sin( phi ) ); 
  npsi = ( 1.0 + sin( psi ) ) / ( 1.0 - sin( psi ) ); 
  f_s = sigma1 - sigma3 * nphi + 2.0 * coh_ * sqrt( nphi ); 
  f_t = sigt_ - sigma3; 
  alphap = sqrt( 1.0 + nphi * nphi ) + nphi; 
  heichi = sigma3 - sigt_ + alphap * ( sigma1 - nphi * sigt_  
					     + 2.0 * coh_ * sqrt( nphi )); 

  K = youngs_modulus_ / ( 3.0 * ( 1 - 2 * poisson_ratio_ ) );
  G = youngs_modulus_ / ( 2.0 * ( 1 + poisson_ratio_ ) );
  a1 = K + ( 4.0 / 3.0 ) * G;
  a2 = K - ( 2.0 / 3.0 ) * G;
    
  if( heichi < 0. ){ 
	  if( f_s < 0. ){          	                  
	    //tens + shear failure    
	    lamda_s = f_s / ( ( a1 - a2 * npsi ) - ( a2 - a1 * npsi ) * nphi ); 
	    sigma1 -= lamda_s * ( a1 - a2 * npsi ); 
	    sigma2 -= lamda_s * a2 *( 1.0 - npsi ); 
	    sigma3 -= lamda_s * ( a2 - a1 * npsi ); 
      yield = 1;
	  }
  } 
  else if( heichi > 0. ){ 
	  if( f_t < 0. ){          	                  
	    //tens failure 
	    lamda_t = f_t / a1;  
	    sigma1 += lamda_t * a2; 
	    sigma2 += lamda_t * a2; 
	    sigma3 += lamda_t * a1; 
	    sigt_ = 0.0;    
      yield = 1;            
	  } 
  }
    
    // if( f_t < 0. ){          	                  
    // 	//tens + shear failure    
    // 	if( heichi < 0. ){ 
    // 	    lamda_s = f_s / ( ( a1 - a2 * npsi ) - ( a2 - a1 * npsi ) * nphi ); 
    // 	    sigma1 -= lamda_s * ( a1 - a2 * npsi ); 
    // 	    sigma2 -= lamda_s * a2 *( 1.0 - npsi ); 
    // 	    sigma3 -= lamda_s * ( a2 - a1 * npsi ); 
    // 	}else{ 
    // 	    //tens failure 
    // 	    lamda_t = f_t / a1;  
    // 	    sigma1 += lamda_t * a2; 
    // 	    sigma2 += lamda_t * a2; 
    // 	    sigma3 += lamda_t * a1; 
    // 	    sigt_ = 0.0;                
    // 	} 
    // }else if( f_s < 0. ){ 
    // 	lamda_s = f_s / ( ( a1 - a2 * npsi ) - ( a2 - a1 * npsi ) * nphi ); 
    // 	sigma1 -= lamda_s * ( a1 - a2 * npsi ); 
    // 	sigma2 -= lamda_s * a2 * ( 1.0 - npsi ); 
    // 	sigma3 -= lamda_s *( a2 - a1 * npsi ); 
    // }

    // if( sigma3 > 0. ) std::cout<<", sigma: "<<sigma1<<", "<<sigma2<<", "<<sigma3<<", "<<"\n";

  if( sigma1 == sigma3 ){ 
	  cs2 = 1.0; 
    si2 = 0.0; 
  } else{ 
	  cs2 = ( fx - fy ) / ( f1 - f2 ); 
	  si2 = 2.0 * fxy / ( f1 - f2 ); 
  } 
 
  if( mohr_flag == 1 ){ 
	  dc2 = ( sigma1 - sigma3 ) * cs2; 
	  dss = sigma1 + sigma3; 
	  S(0) = 0.5 * ( dss + dc2 );
	  S(1) = 0.5 * ( dss - dc2 ); 
	  S(2) = sigma2; 
	  S(3) = 0.5 * ( sigma1 - sigma3 ) * si2; 
	  S(4) = 0.0; 
	  S(5) = 0.0; 
  } else if( mohr_flag == 2 ){  
	  dc2 = ( sigma2 - sigma3 ) * cs2; 
	  dss = sigma2 + sigma3; 
	  S(0) = 0.5 * ( dss + dc2 ); 
	  S(1) = 0.5 * ( dss - dc2 ); 
	  S(2) = sigma1; 
	  S(3) = 0.5 * ( sigma2 - sigma3 ) * si2; 
	  S(4) = 0.0; 
	  S(5) = 0.0; 
  } else{  
	  dc2 = ( sigma1 - sigma2 ) * cs2; 
    dss = sigma1 + sigma2; 
	  S(0) = 0.5 * ( dss + dc2 ); 
	  S(1) = 0.5 * ( dss - dc2 ); 
	  S(2) = sigma3; 
	  S(3) = 0.5 * ( sigma1 - sigma2 ) * si2; 
	  S(4) = 0.0; 
	  S(5) = 0.0; 
  }
  dS = S - stress.transpose();
  VectorD6x1 pde = Eigen::Matrix<double,6,1>::Zero();
  if(yield == 1)
    pde = F - this -> de_.inverse() * dS;
  //epds = sqrt(2/3*(pde(0)*pde(0)+pde(1)*pde(1)+pde(2)*pde(2)+2*pde(3)*pde(3)+2*pde(4)*pde(4)+2*pde(5)*pde(5)));
  epds = epds + pde(1);
  stress = S.transpose();
}

