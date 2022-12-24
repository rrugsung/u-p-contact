mpm::material::ModifiedCamClay::ModifiedCamClay() {
    setProperty("density", density_);
    setProperty("porosity", porosity_);
    setProperty("permeability", permeability_);
    setProperty("youngModulus", youngs_modulus_);
    setProperty("poissonRatio", poisson_ratio_);
    setProperty("e_ref", e_ref_);
    setProperty("p_ref", p_ref_);
    setProperty("ocr", ocr_);
    setProperty("pc0", pc0_);
    setProperty("m", m_);
    setProperty("lambda", lambda_);
    setProperty("kappa", kappa_);

    // initialise parameters
    PDS_ = Eigen::Matrix<double,6,1>::Zero();
    epds_ = 0.;

    dpvstrain_ = 0;
    dpdstrain_ = 0;
    pvstrain_ = 0;
    pdstrain_ = 0;
    
    pc_ = pc0_;
    
    delta_phi_ = 0;
    m_theta_ = m_;
    f_function_ = 0;
    subloading_r_ = 1;

    e0_ = e_ref_ - lambda_ * log(pc0_/ ocr_ / p_ref_) - kappa_ * log(ocr_);
    void_ratio_ = e0_;
    // Bulk modulus
    bulk_modulus_ = youngs_modulus_ / (3.0 * (1. - 2. * poisson_ratio_));
    // Shear modulus
    shear_modulus_ = youngs_modulus_ / (2.0 * (1 + poisson_ratio_));
    this->compute_elastic_stiffness_matrix();
}

void mpm::material::ModifiedCamClay::compute_elastic_stiffness_matrix() {
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

void mpm::material::ModifiedCamClay::compute_elastic_tensor() {
  // Compute elastic modulus based on stress status
  if (p_ > 1e-08) {
    // Bulk modulus
    bulk_modulus_ =
        (1 + void_ratio_) / kappa_ * p_;
    // Shear modulus
    shear_modulus_ = 3 * bulk_modulus_ *
                                        (1 - 2 * poisson_ratio_) /
                                        (2 * (1 + poisson_ratio_));
  }
  // Components in stiffness matrix
  const double G = shear_modulus_;
  const double a1 = bulk_modulus_ + (4.0 / 3.0) * G;
  const double a2 = bulk_modulus_ - (2.0 / 3.0) * G;
  //std::cout << "bulk modulus: " << bulk_modulus_ << "\n";
  //std::cout << "shear_modulus: " << shear_modulus_ << "\n";

  // Compute elastic stiffness matrix
  de_(0,0)=a1;    de_(0,1)=a2;    de_(0,2)=a2;    de_(0,3)=0;    de_(0,4)=0;    de_(0,5)=0;
  de_(1,0)=a2;    de_(1,1)=a1;    de_(1,2)=a2;    de_(1,3)=0;    de_(1,4)=0;    de_(1,5)=0;
  de_(2,0)=a2;    de_(2,1)=a2;    de_(2,2)=a1;    de_(2,3)=0;    de_(2,4)=0;    de_(2,5)=0;
  de_(3,0)= 0;    de_(3,1)= 0;    de_(3,2)= 0;    de_(3,3)=G;    de_(3,4)=0;    de_(3,5)=0;
  de_(4,0)= 0;    de_(4,1)= 0;    de_(4,2)= 0;    de_(4,3)=0;    de_(4,4)=G;    de_(4,5)=0;
  de_(5,0)= 0;    de_(5,1)= 0;    de_(5,2)= 0;    de_(5,3)=0;    de_(5,4)=0;    de_(5,5)=G;
}

void mpm::material::ModifiedCamClay::compute_stress_invariants(const VectorD6x1 stress) {
  p_ = -mpm::material::p(stress);
  q_ = mpm::material::q(stress);
  theta_ = mpm::material::lode_angle(stress);
}

mpm::material::ModifiedCamClay::FailureState
    mpm::material::ModifiedCamClay::compute_yield_state() {
  double p = p_;
  double q = q_;
  double pc = pc_;
  double m_theta = m_theta_;
  double subloading_r = subloading_r_;
  // Initialise yield status (0: elastic, 1: yield)
  auto yield_type = FailureState::Elastic;
  // Compute yield functions
  f_function_ = std::pow(q / m_theta, 2) + p * (p - subloading_r * pc);
  // Tension failure
  if (f_function_ > 1E-08)
    yield_type = FailureState::Yield;
  return yield_type;
}

//! Compute dF/dmul
void mpm::material::ModifiedCamClay::compute_df_dmul(double* df_dmul) {
  double p = p_;
  double q = q_;
  double m_theta = m_theta_;
  double pc = pc_;
  double mul = delta_phi_;
  double e_b = bulk_modulus_;
  double e_s = shear_modulus_;
  double void_ratio = void_ratio_;
  // Compute dF / dp
  double df_dp = 2 * p - pc;
  // Compute dF / dq
  double df_dq = 2 * q / std::pow(m_theta, 2);
  // Compute dF / dpc
  double df_dpc = -p;
  // Upsilon
  double upsilon = (1 + void_ratio) / (lambda_ - kappa_);
  // A_den
  double a_den = 1 + (2 * e_b + upsilon * pc) * mul;
  // Compute dp / dmul
  double dp_dmul = -e_b * (2 * p - pc) / a_den;
  // Compute dpc / dmul
  double dpc_dmul = upsilon * (pc) * (2 * p - pc) / a_den;
  // Compute dq / dmul
  double dq_dmul = -q / (mul + std::pow(m_theta, 2) / (6 * e_s));
  // Compute dF / dmul
  *df_dmul = (df_dp * dp_dmul) + (df_dq * dq_dmul) + (df_dpc * dpc_dmul);
}

//! Compute dg/dpc
void mpm::material::ModifiedCamClay::compute_dg_dpc(
    const double pc_n, const double p_trial,
    double* g_function, double* dg_dpc) {
  // Upsilon
  const double upsilon =
      (1 + void_ratio_) / (lambda_ - kappa_);
  // Exponential index
  double e_index =
      upsilon * delta_phi_ *
      (2 * p_trial - pc_) /
      (1 + 2 * delta_phi_ * bulk_modulus_);
  // Compute consistency parameter function
  (*g_function) = pc_n * exp(e_index) - pc_;
  // Compute dG / dpc
  (*dg_dpc) = pc_n * exp(e_index) * (-upsilon * delta_phi_ / (1 + 2 * delta_phi_ * bulk_modulus_)) - 1;
}

//Compute Stress
void mpm::material::ModifiedCamClay::compute_stress(const Eigen::Matrix<double,1,dof>& dstrain, Eigen::Matrix<double,1,6>& stress, double& epds, double& pc, double& void_ratio) {
  // Get equivalent plastic deviatoric strain
  VectorD6x1 F = Eigen::Matrix<double,6,1>::Zero();
  VectorD6x1 dS = Eigen::Matrix<double,6,1>::Zero();  
  VectorD6x1 S = stress.transpose();
  //std::cout << "Initial Stress: " << S << "\n";
  F(0) = dstrain(0);
  F(1) = dstrain(1);
  F(3) = dstrain(2);
  //F(0) = 0.0005;
  //F(1) = 0.0005;
  //F(2) = -0.001;

  pc_ = pc;
  void_ratio_ = void_ratio;
  //std::cout << "pc0: " << pc_ << "\n";
  //std::cout << "void_ratio: " << void_ratio_ << "\n";


  // Tolerance for yield function
  const double Ftolerance = 1.E-5;
  // Tolerance for preconsolidation function
  const double Gtolerance = 1.E-5;
  // Maximum iteration step number
  const int itrstep = 100;
  // Maximum subiteration step number
  const int substep = 100;
  // initial mean stress
  p_ = -(stress(0)+stress(1)+stress(2))/3;
  //std::cout << "p': " << p_ << "\n";
  this -> compute_elastic_tensor();
  //std::cout <<"dstrain" << F << "\n";
  //std::cout <<"de_" << de_ << "\n";  
  //-------------------------------------------------------------------------
  // Elastic step
  // Compute trial stress
  VectorD6x1 trial_stress = S + (this->de_ *F);
  // initialise vector n
  VectorD6x1 n_trial = VectorD6x1::Zero();
  // compute trial stress invariants
  this->compute_stress_invariants(trial_stress);
  // compute deviatoric stress tensor
  if (q_ > 1E-08)
    n_trial = mpm::material::deviatoric_stress(trial_stress) / q_;
  // check yield status
  auto yield_type = this->compute_yield_state();
  // Return the updated stress in elastic state
  //std::cout << "yield_type: " << yield_type << "\n";
  if(yield_type == FailureState::Elastic) {
    stress = trial_stress.transpose();
    //std::cout << "stress " << stress << "\n";
    // Update void_ratio
    //void_ratio +=
    //  ((dstrain(0) + dstrain(1) + dstrain(2)) * (1 + e0_));
  }
  //---------------------------------------------------------------------------
  // Plastic step
  else {
    // Counters for interations
    int counter_f = 0;
    int counter_g = 0;
    // Initialise consistency parameter
    delta_phi_ = 0;
    // Volumetric trial stress
    const double p_trial = p_;
    // Deviatoric trial stress
    const double q_trial = q_;
    // M_theta of trial stress
    const double m_theta_trial = m_theta_;
    // Preconsolidation pressure of last step
    const double pc_n = pc_;
    // Initialise dF / dmul
    double df_dmul = 0;
    // Initialise updated stress
    VectorD6x1 updated_stress = trial_stress;
    // Iteration for consistency parameter
    while (std::fabs(f_function_) > Ftolerance && counter_f < itrstep) {
      // Get back the m_theta of trial_stress
      m_theta_ = m_theta_trial;
      // Compute dF / dmul
      this->compute_df_dmul(&df_dmul);
      // Update consistency parameter
      delta_phi_ -= f_function_ / df_dmul;
      // Initialise G and dG / dpc
      double g_function = 0;
      double dg_dpc = 0;
      // Compute G and dG / dpc
      this->compute_dg_dpc(pc_n, p_trial, &g_function, &dg_dpc);
      // Subiteraction for preconsolidation pressure
      while (std::fabs(g_function) > Gtolerance && counter_g < substep) {
        // Update preconsolidation pressure
        pc_ -= g_function / dg_dpc;
        // Update G and dG / dpc
        this->compute_dg_dpc(pc_n, p_trial, &g_function, &dg_dpc);
        // Counter subiteration step
        ++counter_g;
      }
      // Update mean pressure p
      p_ = (p_trial + bulk_modulus_ * delta_phi_ * pc_) /
          (1 + 2 * bulk_modulus_ * delta_phi_);
      // Update deviatoric stress q
      // Equation(3.10b)
      q_ = q_trial / (1 + 6 * shear_modulus_ * delta_phi_ /
                           std::pow(m_theta_, 2));
      // Compute incremental plastic volumetic strain
      // Equation(2.8)
      dpvstrain_ = delta_phi_ * (2 * p_ - pc_);
      // Compute plastic deviatoric strain
      dpdstrain_ = delta_phi_ * (std::sqrt(6) * q_ / std::pow(m_theta_, 2));
      // Update yield function
      yield_type = this->compute_yield_state();
      // Counter iteration step
      ++counter_f;
    }
    // Update plastic strain
    pvstrain_ += dpvstrain_;
    pdstrain_ += dpdstrain_;
    // Update stress
    updated_stress = q_ * n_trial;
    for (int i = 0; i < 3; ++i) updated_stress(i) -= p_;
    // Update void_ratio
    //void_ratio_ +=
    //  ((dstrain(0) + dstrain(1) + dstrain(2)) * (1 + e0_));

    stress = updated_stress.transpose();
    //std::cout << "stress " << stress << "\n";
    epds = pdstrain_;
    //void_ratio = void_ratio_;
    pc = pc_;
    //std::cout << "pc " << pc << "\n";

  }
}

