mpm::Node::Node(std::string& input_line, const unsigned& id) {
  nid_ = id;
  std::istringstream input(input_line);
  double third;
  input >> ncoord_(0) >> ncoord_(1);
  input >> third;

  nphase_ = 1;
  penalty_factor_ = 0.0;

  nsolid_mass_     = 0.;
  nwater_mass_     = 0.;
  nmixture_mass_   = 0.;
  nintrmd_mass_ = 0.;

  nsolid_momentum_ = Eigen::Matrix<double,1,dim>::Zero();
  nwater_momentum_ = Eigen::Matrix<double,1,dim>::Zero();
  npressure_ = 0.;
  npressure_increment_ = 0;

  nsolid_initial_velocity_ = Eigen::Matrix<double,1,dim>::Zero();
  nwater_initial_velocity_ = Eigen::Matrix<double,1,dim>::Zero();
  nsolid_final_velocity_ = Eigen::Matrix<double,1,dim>::Zero();
  nwater_final_velocity_ = Eigen::Matrix<double,1,dim>::Zero();

  nsolid_int_acceleration_ = Eigen::Matrix<double,1,dim>::Zero();
  nsolid_int_velocity_ = Eigen::Matrix<double,1,dim>::Zero();
  nwater_int_velocity_ = Eigen::Matrix<double,1,dim>::Zero();
  
  nsolid_final_acceleration_ = Eigen::Matrix<double,1,dim>::Zero();
  nwater_final_acceleration_ = Eigen::Matrix<double,1,dim>::Zero();  

 
  nmixture_body_force_ = Eigen::Matrix<double,1,dim>::Zero();
  nwater_body_force_ = Eigen::Matrix<double,1,dim>::Zero();
  nmixture_trac_force_ = Eigen::Matrix<double,1,dim>::Zero();
  nwater_trac_force_ = Eigen::Matrix<double,1,dim>::Zero();
  nmixture_int_force_ = Eigen::Matrix<double,1,dim>::Zero();
  nwater_int_force_ = Eigen::Matrix<double,1,dim>::Zero();

  KS2_force_ = Eigen::Matrix<double,1,dim>::Zero();
  KS2_forces_= Eigen::Matrix<double,numMats,dim>::Zero();

  nwater_damping_ = Eigen::Matrix<double,1,dim>::Zero();
  nmixture_damping_ = Eigen::Matrix<double,1,dim>::Zero();

  nsolid_constraints_.clear();
  nwater_constraints_.clear();
  nsolid_accel_constraints_.clear();
  nwater_accel_constraints_.clear();
  std::get<0>(npressure_constraints_) = 0;
  std::get<1>(npressure_constraints_) = 0.;
  
  undrained_status_ = false;
  moving_undrained_ = false;
  std::get<0>(undrained_node_) = 0;
  std::get<1>(undrained_node_) = 0;

  assigned_ = 0; // needs to improved
  free_pressure_ = false;
}


void mpm::Node::initialise_node() {
  nphase_ = 1;
  material_ids_.clear();
  distances0_.clear();
  distances1_.clear();

  penalty_factor_ = 0;

  nsolid_mass_     = 0.;
  nwater_mass_     = 0.;
  nmixture_mass_ = 0.;
  nintrmd_mass_ = 0.;

  nsolid_masses_ = Eigen::Matrix<double,numMats,1>::Zero();
  nwater_masses_ = Eigen::Matrix<double,numMats,1>::Zero();
  nmixture_masses_ = Eigen::Matrix<double,numMats,1>::Zero();
  nintrmd_masses_ = Eigen::Matrix<double,numMats,1>::Zero();

  nsolid_momentum_ = Eigen::Matrix<double,1,dim>::Zero();
  nwater_momentum_ = Eigen::Matrix<double,1,dim>::Zero();
  npressure_ = 0.;
  npressure_increment_ = 0.;

  nsolid_momenta_ = Eigen::Matrix<double,numMats,dim>::Zero();
  npressures_ = Eigen::Matrix<double,numMats,1>::Zero();
  npressure_increments_ = Eigen::Matrix<double,numMats,1>::Zero();

  nsolid_initial_velocity_ = Eigen::Matrix<double,1,dim>::Zero();
  nwater_initial_velocity_ = Eigen::Matrix<double,1,dim>::Zero();
  nsolid_final_velocity_ = Eigen::Matrix<double,1,dim>::Zero();
  nwater_final_velocity_ = Eigen::Matrix<double,1,dim>::Zero();

  nsolid_initial_velocities_ = Eigen::Matrix<double,numMats,dim>::Zero();
  nsolid_final_velocities_ = Eigen::Matrix<double,numMats,dim>::Zero();

  nsolid_int_acceleration_ = Eigen::Matrix<double,1,dim>::Zero();
  nsolid_int_velocity_ = Eigen::Matrix<double,1,dim>::Zero();
  nwater_int_velocity_ = Eigen::Matrix<double,1,dim>::Zero();

  nsolid_int_accelerations_ = Eigen::Matrix<double,numMats,dim>::Zero();
  nsolid_int_velocities_ = Eigen::Matrix<double,numMats,dim>::Zero();
  nwater_int_velocities_ = Eigen::Matrix<double,numMats,dim>::Zero();
  
  nsolid_final_acceleration_ = Eigen::Matrix<double,1,dim>::Zero(); 
  nsolid_final_accelerations_ = Eigen::Matrix<double,numMats,dim>::Zero(); 
 
  nmixture_body_force_ = Eigen::Matrix<double,1,dim>::Zero();
  nwater_body_force_ = Eigen::Matrix<double,1,dim>::Zero();
  nmixture_trac_force_ = Eigen::Matrix<double,1,dim>::Zero();
  nwater_trac_force_ = Eigen::Matrix<double,1,dim>::Zero();
  nmixture_int_force_ = Eigen::Matrix<double,1,dim>::Zero();
  nwater_int_force_ = Eigen::Matrix<double,1,dim>::Zero();

  nmixture_body_forces_ = Eigen::Matrix<double,numMats,dim>::Zero();
  nwater_body_forces_ = Eigen::Matrix<double,numMats,dim>::Zero();
  nmixture_trac_forces_ = Eigen::Matrix<double,numMats,dim>::Zero();
  nwater_trac_forces_ = Eigen::Matrix<double,numMats,dim>::Zero();
  nmixture_int_forces_ = Eigen::Matrix<double,numMats,dim>::Zero();
  nwater_int_forces_ = Eigen::Matrix<double,numMats,dim>::Zero();

  KS2_force_ = Eigen::Matrix<double,1,dim>::Zero();
  KS2_forces_= Eigen::Matrix<double,numMats,dim>::Zero();

  nwater_damping_ = Eigen::Matrix<double,1,dim>::Zero();
  nmixture_damping_ = Eigen::Matrix<double,1,dim>::Zero();
    
  // if there is a moving free surface, the pressure boundary condition
  // is updated at each time step
  if(mpm::misc::freeSurface) {
    if(free_pressure_) {
      std::get<0>(npressure_constraints_) = 0;
      std::get<1>(npressure_constraints_) = 0.;
    }
  }
  if(moving_undrained_) {
    undrained_status_ = false;
    moving_undrained_ = false;
  }
  assigned_ = 0;
  free_pressure_ = false;
  pressure_status_ = true;
  contact_status_ = false;
}

void mpm::Node::compute_nodal_initial_velocity() {
// Compute COM nodal initial velocity
  if (std::fabs(nsolid_mass_) > 1.0E-18)
    nsolid_initial_velocity_ = nsolid_momentum_ / nsolid_mass_;

  for (unsigned i = 0; i < dim; i++) 
    this->check_double_precision(nsolid_initial_velocity_(i));
  
  // apply velocity constraints
  if(nsolid_constraints_.size()) {
    for (unsigned j = 0; j < nsolid_constraints_.size(); j++) {
      unsigned direction = std::get<0>(nsolid_constraints_.at(j));
      nsolid_initial_velocity_(direction) = std::get<1>(nsolid_constraints_.at(j));
    }
  }
  
// Compute nodal initial velocities for each material body
  for (unsigned i=0;  i< numMats; i++) {
    if (std::fabs(nsolid_masses_(i)) > 1.0E-18)
      nsolid_initial_velocities_.row(i) = nsolid_momenta_.row(i) / nsolid_masses_(i);

    for (unsigned j = 0; j < dim; j++)
      this->check_double_precision(nsolid_initial_velocities_(j,i));
      
      // apply velocity constraints
    if(nsolid_constraints_.size()) {
      for (unsigned j = 0; j < nsolid_constraints_.size(); j++) {
          unsigned direction = std::get<0>(nsolid_constraints_.at(j));
          nsolid_initial_velocities_(i,direction) = std::get<1>(nsolid_constraints_.at(j));
      }
    }
  }
  // Apply Contact Mechanics at the interface
  if (contact_status_ == true) {
    for (unsigned mitr = 0; mitr < numMats; mitr ++) {
      Eigen::Matrix<double, 1, dim> relative_velocity = nsolid_initial_velocity_ - nsolid_initial_velocities_.row(mitr);
      double velocity_normal = relative_velocity.dot(normal_unit_vectors_.row(mitr));
      if (std::fabs(velocity_normal) < 1E-18)
        velocity_normal = 0;
      Eigen::Matrix<double,1,dim> normal_correction = velocity_normal*normal_unit_vectors_.row(mitr);
      nsolid_initial_velocities_.row(mitr) = nsolid_initial_velocities_.row(mitr) + normal_correction;
      for (unsigned i=0; i<dim; i++) {
        if (std::fabs(nsolid_initial_velocities_(0,i)) < 1E-18)
          nsolid_initial_velocities_(0,i) = 0.0;
      }
    }
  }
}

void mpm::Node::adjust_nodal_initial_velocity () {
  if (contact_status_ == true) {
      Eigen::Matrix<double, 1, dim> relative_velocity = nsolid_initial_velocity_ - nsolid_initial_velocities_.row(0);
    double velocity_normal = relative_velocity.dot(normal_unit_vectors_.row(0));
    if (std::fabs(velocity_normal) < 1E-18)
      velocity_normal = 0;
    Eigen::Matrix<double,1,dim> normal_correction = velocity_normal*normal_unit_vectors_.row(0);
    nsolid_initial_velocities_.row(0) = nsolid_initial_velocities_.row(0) + normal_correction;
    for (unsigned i=0; i<dim; i++) {
      if (std::fabs(nsolid_initial_velocities_(0,i)) < 1E-18)
        nsolid_initial_velocities_(0,i) = 0.0;
    }
  }
}

void mpm::Node::matrix_test() {
  if (nid_ == 3 || nid_ == 303 || nid_ == 309 || nid_ == 315) {
    std::cout << "node_id: " << nid_ << " \n";
    std::cout << "nsolid_masses_: \n" << nsolid_masses_ << " \n";
    std::cout << "nsolid_mass_: \n" << nsolid_mass_ << " \n";
    std::cout << "nwater_masses_: \n" << nwater_masses_ << " \n";
    std::cout << "nwater_mass_: \n" << nwater_mass_ << "\n";
    std::cout << "nintrmd_masses_: \n" << nintrmd_masses_ << "\n";
    std::cout << "mixture_masses_: \n" << nmixture_masses_ << " \n";
    std::cout << "mixture_mass_: \n" << nmixture_mass_ << " \n\n";
    std::cout << "nsolid_momentum: \n" << nsolid_momentum_ << "\n";
    std::cout << "nsolid_momenta: \n" << nsolid_momenta_ << "\n";
    //std::cout << "nmixture_body_forces_: \n" << nmixture_body_forces_ << " \n";
    //std::cout << "nmixture_body_force_: \n" << nmixture_body_force_ << " \n";
    //std::cout << "nwater_body_forces_: \n" << nwater_body_forces_ << " \n";
    //std::cout << "nwater_body_force_: \n" << nwater_body_force_ << " \n";
    std::cout << "nmixture_trac_forces_: \n" << nmixture_trac_forces_ << " \n";
    std::cout << "nmixture_trac_force_: \n" << nmixture_trac_force_ << " \n";
    std::cout << "nmixture_int_forces_: \n" << nmixture_int_forces_ << " \n";
    std::cout << "nmixture_int_force_: \n" << nmixture_int_force_ << " \n";
    std::cout << "nwater_int_forces_: \n" << nwater_int_forces_ << " \n";
    std::cout << "nwater_int_force_: \n" << nwater_int_force_ << " \n\n";
    std::cout << "KS2 Force_: \n" << KS2_force_ << " \n\n";
    std::cout << "nsolid_initial_velocities_: \n" << nsolid_initial_velocities_ << " \n";
    std::cout << "nsolid_intitial_velocity_: \n" << nsolid_initial_velocity_ << " \n";  
    std::cout << "nsolid_int_accelerations_: \n" << nsolid_int_accelerations_ << " \n";
    std::cout << "nsolid_int_acceleration_: \n" << nsolid_int_acceleration_ << " \n";  
    std::cout << "nsolid_int_velocities_: \n" << nsolid_int_velocities_ << " \n";
    std::cout << "nsolid_int_velocity_: \n" << nsolid_int_velocity_ << " \n";    
    std::cout << "nwater_int_velocities_: \n" << nwater_int_velocities_ << " \n";
    std::cout << "nwater_int_velocity_: \n" << nwater_int_velocity_ << " \n";
    std::cout << "nsolid_final_accelerations_: \n" << nsolid_final_accelerations_ << " \n";
    std::cout << "nsolid_final_acceleration_: \n" << nsolid_final_acceleration_ << " \n"; 
    std::cout << "unit_normal_vector: \n" << normal_unit_vectors_ <<"\n";
    std::cout << "penalty_factor: \n" << penalty_factor_ << "\n";
    std::cout << "nsolid_final_velocities_: \n" << nsolid_final_velocities_ << " \n";
    std::cout << "nsolid_final_velocity_: \n" << nsolid_final_velocity_ << " \n";
    std::cout << "contact_status_: \n" << contact_status_ << "\n\n";  
    //std::cout << "relative_velocities: \n" << nsolid_relative_velocities_ << "\n\n";
    //std::cout << "ncoord_: \n" << ncoord_ << "\n\n";
    std::cout << "=====================================================" << "\n\n";
  }
}

void mpm::Node::write_para_data(std::ostream& oFile) {
    oFile << "=====================================================" << "\n";
    oFile << "node_id: " << nid_ << " \n";
    //oFile << "Interface: " << material_ids_.size() << "\n\n";
    //oFile << "intermediate_masses_: \n" << nintrmd_masses_ << " \n";
    oFile << "water_masses_: \n" << nwater_masses_ << " \n";
    oFile << "mixture_masses_: \n" << nmixture_masses_ << " \n";
    oFile << "mixture_body_forces_: \n" << nmixture_body_forces_ << " \n";
    oFile << "nwater_body_forces_: \n" << nwater_body_forces_ << "\n";
    oFile << "nmixture_trac_forces_: \n" << nmixture_trac_forces_ << " \n";
    oFile << "nmixture_int_forces_: \n" << nmixture_int_forces_ << " \n";
    //oFile << "nwater_trac_forces_: \n" << nwater_trac_forces_ << " \n";
    oFile << "nwater_int_forces_: \n" << nwater_int_forces_ << " \n";
    oFile << "nwater_int_velocities_: \n" << nwater_int_velocities_ << " \n";
    oFile << "unit_normal_vector: \n" << normal_unit_vectors_ <<"\n";
    oFile << "npressures_: \n" << npressures_ << " \n";
    oFile << "pressure_increment: \n" << npressure_increment_ << " \n";
    oFile << "pressure_increments: \n" << npressure_increments_ << " \n";
    oFile << "KS2 Force_: \n" << KS2_force_ << " \n";
    oFile << "KS2 Forces_: \n" << KS2_forces_ << " \n\n";
    oFile << "penalty factor: \n" << penalty_factor_ << "\n";
    oFile << "Contact Status: \n" << contact_status_ << "\n";
    oFile << "Free Surface: \n" << free_pressure_ << "\n";
    oFile << "nsolid_initial_velocity: \n" << nsolid_initial_velocity_ << "\n";
    oFile << "nsolid_initial_velocities_: \n" << nsolid_initial_velocities_ << " \n";
    oFile << "nsolid_int_acceleration_: \n" << nsolid_int_acceleration_ << " \n";
    oFile << "nsolid_int_accelerations_: \n" << nsolid_int_accelerations_ << " \n";
    oFile << "nsolid_int_velocity_: \n" << nsolid_int_velocity_ << " \n";
    oFile << "nsolid_int_velocities_: \n" << nsolid_int_velocities_ << " \n"; 
    oFile << "nsolid_final_acceleration_: \n" << nsolid_final_acceleration_ << " \n";
    oFile << "nsolid_final_accelerations_: \n" << nsolid_final_accelerations_ << " \n";
    oFile << "nsolid_final_velocity: \n" << nsolid_final_velocity_ << " \n";
    oFile << "nsolid_final_velocities: \n" << nsolid_final_velocities_ << " \n\n";
    oFile << "=====================================================" << "\n\n";
}

void mpm::Node::identify_contact_nodes() {
  if (material_ids_.size() == 2){
    //std::vector<double> d0;
    //std::vector<double> d1;
    //std::copy(normal_distances0_.begin(), normal_distances0_.end(), std::back_inserter(d0));
    //std::copy(normal_distances1_.begin(), normal_distances1_.end(), std::back_inserter(d1));
    double nd_pipe;
    double td_pipe;
    double nd_soil;
    double td_soil;
    double d_alpha;
    double d_beta;
    double d_RL;
    if (distances0_.size() && distances1_.size()) {
      for (auto i=0; i < distances0_.size(); i++) {
        nd_soil = std::get<0>(distances0_.at(i));
        td_soil = std::get<1>(distances0_.at(i));
        for (auto j=0; j< distances1_.size(); j++) {
          nd_pipe = std::get<0>(distances1_.at(j));
          td_pipe = std::get<1>(distances1_.at(j));
          d_RL = nd_soil - nd_pipe;
          double td_RL = std::fabs(td_soil - td_pipe);
          if (d_RL < 1E-04 && td_RL < 0.01) {
            contact_status_ = true;
            penalty_factor_ = 1.0;
            //d_alpha = nd_pipe;
            //d_beta = nd_soil;
            //std::cout << ncoord_ << "\n";
            //std::cout << "distance to pipe: " << d_alpha << "\n"; 
            //std::cout << "distance to soil: " << d_beta << "\n"; 
            //std::cout << "distance: " << d_RL << "\n";
            //std::cout << "t_distance: " << td_RL << "\n";
          }    
        }
      }
      //double d_alpha = *std::max_element(d1.begin(), d1.end());
      //double d_beta = *std::min_element(std::begin(d0), std::end(d0));
      //double d_RL = d_beta-d_alpha;
      //if(contact_status_ == true) {
        //std::cout << "contact node: " << nid_ << "\n";
        //penalty_factor_ = pow((-d_RL)/4E-05,3);
        //std::cout<< "Pf: " << penalty_factor_ << "\n";
        //penalty_factor_ = 1;
        //if (penalty_factor_ > 1)
        //  penalty_factor_ = 1.0;
      //}
    }
  }
}

void mpm::Node::compute_nodal_pressure() {   
  for (unsigned i=0; i<numMats; i++){
    if (std::fabs(nwater_masses_(i)) > 1.E-18)
      npressures_(i) = npressures_(i) / nwater_masses_(i);
    this->check_double_precision(npressures_(i));
    this->apply_pressure_constraint();
  }
}

void mpm::Node::compute_nodal_pressure_increment() {   
  for (unsigned i=0; i<numMats; i++){
    if (std::fabs(nwater_masses_(i)) > 1.E-18)
      npressure_increments_(i) = npressure_increments_(i) / nwater_masses_(i);
    this->check_double_precision(npressure_increments_(i));
    this->apply_pressure_constraint();
  }
}

void mpm::Node::compute_intermediate_solid_acceleration(const double &dt) {
  // Compute nodal int accelerations and velocities for each material body
  for (unsigned i=0; i < numMats; i++) {
    if (std::fabs(nmixture_masses_(i)) > 1E-18) {
      Eigen::Matrix<double, 1, dim> force = Eigen::Matrix<double, 1, dim>::Zero();
      force = nmixture_body_forces_.row(i) + nmixture_int_forces_.row(i) + nmixture_trac_forces_.row(i);
      for (unsigned j=0; j < dim; j++) {
        if (std::fabs(force(j)) < 1E-15) {
          force(j) = 0;
        }
      }
      nsolid_int_accelerations_.row(i) = force / nmixture_masses_(i);
      nsolid_int_velocities_.row(i) = nsolid_initial_velocities_.row(i) + nsolid_int_accelerations_.row(i)*dt;
    }  
  }
  nsolid_int_accelerations_(1,0) = 0;
  nsolid_int_velocities_(1,0) = 0;

  // Compute COM nodal int accelerations and velocities
  if (std::fabs(nmixture_mass_ > 1E-18)) {
      Eigen::Matrix<double, 1, dim> force = Eigen::Matrix<double, 1, dim>::Zero();
      force = nmixture_body_force_ + nmixture_int_force_ + nmixture_trac_force_;
      for (unsigned j=0; j < dim; j++) {
        if (std::fabs(force(j)) < 1E-15) {
          force(j) = 0;
        }
      }
      nsolid_int_acceleration_ = force / nmixture_mass_;
      nsolid_int_velocity_ = nsolid_initial_velocity_ + nsolid_int_acceleration_*dt;
  }

  // Apply contact mechancis at interface
  if (contact_status_ == true) {
    nsolid_int_velocity_(0) = 0;
    nsolid_int_acceleration_(0) = 0;
    for (unsigned mitr = 0; mitr < numMats; mitr ++) {
      // Adjust velocities
      Eigen::Matrix<double, 1, dim> relative_velocity = nsolid_int_velocity_ - nsolid_int_velocities_.row(mitr);
      double velocity_normal = relative_velocity.dot(normal_unit_vectors_.row(mitr));
      if (std::fabs(velocity_normal) < 1E-18)
        velocity_normal = 0;
      Eigen::Matrix<double,1,dim> normal_correction = velocity_normal*normal_unit_vectors_.row(mitr);
      nsolid_int_velocities_.row(mitr) = nsolid_int_velocities_.row(mitr) + normal_correction;
      for (unsigned i=0; i<dim; i++) {
        if (std::fabs(nsolid_int_velocities_(0,i)) < 1E-18)
          nsolid_int_velocities_(0,i) = 0.0;
      }
      // Adjust accelerations
      Eigen::Matrix<double, 1, dim> relative_acceleration = nsolid_int_acceleration_ - nsolid_int_accelerations_.row(mitr);
      double acceleration_normal = relative_acceleration.dot(normal_unit_vectors_.row(mitr));
      if (std::fabs(acceleration_normal) < 1E-18)
        acceleration_normal = 0;
      Eigen::Matrix<double,1,dim> normal_acc_correction = acceleration_normal*normal_unit_vectors_.row(mitr);
      nsolid_int_accelerations_.row(mitr) = nsolid_int_accelerations_.row(mitr) + normal_acc_correction;
      for (unsigned i=0; i<2; i++) {
        if (std::fabs(nsolid_int_accelerations_(mitr,i)) < 1E-18) {
          nsolid_int_accelerations_(mitr,i) = 0;
        }
      }
    }
  }
  // apply velocity constraints 
  for(auto mitr = 0; mitr < numMats; mitr++){
    if(nsolid_constraints_.size()) {
      for (unsigned j = 0; j < nsolid_constraints_.size(); j++) {
        unsigned direction = std::get<0>(nsolid_constraints_.at(j));
        nsolid_int_velocities_(mitr,direction) = std::get<1>(nsolid_constraints_.at(j));
        nsolid_int_accelerations_(mitr,direction) = 0;
      }
    }
  }
}

void mpm::Node::compute_final_solid_acceleration(const double &dt) {
  // Compute nodal final accelerations and velocities for each material body
  for (unsigned i=0; i < numMats; i++) {
    if (std::fabs(nmixture_masses_(i)) > 1E-18) {
      Eigen::Matrix<double, 1, dim> force = Eigen::Matrix<double, 1, dim>::Zero();
      force = nmixture_body_forces_.row(i) + nmixture_int_forces_.row(i) + nmixture_trac_forces_.row(i) - KS2_forces_.row(i);
      for (unsigned j=0; j < dim; j++) {
        if (std::fabs(force(j)) < 1E-15) {
          force(j) = 0;
        }
      }
      nsolid_final_accelerations_.row(i) = force / nmixture_masses_(i);
      nsolid_final_velocities_.row(i) = nsolid_initial_velocities_.row(i) + nsolid_final_accelerations_.row(i)*dt;
    }  
  }
  nsolid_final_velocities_(1,0) = 0;

  // Compute COM nodal final accelerations and velocities
  if (std::fabs(nmixture_mass_ > 1E-18)) {
      Eigen::Matrix<double, 1, dim> force = Eigen::Matrix<double, 1, dim>::Zero();
      force = nmixture_body_force_ + nmixture_int_force_ + nmixture_trac_force_ - KS2_force_;
      for (unsigned j=0; j < dim; j++) {
        if (std::fabs(force(j)) < 1E-15) {
          force(j) = 0;
        }
      }
      nsolid_final_acceleration_ = force / nmixture_mass_;
      nsolid_final_velocity_ = nsolid_initial_velocity_ + nsolid_final_acceleration_*dt;
  }

  // Apply contact mechancis at interface
  if (contact_status_ == true) {
    nsolid_final_velocity_(0) = 0;
    nsolid_final_acceleration_(0) = 0;
    for (unsigned mitr = 0; mitr < numMats; mitr ++) {
      // Adjust velocities
      Eigen::Matrix<double, 1, dim> relative_velocity = nsolid_final_velocity_ - nsolid_final_velocities_.row(mitr);
      double velocity_normal = relative_velocity.dot(normal_unit_vectors_.row(mitr));
      if (std::fabs(velocity_normal) < 1E-18)
        velocity_normal = 0;
      Eigen::Matrix<double,1,dim> normal_correction = velocity_normal*normal_unit_vectors_.row(mitr);
      nsolid_final_velocities_.row(mitr) = nsolid_final_velocities_.row(mitr) + normal_correction;
      for (unsigned i=0; i<dim; i++) {
        if (std::fabs(nsolid_final_velocities_(0,i)) < 1E-18)
          nsolid_final_velocities_(0,i) = 0.0;
      }
      // Adjust accelerations
      Eigen::Matrix<double, 1, dim> relative_acceleration = nsolid_final_acceleration_ - nsolid_final_accelerations_.row(mitr);
      double acceleration_normal = relative_acceleration.dot(normal_unit_vectors_.row(mitr));
      if (std::fabs(acceleration_normal) < 1E-18)
        acceleration_normal = 0;
      Eigen::Matrix<double,1,dim> normal_acc_correction = acceleration_normal*normal_unit_vectors_.row(mitr);
      nsolid_final_accelerations_.row(mitr) = nsolid_final_accelerations_.row(mitr) + normal_acc_correction;
      for (unsigned i=0; i<2; i++) {
        if (std::fabs(nsolid_final_accelerations_(mitr,i)) < 1E-18) {
          nsolid_final_accelerations_(mitr,i) = 0;
        }
      }
    }
  }

  // apply velocity constraints 
  for(auto mitr = 0; mitr < numMats; mitr++){
    if(nsolid_constraints_.size()) {
      for (unsigned j = 0; j < nsolid_constraints_.size(); j++) {
        unsigned direction = std::get<0>(nsolid_constraints_.at(j));
        nsolid_final_velocities_(mitr,direction) = std::get<1>(nsolid_constraints_.at(j));
        nsolid_final_accelerations_(mitr,direction) = 0;
      }
    }
  }
}

void mpm::Node::compute_intermediate_water_velocity() {
  if (std::fabs(nwater_masses_(0)) > 1.E-18) {
    Eigen::Matrix<double, 1, dim> force = Eigen::Matrix<double, 1, dim>::Zero();
    force = (nwater_body_forces_.row(0) + nwater_int_forces_.row(0) + nwater_trac_forces_.row(0) 
    - (nintrmd_masses_.row(0) * nsolid_int_accelerations_.row(0)));
    for (unsigned i=0; i<dim; i++) {
      if (std::fabs(force(i)) < 1E-12) {
        force(i) = 0;
      }
    }
    nwater_int_velocities_.row(0) = force / nwater_masses_(0);
    for (unsigned i=0; i<dim; i++) {
      if (std::fabs(nwater_int_velocities_(0,i)) < 1E-18)
        nwater_int_velocities_(0,i) = 0.0;
    }
    if (contact_status_ == true) {
      nwater_int_velocities_(0,1) = 0;
      double tolerance = 1.E-18;
      Eigen::Matrix<double,1,dim> normal_unit_vector = normal_unit_vectors_.row(1);
      Eigen::Matrix<double,1,dim> relative_velocity = -nwater_int_velocities_.row(0);
      double velocity_normal = relative_velocity.dot(normal_unit_vector);
      if (std::fabs(velocity_normal) < tolerance)
        velocity_normal = 0;
      if (velocity_normal > 0){
        Eigen::Matrix<double,1,dim> corrections = Eigen::Matrix<double,1,dim>::Zero();
        Eigen::Matrix<double,1,dim> normal_correction = velocity_normal * normal_unit_vector;
      // Update the velocity with the computed corrections
        corrections = normal_correction;
        nwater_int_velocities_.row(0) = nwater_int_velocities_.row(0) + corrections;
        for (auto j=0; j<dim; j++){
          if (std::abs(nwater_int_velocities_(j)) < tolerance)
              nwater_int_velocities_(0,j) = 0.0;
        }
      }
    }
    if(nwater_constraints_.size()) {
      for (unsigned i = 0; i < nwater_constraints_.size(); i++) {
        unsigned direction = std::get<0>(nwater_constraints_.at(i));
        nwater_int_velocities_(0,direction) = std::get<1>(nwater_constraints_.at(i));
      }
    }
    if(undrained_status_) {
      unsigned direction = std::get<0>(undrained_node_);
      nwater_int_velocities_(0,direction) = 0.;
    }
  }
}

void mpm::Node::compute_multimaterial_relative_velocities(){
  for (auto mitr = 0; mitr < numMats; mitr++){
    nsolid_relative_velocities_.row(mitr) = nsolid_final_velocities_.row(mitr) - nsolid_final_velocity_;
  }
}

//void mpm::Node::compute_multimaterial_normal_unit_vectors() {
//  Eigen::Matrix<double,1,dim> largest_domain_gradient = Eigen::Matrix<double,1,dim>::Zero();
//  double max_magnitude = 0;
//  unsigned mat_id_largest = 0;
//  for (unsigned i = 0; i < numMats; i++){
//    Eigen::Matrix<double,1,dim> current_domain_gradient = domain_gradients_.row(i);
//    if (current_domain_gradient.norm() >= max_magnitude){
//      max_magnitude = current_domain_gradient.norm();
//      largest_domain_gradient = current_domain_gradient;
//      mat_id_largest = i;
//    }
//  }
//  for (unsigned i = 0; i < numMats; i++) {
//    Eigen::Matrix<double,1,dim> normal_unit_vector = Eigen::Matrix<double,1,dim>::Zero();
//    if (largest_domain_gradient.norm() > 1.E-18)
//      normal_unit_vector = largest_domain_gradient.normalized();
//    if (mat_id_largest != i)
//      normal_unit_vector = -1 * normal_unit_vector;
    //if (nid_ == 17 || nid_ == 18) {
    // 1d con
//    if(nid_ == 306 || nid_ == 307 || nid_ == 308 || nid_ == 309) {  
//      if (i == 0) {
//        normal_unit_vector(0) = 0;
//        normal_unit_vector(1) = 1;
//      }
//      else if (i == 1) {
//        normal_unit_vector(0) = 0;
//        normal_unit_vector(1) = -1;
//      }
//    }
//    normal_unit_vectors_.row(i) = normal_unit_vector;
//  }
//}

void mpm::Node::compute_multimaterial_normal_unit_vectors(double coord) {
  double x = ncoord_(0);
  double y = ncoord_(1)-coord;
  double r = sqrt(pow(x,2)+pow(y,2)); 
  Eigen::Matrix<double,1,dim> normal_unit_vector = Eigen::Matrix<double,1,dim>::Zero();
  normal_unit_vector(0) = x/r;
  normal_unit_vector(1) = y/r;
  //normal_unit_vector(0) = 0;
  //normal_unit_vector(1) = -1;
  normal_unit_vectors_.row(1) = normal_unit_vector;
  normal_unit_vectors_.row(0) = -normal_unit_vector;
}

void mpm::Node::apply_contact_mechanics(const double& dt){
  if (contact_status_ == true) {
    double tolerance = 1.E-15;
    unsigned i = 0;
    Eigen::Matrix<double,1,dim> normal_unit_vector = normal_unit_vectors_.row(i);
    Eigen::Matrix<double,1,dim> relative_velocity = nsolid_relative_velocities_.row(i);
    double velocity_normal = relative_velocity.dot(normal_unit_vector);
    if (std::fabs(velocity_normal) < tolerance)
      velocity_normal = 0;
    Eigen::Matrix<double,1,dim> corrected_velocity = nsolid_final_velocities_.row(i);
    if (std::fabs(velocity_normal) > 0){
      Eigen::Matrix<double,1,dim> corrections = Eigen::Matrix<double,1,dim>::Zero();
      double cross_product =
              relative_velocity(0) * normal_unit_vector(1) -
              relative_velocity(1) * normal_unit_vector(0);
      double mu = 0;
      // Compute the normal and tangential corrections
      Eigen::Matrix<double,1,dim> normal_correction = -velocity_normal * normal_unit_vector;
      Eigen::Matrix<double,1,dim> tangent_correction = Eigen::Matrix<double,1,dim>::Zero();
      //tangent_correction(0) = normal_unit_vector(1) * cross_product;
      //tangent_correction(1) = -normal_unit_vector(0) * cross_product;
      //tangent_correction = -mu * velocity_normal / std::abs(cross_product) *
      //                       tangent_correction;

      // Update the velocity with the computed corrections
      corrections = (normal_correction + tangent_correction);
      corrected_velocity = corrected_velocity + corrections;
      for (auto j=0; j<dim; j++){
        if (std::abs(corrected_velocity(j)) < tolerance)
            corrected_velocity(j) = 0.0;
      }
      nsolid_final_velocities_.row(i) = corrected_velocity;
      Eigen::Matrix<double,1,dim> corrected_acc = corrections / dt;
      nsolid_final_accelerations_.row(i) = nsolid_final_accelerations_.row(i) + corrected_acc;
    }
    //}
  }
}

void mpm::Node::check_double_precision(double& value) {
  if (std::fabs(value) < 1.E-18)
    value = 0.;
}

void mpm::Node::apply_solid_velocity_constraint(Eigen::Matrix<double,1,dim> velocity) {
    if(nsolid_constraints_.size()) {
        for (unsigned i = 0; i < nsolid_constraints_.size(); i++) {
            unsigned direction = std::get<0>(nsolid_constraints_.at(i));
            velocity(direction) = std::get<1>(nsolid_constraints_.at(i));
        }
    }
}

void mpm::Node::apply_water_velocity_constraint(Eigen::Matrix<double,1,dim> velocity) {
  if(nwater_constraints_.size()) {
    for (unsigned i = 0; i < nwater_constraints_.size(); i++) {
      unsigned direction = std::get<0>(nwater_constraints_.at(i));
      velocity(direction) = std::get<1>(nwater_constraints_.at(i));
    }
  }
}

void mpm::Node::apply_solid_acceleration_constraint(Eigen::Matrix<double,1,dim> acceleration) {
  if(nsolid_accel_constraints_.size()) {
    for (unsigned i = 0; i < nsolid_accel_constraints_.size(); i++) {
      unsigned direction = std::get<0>(nsolid_accel_constraints_.at(i));
      acceleration(direction) = 0.0;
    }
  }
}

void mpm::Node::apply_water_acceleration_constraint(Eigen::Matrix<double,1,dim> &acceleration) {
    if(nwater_accel_constraints_.size()) {
        for (unsigned i = 0; i < nwater_accel_constraints_.size(); i++) {
            unsigned direction = std::get<0>(nwater_accel_constraints_.at(i));
            acceleration(direction) = 0.0;
        }
    }
}


void mpm::Node::apply_pressure_constraint() {
    if (std::get<0>(npressure_constraints_)) {
      npressure_ = std::get<1>(npressure_constraints_);
      for (unsigned i=0; i<numMats; i++) {
        npressures_(i) = std::get<1>(npressure_constraints_);
      }  
    }
}

void mpm::Node::apply_undrained_conditions(Eigen::Matrix<double,1,dim> water_relative_velocity) {
  if(undrained_status_) {
    unsigned direction = std::get<0>(undrained_node_);
    water_relative_velocity(direction) = 0.;
  }
}

