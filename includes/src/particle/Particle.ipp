mpm::Particle::Particle(const unsigned& id, const unsigned& material_id, const Eigen::Matrix<double, 1, dim>& spacing,const double& vol, const unsigned& nphase) {
  id_ = id;
  mat_id_ = material_id;
  spacing_ = spacing;
  volume_ = vol;
  nphase_ = nphase;
  spacing_(0) = spacing_(1)=sqrt(volume_);
  //particle_radius_ = sqrt(volume_/3.14159);
  particle_radius_ = sqrt(volume_/4);
  gravity_ = Eigen::Matrix<double, 1, dim>::Zero();
  if (mpm::misc::gravity)
    gravity_(dim -1) = -9.81;

  solid_traction_.clear();
  water_traction_.clear();

  solid_mass_ = 0.;
  water_mass_ = 0.;

  displacement_ = Eigen::Matrix<double, 1, dim>::Zero();

  solid_velocity_ = Eigen::Matrix<double, 1, dim>::Zero();
  water_velocity_ = Eigen::Matrix<double, 1, dim>::Zero();

  solid_acceleration_ = Eigen::Matrix<double, 1, dim>::Zero();
  water_acceleration_ = Eigen::Matrix<double, 1, dim>::Zero();

  pore_pressure_ = 0.;
  solid_stress_ = Eigen::Matrix<double, 1, 6>::Zero();
  solid_strain_ = Eigen::Matrix<double, 1, dof>::Zero();
  solid_strain_rate_ = Eigen::Matrix<double, 1, dof>::Zero();
  solid_centre_strain_rate_ = Eigen::Matrix<double, 1, dof>::Zero();

  element_ = NULL;

  if (dim == 2) {
    m(0) = 1;
    m(1) = 1;
    m(2) = 0;
  }
  plastic_strain_ = 0.;
  for (auto i=0; i<dim; i++){
    particle_size_(i) = sqrt(volume_);
  }
}

void mpm::Particle::print_id(){
  std::cout << id_ << "\n";
}
void mpm::Particle::initialise_particle() {
    solid_strain_rate_ = Eigen::Matrix<double, 1, dof>::Zero();
    solid_centre_strain_rate_ = Eigen::Matrix<double, 1, dof>::Zero();
    du_ = Eigen::Matrix<double,1,dim>::Zero();
    pressure_increment_ = 0;
}

void mpm::Particle::set_initial_pore_pressure(double &initial_pressure) {
    pore_pressure_ = initial_pressure;
}


void mpm::Particle::set_initial_stress(const Eigen::Matrix<double, 1, 6>& initial_stress) { 
    solid_stress_ = initial_stress;
}

void mpm::Particle::assign_velocity() {
  solid_velocity_(0) = 0;
  solid_velocity_(1) = -2;
}

void mpm::Particle::set_solid_surface_traction(const unsigned &direction, const double &solidtraction) {
  std::tuple<unsigned, double> s_force = std::make_tuple(direction,solidtraction);
  solid_traction_.push_back(s_force);
}

void mpm::Particle::set_water_surface_traction(const unsigned &direction, const double& watertraction) {
  std::tuple<unsigned, double> w_force = std::make_tuple(direction,watertraction);
  water_traction_.push_back(w_force);
}


void mpm::Particle::set_material(const std::vector<mpm::material::MaterialBase*> materials) {
  material_ptrs_ = materials;
  material_ = materials.at(mat_id_);
  solid_grain_density_ = material_ -> give_density();
  water_grain_density_ = 1000.0;
  porosity_ = material_->give_porosity();
  permeability_ = material_->give_permeability();
  //pc_ = material_ -> give_pc();
  //void_ratio_ = material_ -> give_void_ratio();
  void_ratio_ = porosity_/(1 - porosity_);

  compute_mass();
}


void mpm::Particle::compute_mass() {
  solid_mass_ = (1 - porosity_) * solid_grain_density_ * volume_;
  water_mass_ = porosity_ * water_grain_density_ * volume_;
  if (water_mass_ < 1.0E-15)
    water_mass_ = 0.0;
}

void mpm::Particle::update_nodal_phase() {
  for (unsigned i = 0; i <numNodes; i++){
      unsigned nodal_phase = 2;
      nodes_(i)->assign_nodal_phase(nodal_phase);
  }
}

void mpm::Particle::append_material_id() {
  for (unsigned i = 0; i < gimpNodes; i++){
    if (shape_fun_(i) > 1E-18)
      gimp_nodes_(i) -> append_material_id(mat_id_);
  }
  element_ -> append_material_id(mat_id_);
}

void mpm::Particle::append_particle_coordinate_to_element(){
  element_ -> append_particle_coordinate(coord_);
}

void mpm::Particle::map_mass_to_nodes(bool contact) {
  double node_solid_mass = 0.;
  double node_water_mass = 0.;
  double node_mixture_mass = 0.;
  double node_intrmd_mass = 0.;
  for (unsigned i = 0; i < gimpNodes; i++) {
    node_solid_mass = solid_mass_ * shape_fun_(i);
    node_water_mass = water_mass_ * shape_fun_(i);
    node_mixture_mass = (solid_mass_ + water_mass_) * shape_fun_(i);
    node_intrmd_mass = water_mass_ * shape_fun_(i) * permeability_ / (porosity_ * 9.81);
    gimp_nodes_(i)->assign_nodal_masses(node_solid_mass, node_water_mass, node_mixture_mass);
    gimp_nodes_(i)->assign_nodal_intrmd_mass(node_intrmd_mass);
    if (contact) {
      gimp_nodes_(i)->assign_multimaterial_nodal_masses(mat_id_, node_solid_mass, node_water_mass, node_mixture_mass);
      gimp_nodes_(i)->assign_multimaterial_nodal_intrmd_masses(mat_id_, node_intrmd_mass);
    }
  }
}

void mpm::Particle::check_preconsolidation_pressure() {
  double M = 1;
  double p = -(solid_stress_(0)+solid_stress_(1)+solid_stress_(2))/3;
  double q = sqrt(3*(pow((solid_stress_(0)-solid_stress_(1)),2)
              + pow((solid_stress_(1)-solid_stress_(2)),2)
              + pow((solid_stress_(2)-solid_stress_(0)),2)));
  double pc = p + pow((q/M),2)/p;
  //if (pc_ < pc)
  pc_ = pc;

  //void_ratio_ = 3.3 - 0.26*log(pc_/1000);
  //porosity_ = void_ratio_ / (1 + void_ratio_);
}

void mpm::Particle::map_multimaterial_domain_gradients(){
  for (unsigned i=0; i < gimpNodes; i++){
    Eigen::Matrix<double, 1, dim> gradient;
    gradient = volume_ * grad_shape_fun_.col(i);
    gimp_nodes_(i)->assign_multimaterial_domain_gradients(mat_id_, gradient);
  }
}

void mpm::Particle::append_normal_distance_to_nodes() {
  for (unsigned i=0; i<gimpNodes; i++) {
    if(shape_fun_(i) > 0) {
      std::set<unsigned> mat_set = gimp_nodes_(i) -> give_node_material_ids();
      if (mat_set.size() == 2) {
        //n: unit normal vector, t: unit tangential vector.
        Eigen::Matrix<double,1,dim> n = gimp_nodes_(i)-> give_node_multimaterial_normal_vectors(1);
        Eigen::Matrix<double,1,dim> t;
        t(0) = n(1);
        t(1) = -n(0);
        Eigen::Matrix<double,1,dim> ncoord = gimp_nodes_(i) -> give_node_coordinates();
        double coord_normal = (coord_).dot(n);
        double tangential_distance = (coord_-ncoord).dot(t);
        // Temporary
        Eigen::Matrix<double,1,dim> Fp;
        Fp(0) = 1/(1 + solid_strain_(0));
        Fp(1) = 1/(1 + solid_strain_(1));
        Eigen::Matrix<double,1,dim> n0;
        n0(0) = Fp(0)*n(0);
        n0(1) = Fp(1)*n(1);
        double Rp = particle_radius_/sqrt(pow(n0(0),2)+pow(n0(1),2));
        double normal_distance;
        if (mat_id_ == 0)
          normal_distance = coord_normal - Rp;
        else if (mat_id_ == 1)
          normal_distance = coord_normal + Rp;
        //if (id_ == 102084 ||id_ == 102085 ||id_ == 102086 || id_ == 102408 ||id_ == 102409 ||id_ == 102410 ||id_ == 102411) {
        //  std::cout << id_ << "\n";
        //  std::cout << "n: " << n << "\n";
        //  std::cout << "coord_normal: " << coord_normal << "\n";
        //  std::cout << "R: " << Rp << "\n";
        //  std::cout << "normal_distance: " << normal_distance << "\n";
        //}
        gimp_nodes_(i)->append_node_normal_distance(mat_id_,normal_distance,tangential_distance);
      }
    }
  }
}

void mpm::Particle::compute_penalty_factor() {
  //Compute radius of particle
  double R = sqrt(volume_/3.14159);
  Eigen::Matrix<double,1,dim> elem_length = element_->give_element_length();
  double d = elem_length(0)-2*R;
  // particle spacing
  double ps = d;
  std::set<unsigned> mat_set = element_ -> give_element_material_ids();
  if (mat_set.size() == 2){
    std::vector<Eigen::Matrix<double,1,dim>> pcoord_set = element_ -> give_element_particle_coordinates();
    unsigned num_particles = pcoord_set.size();
    for (auto i=0; i < num_particles; i++){
      Eigen::Matrix<double,1,dim> pcoord = pcoord_set[i];
      double x1 = pcoord(0) - coord_(0);
      double x2 = pcoord(1) - coord_(1);
      double x = sqrt(pow(x1,2) + pow(x2,2))-2*R;
      if (x < 0)
        x = 0;
      if (x < ps)
        ps = x;
    }
  }
  for (unsigned i=0; i< numNodes; i++) {
    double temp_factor;
    double factor = nodes_(i)->give_node_penalty_factor();
    std::set<unsigned> material_ids = nodes_(i)->material_ids_;
    if (material_ids.size() == 2) {
      Eigen::Matrix<double,1,dim> ncoord = nodes_(i)->give_node_coordinates();
      double ns1 = ncoord(0) - coord_(0);
      double ns2 = ncoord(1) - coord_(1);
      double ns = sqrt(pow(ns1,2) + pow(ns2,2))-R;
      if (ns < 0)
        ns = 0;
      double s = std::min(ns,ps);
      temp_factor = 1- pow(((s)/d),3);
      //temp_factor = 1;
      if(temp_factor < 0) {
        temp_factor = 0;
      }
      if(temp_factor > 0.999) {
        temp_factor = 1;
      }
      if(temp_factor > factor){
        nodes_(i)->assign_node_penalty_factor(temp_factor);
      }
    } 
  }
}

void mpm::Particle::map_multimaterial_masses_to_nodes() {
  double node_solid_mass = 0.;
  double node_water_mass = 0.;
  double node_mixture_mass = 0.;
  double node_intrmd_mass = 0.;
  for (unsigned i = 0; i < gimpNodes; i++) {
    node_solid_mass = solid_mass_ * shape_fun_(i);
    node_water_mass = water_mass_ * shape_fun_(i);
    node_mixture_mass = (solid_mass_ + water_mass_) * shape_fun_(i);
    node_intrmd_mass = water_mass_ * shape_fun_(i) * permeability_ / (porosity_ * 9.81);
    gimp_nodes_(i)->assign_multimaterial_nodal_masses(mat_id_, node_solid_mass, node_water_mass, node_mixture_mass);
    gimp_nodes_(i)->assign_multimaterial_nodal_intrmd_masses(mat_id_, node_intrmd_mass);
  }
}

void mpm::Particle::map_sp_mass_to_nodes(bool contact) {
  double node_solid_mass = 0.;
  double node_water_mass = 0.;
  double node_mixture_mass = 0.;
  double node_intrmd_mass = 0.;
  for (unsigned i = 0; i < gimpNodes; i++) {
    node_solid_mass = solid_mass_ * shape_fun_(i);
    node_mixture_mass = solid_mass_ * shape_fun_(i);
    gimp_nodes_(i)->assign_nodal_masses(node_solid_mass, node_water_mass, node_mixture_mass);
    gimp_nodes_(i)->assign_nodal_intrmd_mass(node_intrmd_mass);
    if (contact)
    {
      gimp_nodes_(i)->assign_multimaterial_nodal_masses(mat_id_, node_solid_mass, node_water_mass, node_mixture_mass);
      gimp_nodes_(i)->assign_multimaterial_nodal_intrmd_masses(mat_id_, node_intrmd_mass);
    }
  }
}

void mpm::Particle::map_momentum_to_nodes(bool contact) {
  Eigen::Matrix<double, 1, dim> node_solid_momentum = Eigen::Matrix<double, 1, dim>::Zero();
  for (unsigned i = 0; i < gimpNodes; i++) {
    node_solid_momentum = solid_mass_ * solid_velocity_ * shape_fun_(i);
    gimp_nodes_(i)->assign_nodal_momentum(node_solid_momentum);
    if (contact)
      gimp_nodes_(i)->assign_multimaterial_nodal_momenta(mat_id_, node_solid_momentum);
  }
}


void mpm::Particle::map_pore_pressure_to_nodes() {
  double node_pore_pressure = 0.;
  for (unsigned i = 0; i < gimpNodes; i++) {
    node_pore_pressure = water_mass_ * pore_pressure_ * shape_fun_(i);
    gimp_nodes_(i)->assign_nodal_multimaterial_pore_pressures(mat_id_,node_pore_pressure);
  }
}

void mpm::Particle::map_pore_pressure_from_nodes() {
  double temp_pore_pressure = 0.;
  for (unsigned i = 0; i < gimpNodes; i++) {
    double node_pore_pressure = gimp_nodes_(i)->give_node_pressure(mat_id_);
    temp_pore_pressure += shape_fun_(i) * node_pore_pressure;
    }
  if(std::fabs(temp_pore_pressure) < 1.0E-18)
    temp_pore_pressure = 0.;
  pore_pressure_ = temp_pore_pressure;
}

void mpm::Particle::map_pressure_increment_to_nodes() {
  double node_pressure_increment = 0.;
  for (unsigned i = 0; i < numNodes; i++) {
    node_pressure_increment = water_mass_ * pressure_increment_ * shape_fun_(i);
    nodes_(i)->assign_nodal_multimaterial_pore_pressures(mat_id_, node_pressure_increment);
  }
}

void mpm::Particle::map_pressure_increment_from_nodes() {
  double temp_pressure_increment = 0.;
  for (unsigned i = 0; i < numNodes; i++) {
    double node_pressure_increment = nodes_(i)->give_node_pressure(mat_id_);
    temp_pressure_increment += shape_fun_(i) * node_pressure_increment;
  }
  if(std::fabs(temp_pressure_increment) < 1.0E-18)
    temp_pressure_increment = 0.;
  pressure_increment_ = temp_pressure_increment;
  pore_pressure_ += pressure_increment_;
}

void mpm::Particle::assign_body_force_to_nodes(bool contact, const double &time) {
  for (unsigned i = 0; i < gimpNodes; i++) {
    Eigen::Matrix<double,1,dim> node_mixture_body_force = shape_fun_(i) * (solid_mass_ + water_mass_) * gravity_ ;
    Eigen::Matrix<double,1,dim> node_water_body_force = shape_fun_(i) * water_mass_ * gravity_ * permeability_ / (porosity_ * 9.81);
    gimp_nodes_(i)->assign_body_force(node_mixture_body_force, node_water_body_force);
    if (contact)
      gimp_nodes_(i)->assign_multimaterial_body_forces(mat_id_, node_mixture_body_force, node_water_body_force);
  }
}

void mpm::Particle::assign_sp_body_force_to_nodes(bool contact) {
  Eigen::Matrix<double,1,dim> node_water_body_force = Eigen::Matrix<double,1,dim>::Zero();
  for (unsigned i=0; i <gimpNodes; i++) {
    Eigen::Matrix<double,1,dim> node_mixture_body_force = shape_fun_(i) * solid_mass_ * gravity_ ;
    gimp_nodes_(i)->assign_body_force(node_mixture_body_force, node_water_body_force);
    if (contact)
      gimp_nodes_(i)->assign_multimaterial_body_forces(mat_id_, node_mixture_body_force, node_water_body_force);
  }
}

void mpm::Particle::assign_traction_force_to_nodes(bool contact, const double &time) {
  double pi = 22./7.;
  Eigen::Matrix<double,1,dim> ntotal_traction = Eigen::Matrix<double,1,dim>::Zero();
  Eigen::Matrix<double,1,dim> nwater_traction = Eigen::Matrix<double,1,dim>::Zero();
  for (unsigned i = 0; i < gimpNodes; i++) {
    nwater_traction = Eigen::Matrix<double,1,dim>::Zero();
    for (unsigned k = 0; k < water_traction_.size(); k++) {
      unsigned direction = std::get<0>(water_traction_.at(k));
      double traction = std::get<1>(water_traction_.at(k));
      nwater_traction(direction) -= beta_scalar_ * shape_fun_(i) * porosity_ * traction * volume_ / spacing_(direction);
    }
    nwater_traction = nwater_traction * permeability_ /(porosity_ * 9.81);
    ntotal_traction = Eigen::Matrix<double,1,dim>::Zero();
    for (unsigned j = 0; j < solid_traction_.size(); j++) {
      unsigned direction = std::get<0>(solid_traction_.at(j));
      double total_traction =std::get<1>(solid_traction_.at(j));
      if (time <= 1){
        double x = time;
        total_traction = total_traction*(pow(x,3)*(10 - 15*x +6*pow(x,2)));
      }
      ntotal_traction(direction) += shape_fun_(i) * total_traction * volume_ / spacing_(direction);
    } 
    gimp_nodes_(i)->assign_traction_force(ntotal_traction, nwater_traction);
    if (contact)
      gimp_nodes_(i)->assign_multimaterial_traction_forces(mat_id_, ntotal_traction, nwater_traction);
  }
}

void mpm::Particle::assign_sp_traction_force_to_nodes(bool contact, const double &time) {
  double pi = 22./7.;
  Eigen::Matrix<double,1,dim> ntotal_traction = Eigen::Matrix<double,1,dim>::Zero();
  Eigen::Matrix<double,1,dim> nwater_traction = Eigen::Matrix<double,1,dim>::Zero();
  for (unsigned i = 0; i < gimpNodes; i++) {  
    ntotal_traction = Eigen::Matrix<double,1,dim>::Zero();
    for (unsigned j = 0; j < solid_traction_.size(); j++) {
      unsigned direction = std::get<0>(solid_traction_.at(j));
      double total_traction =std::get<1>(solid_traction_.at(j));
    //}
    //  if (time <= 1){
    //     double x = time;
    //     total_traction = total_traction*(pow(x,3)*(10 - 15*x + 6*pow(x,2)));
    //  }
    //  else
    //    total_traction = total_traction;
     //   
      ntotal_traction(direction) += shape_fun_(i) * total_traction * volume_ / spacing_(direction);
    } 
    gimp_nodes_(i)->assign_traction_force(ntotal_traction, nwater_traction);
    if (contact)
      gimp_nodes_(i)->assign_multimaterial_traction_forces(mat_id_, ntotal_traction, nwater_traction);
  }
}

void mpm::Particle::assign_internal_force_to_nodes(bool contact) {
  Eigen::Matrix<double,1,dim> ntotal_int_force = Eigen::Matrix<double,1,dim>::Zero();
  Eigen::Matrix<double,1,dim> nwater_int_force = Eigen::Matrix<double,1,dim>::Zero();
  Eigen::Matrix<double,1,dof> total_temp_stress;
  if (dim == 2) {
    total_temp_stress(0) = solid_stress_(0) - beta_scalar_ * pore_pressure_;
    total_temp_stress(1) = solid_stress_(1) - beta_scalar_ * pore_pressure_;
    total_temp_stress(2) = solid_stress_(3);
  }
  //std::cout << volume_ << "\n";
  //if (id_ == 596) {
    //std::cout << total_temp_stress << "\n";
  //}
  Eigen::Matrix<double, dim, 1> int_force_t, int_force_w;
  for (unsigned i = 0; i < gimpNodes; i++) {
    Eigen::Matrix<double,dof,dim> Bi = B_.at(i);
    int_force_t = - volume_ * (Bi.transpose() * total_temp_stress.transpose());
    int_force_w = beta_scalar_ * volume_ * grad_shape_fun_.col(i) * porosity_ * pore_pressure_;

    ntotal_int_force = int_force_t.transpose();
    nwater_int_force = int_force_w * (permeability_ / (porosity_ * 9.81));
    gimp_nodes_(i)->assign_internal_force(ntotal_int_force, nwater_int_force);
    if (contact)
      gimp_nodes_(i) -> assign_multimaterial_internal_forces(mat_id_, ntotal_int_force, nwater_int_force);
  }
}

void mpm::Particle::assign_sp_internal_force_to_nodes(bool contact) {
  Eigen::Matrix<double,1,dim> ntotal_int_force = Eigen::Matrix<double,1,dim>::Zero();
  Eigen::Matrix<double,1,dim> nwater_int_force = Eigen::Matrix<double,1,dim>::Zero();
  Eigen::Matrix<double,1,dof> total_temp_stress;
  if (dim == 2) {
    total_temp_stress(0) = solid_stress_(0);
    total_temp_stress(1) = solid_stress_(1);
    total_temp_stress(2) = solid_stress_(3);
  }

  Eigen::Matrix<double, dim, 1> int_force_t;
  for (unsigned i = 0; i < gimpNodes; i++) {
    Eigen::Matrix<double,dof,dim> Bi = B_.at(i);
    int_force_t = - volume_ * (Bi.transpose() * total_temp_stress.transpose());
    ntotal_int_force = int_force_t.transpose();
    gimp_nodes_(i)->assign_internal_force(ntotal_int_force, nwater_int_force);
    if (contact)
      gimp_nodes_(i)->assign_multimaterial_internal_forces(mat_id_, ntotal_int_force, nwater_int_force);
  }
}

void mpm::Particle::compute_solid_strain_rate(bool contact) {
  Eigen::Matrix<double, dim, 1> node_solid_velocity;
  solid_strain_rate_ = Eigen::Matrix<double,1,dof>::Zero();
  for (unsigned i = 0; i < gimpNodes; i++) {
    Eigen::Matrix<double,1,dim> velocity = Eigen::Matrix<double,1,dim>::Zero();
    Eigen::Matrix<double,dof,dim> Bi = B_.at(i);
    //Eigen::Matrix<double,dof,dim> BBari = BBar_.at(i);
    if (contact)
      velocity = gimp_nodes_(i)->give_node_solid_multimaterial_final_velocity(mat_id_);
    else
      velocity = gimp_nodes_(i)->give_node_solid_final_velocity();
    node_solid_velocity = velocity.transpose();
    solid_strain_rate_ += (Bi * node_solid_velocity);
    //solid_strain_rate_ += (BBari * node_solid_velocity);
  }  
}


void mpm::Particle::compute_solid_strain_rate_at_elem_center(bool contact) {
  Eigen::Matrix<double, dim, 1> node_solid_velocity;
  solid_centre_strain_rate_ = Eigen::Matrix<double,1,dof>::Zero();
  for (unsigned i = 0; i < numNodes; i++) {
    Eigen::Matrix<double,dof,dim> Bi_centre = BCentre_.at(i);
    Eigen::Matrix<double,1,dim> velocity = nodes_(i)->give_node_solid_multimaterial_final_velocity(mat_id_);
    node_solid_velocity = velocity.transpose();
    solid_centre_strain_rate_ += (Bi_centre * node_solid_velocity);
  }
}


void mpm::Particle::compute_solid_strain(const double dt) {
  solid_strain_ += (solid_strain_rate_ * dt);
}


void mpm::Particle::compute_solid_stress(const double dt) {
  Eigen::Matrix<double,1,dof> solid_dstrain = (dt * solid_strain_rate_);
  material_->compute_stress(solid_dstrain, solid_stress_, plastic_strain_, pc_, void_ratio_);
}

void mpm::Particle::update_contact_velocity_and_position(const double& dt) {
  Eigen::Matrix<double,1,dim> temp_solid_final_acceleration = Eigen::Matrix<double,1,dim>::Zero();
  Eigen::Matrix<double,1,dim> temp_solid_final_velocity = Eigen::Matrix<double,1,dim>::Zero();
  Eigen::Matrix<double,1,dim> temp_solid_inc_velocity = Eigen::Matrix<double,1,dim>::Zero();
  for (unsigned i = 0; i < gimpNodes; i++) {
    Eigen::Matrix<double,1,dim> nsolid_acceleration = gimp_nodes_(i)->give_node_solid_multimaterial_final_acceleration(mat_id_);
    Eigen::Matrix<double,1,dim> nsolid_final_velocity = gimp_nodes_(i)->give_node_solid_multimaterial_final_velocity(mat_id_);
    Eigen::Matrix<double,1,dim> nsolid_inc_velocity = gimp_nodes_(i)->give_node_solid_multimaterial_increment_velocity(mat_id_);
    temp_solid_final_acceleration += shape_fun_(i) * nsolid_acceleration;
    temp_solid_final_velocity += shape_fun_(i) * nsolid_final_velocity;
    temp_solid_inc_velocity += shape_fun_(i) * nsolid_inc_velocity;
  }
  for (unsigned j = 0; j < dim; j++) {
    if (std::fabs(temp_solid_final_acceleration(j)) < 1.0E-18)
      temp_solid_final_acceleration(j) = 0.;
    if (std::fabs(temp_solid_final_velocity(j)) < 1.0E-18)
      temp_solid_final_velocity(j) = 0.;
  }

  coord_ += (dt * temp_solid_final_velocity);
  displacement_ += (dt * temp_solid_final_velocity);
  //solid_velocity_ += (dt * temp_solid_final_acceleration);
  solid_velocity_ = temp_solid_final_velocity;
  du_ = dt*temp_solid_final_velocity;
}

void mpm::Particle::update_velocity_and_position(const double& dt,
						 const double& time) {
  Eigen::Matrix<double,1,dim> temp_solid_final_acceleration = Eigen::Matrix<double,1,dim>::Zero();
  Eigen::Matrix<double,1,dim> temp_solid_final_velocity = Eigen::Matrix<double,1,dim>::Zero();
  Eigen::Matrix<double,1,dim> temp_solid_inc_velocity = Eigen::Matrix<double,1,dim>::Zero();
  for (unsigned i = 0; i < gimpNodes; i++) {
    Eigen::Matrix<double,1,dim> nsolid_acceleration = gimp_nodes_(i)->give_node_solid_final_acceleration();
    Eigen::Matrix<double,1,dim> nsolid_final_velocity = gimp_nodes_(i)->give_node_solid_final_velocity();
    Eigen::Matrix<double,1,dim> nsolid_inc_velocity = gimp_nodes_(i)->give_node_solid_increment_velocity();
    temp_solid_final_acceleration += shape_fun_(i) * nsolid_acceleration;
    temp_solid_final_velocity += shape_fun_(i) * nsolid_final_velocity;
    temp_solid_inc_velocity += shape_fun_(i) * nsolid_inc_velocity;
  }

  for (unsigned j = 0; j < dim; j++) {
    if (std::fabs(temp_solid_final_acceleration(j)) < 1.0E-18)
      temp_solid_final_acceleration(j) = 0.;
    if (std::fabs(temp_solid_final_velocity(j)) < 1.0E-18)
      temp_solid_final_velocity(j) = 0.;
  }

  coord_ += (dt * temp_solid_final_velocity);
  displacement_ += (dt * temp_solid_final_velocity);
  //solid_velocity_ += (dt * temp_solid_final_acceleration);
  solid_velocity_ = temp_solid_final_velocity;
  //solid_velocity_ += temp_solid_inc_velocity;
}

void mpm::Particle::update_pressure() {
  double temp_pressure = 0.;
  if (nphase_ == 2) {
    for (unsigned i = 0; i < gimpNodes; i++) {
      double npressure = gimp_nodes_(i)->give_node_multimaterial_pressure_increments(mat_id_);
      temp_pressure += shape_fun_(i) * npressure;
    }
  }
  if (std::fabs(temp_pressure) < 1.0E-18)
    temp_pressure = 0.;
  pressure_increment_ = temp_pressure;
  pore_pressure_ = pore_pressure_ + pressure_increment_;
}

void mpm::Particle::update_density_and_mass(const double &compressibility) {

}

void mpm::Particle::update_porosity(const double &dt) {
  double dvol_strain = dt * (solid_strain_rate_(0) + solid_strain_rate_(1));
  porosity_ = 1 - ((1 - porosity_) / (1 + dvol_strain));
  void_ratio_ = porosity_/(1-porosity_);
  volume_ = volume_ * (1 + dvol_strain);
  for (auto i=0; i<dim; i++) {
    particle_size_(i) = particle_size_(i) * (1 + solid_strain_rate_(i)*dt);
  }
  water_mass_ = porosity_ * water_grain_density_ * volume_;
}

void mpm::Particle::compute_total_solid_stress() {
  total_solid_stress_(0) = solid_stress_(0)-pore_pressure_;
  total_solid_stress_(1) = solid_stress_(1)-pore_pressure_;
  total_solid_stress_(3) = solid_stress_(3);
}

void mpm::Particle::check_failure() {
  failure_ = 0;
  t_ = (solid_stress_(0) + solid_stress_(1)) * 0.5;
  s_ = std::sqrt(solid_stress_(3)*solid_stress_(3) + std::pow((solid_stress_(1)-solid_stress_(2))/2,2));
  if (s_ > 10000)
    failure_ = 1;
}

void mpm::Particle::write_velocity(std::ostream& oFile) {
    if (dim == 2)
        oFile << solid_velocity_(0) << " " << solid_velocity_(1) << " " << "0" << "\n";
    if (dim == 3)
        oFile << solid_velocity_(0) << " " << solid_velocity_(1) << " " << solid_velocity_(2) << "\n";
}

void mpm::Particle::write_t_s(std::ostream& oFile){
  oFile << shape_fun_(2) << " " << shape_fun_(11) << "\n" << B_.at(2) << "\n" << B_.at(11) << "\n" << particle_size_(1) << "\n" << 
  volume_ << "\n";
}

void mpm::Particle::write_pressure(std::ostream& oFile) {
    oFile << pore_pressure_ << "\n";
}


void mpm::Particle::write_stress(std::ostream& oFile) {
    if (dim == 2) {
      oFile << solid_stress_(0) << " " << solid_stress_(1) << " " << solid_stress_(3) << "\n";
    }
    if (dim == 3)
        oFile << solid_stress_(0) << " " << solid_stress_(1) << " " << solid_stress_(2) << " " << solid_stress_(3) << " " << solid_stress_(4) << " " << solid_stress_(5) << "\n";
}

void mpm::Particle::write_total_stress(std::ostream& oFile) {
    if (dim == 2) {
      oFile << total_solid_stress_(0) << " " << total_solid_stress_(1) << " " << total_solid_stress_(3) << "\n";
    }
    if (dim == 3)
        oFile << solid_stress_(0) << " " << solid_stress_(1) << " " << solid_stress_(2) << " " << solid_stress_(3) << " " << solid_stress_(4) << " " << solid_stress_(5) << "\n";
}


void mpm::Particle::write_strain(std::ostream& oFile) {
    if (dim == 2)
        oFile << solid_strain_(0) << " " << solid_strain_(1) << " " << solid_strain_(2) << "\n";
    if (dim == 3)
        oFile << solid_strain_(0) << " " << solid_strain_(1) << " " << solid_strain_(2) << " " << solid_strain_(3) << " " << solid_strain_(4) << " " << solid_strain_(5) << "\n";
}

void mpm::Particle::write_displacement(std::ostream& oFile) {
    if (dim == 2)
        oFile << displacement_(0) << " " << displacement_(1) << " " << 0. << "\n";
}

void mpm::Particle::write_deviatoric_shear_strain(std::ofstream& oFile) {
  oFile << dev_shear_strain_  << "\n";
}

void mpm::Particle::write_equivalent_plastic_strain(std::ofstream& oFile) {
    //oFile << plastic_strain_  << "\n";
    oFile << pc_ << " " << "\n";
}


void mpm::Particle::compute_local_coordinates() {
    Eigen::Matrix<double,1,dim> elem_centre_coord = element_->give_element_centre_coord();
    Eigen::Matrix<double,1,dim> elem_length        = element_ -> give_element_length();
    for (unsigned i = 0; i < dim; i++) {
        xi_(i) = 2. * (coord_(i) - elem_centre_coord(i)) / elem_length(i);
        if (((std::fabs(xi_(i)) > 0.999999) && (std::fabs(xi_(i)) < 1.)) || (std::fabs(xi_(i)) > 1.))
            sign(xi_(i), 1.);
        else if ((std::fabs(xi_(i)) > 0.) && (std::fabs(xi_(i)) < 0.000001))
            xi_(i) = 0.;
    }
}


void mpm::Particle::compute_shape_functions() {
    Eigen::Matrix<double,16,dim> local_nodes = 
    (Eigen::Matrix<double, 16, dim>() << -1., -1.,
                                      1., -1.,
                                      1.,  1.,
                                     -1.,  1.,
                                     -3., -3.,
                                     -1., -3.,
                                      1., -3.,
                                      3., -3.,
                                      3., -1.,
                                      3.,  1.,
                                      3.,  3.,
                                      1.,  3.,
                                     -1.,  3.,
                                     -3.,  3.,
                                     -3.,  1.,
                                     -3., -1.).finished();
    
    //! length of element in local coordinate
    double element_length = 2;
    Eigen::Matrix<double,1,dim> L = element_-> give_element_length();
    //! loop to iterate over nodes
    for (unsigned n = 0; n < gimpNodes; ++n) {
      //! local shape function
      Eigen::Matrix<double,1,dim> sni;
      //! loop to iterate over dimensions
      for (unsigned i = 0; i < dim; ++i) {
        double lp = particle_size_(i)*0.5*element_length/L(i);
        double ni = local_nodes(n,i);
        double npni = xi_(i) - ni; // local particle - local node
        //! Conditional shape function statement see: Bardenhagen 2004
        if (npni <= (-element_length - lp)) {
          sni(i) = 0.;
        } else if ((-element_length - lp) < npni &&
                   npni <= (-element_length + lp)) {
          sni(i) = std::pow(element_length + lp + npni, 2.) /
                   (4. * (element_length * lp));
        } else if ((-element_length + lp) < npni && npni <= -lp) {
          sni(i) = 1. + (npni / element_length);
        } else if (-lp < npni && npni <= lp) {
          sni(i) =
              1. - (((npni * npni) + (lp * lp)) / (2. * element_length * lp));
        } else if (lp < npni && npni <= (element_length - lp)) {
          sni(i) = 1. - (npni / element_length);
        } else if ((element_length - lp) < npni &&
                   npni <= element_length + lp) {
          sni(i) = std::pow(element_length + lp - npni, 2.) /
                   (4. * element_length * lp);
        } else if ((element_length + lp) < npni) {
          sni(i) = 0.;
        } else {
          throw std::runtime_error(
              "GIMP shapefn: Point location outside area of influence");
        }
      }
      shape_fun_(n) = sni(0) * sni(1);  // See: Pruijn, N.S., 2016. Eq(4.30)
      if (shape_fun_(n) < 1E-18)
        shape_fun_(n) = 0;
    }
    //if (id_ == 602) {
    //  std::cout << "Shape functions of particle: " << id_ << ":";
    //  std::cout << "\t" << shape_fun_ << "\n";
      //std::cout << "Particle size: " << particle_size_ << "\n";
    //}
}


void mpm::Particle::compute_global_derivatives_shape_functions() {
    Eigen::Matrix<double,16,dim> local_nodes = 
    (Eigen::Matrix<double, 16, dim>() << -1., -1.,
                                      1., -1.,
                                      1.,  1.,
                                     -1.,  1.,
                                     -3., -3.,
                                     -1., -3.,
                                      1., -3.,
                                      3., -3.,
                                      3., -1.,
                                      3.,  1.,
                                      3.,  3.,
                                      1.,  3.,
                                     -1.,  3.,
                                     -3.,  3.,
                                     -3.,  1.,
                                     -3., -1.).finished();
    
    //! length of element in local coordinate
    double element_length = 2;
    Eigen::Matrix<double,1,dim> L = element_->give_element_length();
    //! loop to iterate over nodes
    for (unsigned n = 0; n < gimpNodes; ++n) {
      //! local shape function
      Eigen::Matrix<double,1,dim> sni;
      //! local grad shape function
      Eigen::Matrix<double,1,dim> dni;
      //! loop to iterate over dimensions
      for (unsigned i = 0; i < dim; ++i) {
        double lp = particle_size_(i)*0.5*element_length/L(i);
        double ni = local_nodes(n,i);
        double npni = xi_(i) - ni; // local particle - local node
        //! Conditional shape function statement see: Bardenhagen 2004
        if (npni <= (-element_length - lp)) {
          sni(i) = 0.;
          dni(i) = 0.;
        } else if ((-element_length - lp) < npni &&
                   npni <= (-element_length + lp)) {

          sni(i) = std::pow(element_length + lp + npni, 2.) /
                   (4. * (element_length * lp));
          dni(i) = (element_length + lp + npni) / (2. * element_length * lp);
        } else if ((-element_length + lp) < npni && npni <= -lp) {
          sni(i) = 1. + (npni / element_length);
          dni(i) = 1. / element_length;
        } else if (-lp < npni && npni <= lp) {
          sni(i) =
              1. - (((npni * npni) + (lp * lp)) / (2. * element_length * lp));
          dni(i) = -(npni / (element_length * lp));
        } else if (lp < npni && npni <= (element_length - lp)) {
          sni(i) = 1. - (npni / element_length);
          dni(i) = -(1. / element_length);
        } else if ((element_length - lp) < npni &&
                   npni <= (element_length + lp)) {
          sni(i) = std::pow(element_length + lp - npni, 2.) /
                   (4. * element_length * lp);
          dni(i) = -((element_length + lp - npni) / (2. * element_length * lp));
        } else if ((element_length + lp) < npni) {
          sni(i) = 0.;
          dni(i) = 0.;
        } else {
          throw std::runtime_error(
              "GIMP grad shapefn: Point location outside area of influence");
        }
      }
      grad_shape_fun_(0,n) = dni(0) * sni(1) * 2 / L(0);
      grad_shape_fun_(1,n) = dni(1) * sni(0) * 2 / L(1);
      // See: Pruijn, N.S., 2016. Eq(4.32)
    }
}


void mpm::Particle::compute_global_derivatives_shape_functions_at_centre() {
    Eigen::Matrix<double,1,dim> L = element_ -> give_element_length();
    if (dim == 2) {
        grad_shape_fun_centre_(0, 0) = -0.5 / L(0);
        grad_shape_fun_centre_(0, 1) =  0.5 / L(0);
        grad_shape_fun_centre_(0, 2) =  0.5 / L(0);
        grad_shape_fun_centre_(0, 3) = -0.5 / L(0);

        grad_shape_fun_centre_(1, 0) = -0.5 / L(1);
        grad_shape_fun_centre_(1, 1) = -0.5 / L(1);
        grad_shape_fun_centre_(1, 2) =  0.5 / L(1);
        grad_shape_fun_centre_(1, 3) =  0.5 / L(1);
    }
    else if (dim == 3) { }
}


void mpm::Particle::compute_B_matrix() {
    Eigen::Matrix<double,dof,dim> Bi_;
    for (unsigned i = 0; i < gimpNodes; i++) {
        if (dim == 2) {
            Bi_(0,0) = grad_shape_fun_(0,i); 
            Bi_(0,1) = 0.;
            Bi_(1,0) = 0.;              
            Bi_(1,1) = grad_shape_fun_(1,i);
            Bi_(2,0) = grad_shape_fun_(1,i); 
            Bi_(2,1) = grad_shape_fun_(0,i);
        }
        else if (dim == 3) {
            Bi_(0,0) = grad_shape_fun_(0,i); 
            Bi_(0,1) = 0.; 
            Bi_(0,2) = 0.;
            Bi_(1,0) = 0.; 
            Bi_(1,1) = grad_shape_fun_(1,i);
            Bi_(1,2) = 0.;
            Bi_(2,0) = 0.; 
            Bi_(2,1) = 0.;
            Bi_(2,2) = grad_shape_fun_(2,i);
            Bi_(3,0) = grad_shape_fun_(1,i); 
            Bi_(3,1) = grad_shape_fun_(0,i); 
            Bi_(3,2) = 0; 
            Bi_(4,0) = 0; 
            Bi_(4,1) = grad_shape_fun_(2,i); 
            Bi_(4,2) = grad_shape_fun_(1,i);
            Bi_(5,0) = grad_shape_fun_(2,i);
            Bi_(5,1) = 0.;
            Bi_(5,2) = grad_shape_fun_(0,i);
        }
        B_.at(i) = Bi_;
      //if (mat_id_ == 0)
      //  std::cout << "B_(" << i << "): " << B_.at(i) << "\n";
    }
}


void mpm::Particle::compute_B_matrix_at_centre() {
    Eigen::Matrix<double,dof,dim> BiCentre_;
    for (unsigned i = 0; i < numNodes; i++) {
        if (dim == 2) {
            BiCentre_(0,0) = grad_shape_fun_centre_(0,i); 
            BiCentre_(0,1) = 0.;
            BiCentre_(1,0) = 0.;              
            BiCentre_(1,1) = grad_shape_fun_centre_(1,i);
            BiCentre_(2,0) = grad_shape_fun_centre_(1,i); 
            BiCentre_(2,1) = grad_shape_fun_centre_(0,i);
        }
        BCentre_.at(i) = BiCentre_;
    }
}

void mpm::Particle::compute_BBar_matrix() {

    Eigen::Matrix<double,dof,dim> BBari_;
    double BBar00, BBar01, BBar10, BBar11;
    BBar00 = BBar01 = BBar10 = BBar11 = 0.;
    const double OneThird = 1.0/3.0;
    if (dim == 2) {
        for (unsigned i = 0; i < numNodes; i ++) {
            BBar00 = OneThird * (2*grad_shape_fun_(0,i)+grad_shape_fun_centre_(0,i));
            BBar01 = OneThird * (grad_shape_fun_centre_(1,i) - grad_shape_fun_(1,i));
            BBar10 = OneThird * (grad_shape_fun_centre_(0,i) - grad_shape_fun_(0,i));
            BBar11 = OneThird * (2*grad_shape_fun_(1,i)+grad_shape_fun_centre_(1,i));

            BBari_(0,0)= BBar00;
            BBari_(0,1)= BBar01;
            BBari_(1,0)= BBar10;
            BBari_(1,1)= BBar11;


            BBari_(2,0)= grad_shape_fun_(1,i);
            BBari_(2,1)= grad_shape_fun_(0,i);
            BBar_.at(i) = BBari_;
        }
    }
}

void mpm::Particle::compute_element_matrix_mp(const double& dt) {
  double L_entry;
  Eigen::Matrix<double, 1, dim> K_L2_entry, K_L3_entry, K_S2_entry;
  double mixture_density = (porosity_ * water_grain_density_) + ((1 - porosity_)* solid_grain_density_);
  double water_density = porosity_ * water_grain_density_;
  double laplace_constant = ((permeability_ / (water_density * 9.81))*((water_density / mixture_density)-porosity_)) - (dt/mixture_density);
  //if(id_ > 49999) {
  
  if (dim == 2) {
    for (unsigned i = 0; i < gimpNodes; i++) {
      for (unsigned j =0; j < gimpNodes; j++) {
        L_entry =  laplace_constant * volume_ * ((grad_shape_fun_(0,i) * grad_shape_fun_(0,j)) + (grad_shape_fun_(1,i) * grad_shape_fun_(1,j)));
	
	      K_L2_entry(0) = volume_ * shape_fun_(i) * grad_shape_fun_(0,j);
        K_L2_entry(1) = volume_ * shape_fun_(i) * grad_shape_fun_(1,j);
	
        K_L3_entry(0) = porosity_ * volume_ * grad_shape_fun_(0,i) * shape_fun_(j); 
        K_L3_entry(1) = porosity_ * volume_ * grad_shape_fun_(1,i) * shape_fun_(j);

        K_S2_entry(0) = volume_ * shape_fun_(i) * grad_shape_fun_(0,j);
        K_S2_entry(1) = volume_ * shape_fun_(i) * grad_shape_fun_(1,j);

        element_->build_Laplace_matrix_mp(i,j,L_entry);  
        element_->build_K_L2_matrix_mp(i,j,K_L2_entry);
        element_->build_K_L3_matrix_mp(i,j,K_L3_entry);
        element_->build_K_S2_matrix_mp(i,j,K_S2_entry);
      }
    }
  }
}


void mpm::Particle::sign(double& variable, double value) {
    if (variable < 0)
        variable = -value;
    if (variable > 0)
        variable = value;
}
