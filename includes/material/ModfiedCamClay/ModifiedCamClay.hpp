/***************************************************************************
                          Material Point Method
                           Shyamini Kularathna
                         Unversity Of Cambridge
File: ModifiedCamClay.hpp
****************************************************************************/

#ifndef MPM_MATERIAL_MODIFIEDCAMCLAY_H
#define MPM_MATERIAL_MODIFIEDCAMCLAY_H

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// c++ header files
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>

#include <boost/math/constants/constants.hpp>

// mpm header files
#include "MaterialBase.hpp"
#include "PropertyParse.hpp"
#include "Constants.hpp" 
#include "Material_utility.hpp"

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

namespace mpm {
    namespace material {
        class ModifiedCamClay;
    }
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class mpm::material::ModifiedCamClay : public mpm::material::MaterialBase {
  const double pi = boost::math::constants::pi<double>();
  const double onethirdpi = 1.047197551;
public:
  static const unsigned dim = mpm::constants::DIM;
  static const unsigned numNodes = mpm::constants::NUMNODES;
  static const unsigned dof = 3 * (dim - 1);

  typedef Eigen::Matrix<double, dof,1> VectorDDOF;
  typedef Eigen::Matrix<double, 6,1>   VectorD6x1;
  typedef Eigen::Matrix<double, 1,6>   VectorD1x6;
  typedef Eigen::Matrix<double, 6, 6> Matrix6x6;

  //! Failure state
  enum FailureState { Elastic = 0, Yield = 1 };

public:
  // CONSTRUCTOR
  ModifiedCamClay();

  //! CREATE THE MATERIAL
  static MaterialBase* create() {
    return new ModifiedCamClay();
  }

  //! GIVE DENSITY
  double give_density() {
    return density_;
  }

    // give porosity
  //! param[out] porosity_
  double give_porosity() {
    return porosity_;
  }

  // give permeability
  //! param[out] permeability_
  double give_permeability() {
    return permeability_;
  }

  // give void ratio
  //! param[out] void_ratio_
  double give_void_ratio() {
    return void_ratio_;
  }

  // give void ratio
  //! param[out] pc_
  double give_pc() {
    return pc_;
  }

  //! compute elastic stiffness matrix
  void compute_elastic_stiffness_matrix();
  void compute_elastic_tensor();
  void change_friction_coefficient() {}

  //! COMPUTE STRESS
  void compute_stress(const Eigen::Matrix<double,1,dof>& dstrain,
		      Eigen::Matrix<double,1,6>& stress, double& epds, double& pc, double& void_ratio);

  //! Compute yield function and yield state
  //! \param[in] 
  //! \retval yield_type Yield type (elastic, shear or tensile)
  mpm::material::ModifiedCamClay::FailureState compute_yield_state();

private:  
  void compute_stress_invariants(const VectorD6x1 stress);
  void compute_df_dmul(double* df_dmul);
  void compute_dg_dpc(const double pc_n, const double p_trial, double* g_function, double* dg_dpc);
  void compute_lambda_trial();


protected:
  Matrix6x6 de_;
  double density_;
  double porosity_;
  double permeability_;

  // plastic deviatori strain vector
  double epds_;
  Eigen::Matrix<double,6,1> PDS_;

  // parameters for yield function
  double p_;
  double q_;
  double theta_;
  double pc_;
  double void_ratio_;
  double delta_phi_;
  double m_theta_;
  double f_function_;
  double subloading_r_;

  // Incremental plastic strain
  // Incremental plastic volumetic strain
  double dpvstrain_;
  // Incremental plastic deviatoric strain
  double dpdstrain_;
  // Plastic volumetic strain
  double pvstrain_;
  // Plastic deviatoric strain
  double pdstrain_;


  double youngs_modulus_{std::numeric_limits<double>::max()};
  //! Bulk modulus
  double bulk_modulus_{std::numeric_limits<double>::max()};
  //! Shear modulus
  double shear_modulus_{std::numeric_limits<double>::max()};
  //! Poisson ratio
  double poisson_ratio_{std::numeric_limits<double>::max()};
  double e_ref_{std::numeric_limits<double>::max()};
  double p_ref_{std::numeric_limits<double>::max()};
  double ocr_{std::numeric_limits<double>::max()};
  double pc0_{std::numeric_limits<double>::max()};
  double m_{std::numeric_limits<double>::max()};
  double lambda_{std::numeric_limits<double>::max()};
  double kappa_{std::numeric_limits<double>::max()};
  double e0_;
};

#include "ModifiedCamClay.ipp"

#endif
