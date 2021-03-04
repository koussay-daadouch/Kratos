//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Bodhinanda Chandra, Ilaria Iaconeta
//
//


#if !defined(KRATOS_PARTICLE_MECHANICS_APPLICATION_VARIABLES_H_INCLUDED )
#define  KRATOS_PARTICLE_MECHANICS_APPLICATION_VARIABLES_H_INCLUDED

// Project includes
#include "includes/define.h"
#include "includes/constitutive_law.h"
#include "includes/variables.h"
#include "includes/deprecated_variables.h"
#include "includes/mat_variables.h"


namespace Kratos
{

    typedef array_1d<double,3> Vector3;
    typedef array_1d<double,6> Vector6;

    // Variables definition

    /* MATERIAL POINT ELEMENTS VARIABLES */
    // Indexing
    KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, int, MP_MATERIAL_ID )
    KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, int, PARTICLES_PER_ELEMENT )
    KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, int, MP_SUB_POINTS)

    // Physical
    KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, double, MP_MASS )
    KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, double, MP_DENSITY )
    KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, double, MP_VOLUME )
    KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, double, MP_RADIUS)
    KRATOS_DEFINE_APPLICATION_VARIABLE(PARTICLE_MECHANICS_APPLICATION, bool, IS_COMPRESSIBLE)
    KRATOS_DEFINE_APPLICATION_VARIABLE(PARTICLE_MECHANICS_APPLICATION, double, MP_TEMPERATURE)

    // Energy
    KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, double, MP_POTENTIAL_ENERGY )
    KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, double, MP_KINETIC_ENERGY )
    KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, double, MP_STRAIN_ENERGY )
    KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, double, MP_TOTAL_ENERGY )
    KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, double, MP_TEMPERATURE)
    KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, double, MP_DAMAGE)

    // Pressure
    KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, double, MP_PRESSURE )
    KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, double, PRESSURE_REACTION )
    KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, double, NODAL_MPRESSURE )
    KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, double, AUX_MASS)

    // Position and kinematics
    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(PARTICLE_MECHANICS_APPLICATION, MP_COORD )
    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(PARTICLE_MECHANICS_APPLICATION, MP_DISPLACEMENT )
    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(PARTICLE_MECHANICS_APPLICATION, MP_VELOCITY )
    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(PARTICLE_MECHANICS_APPLICATION, MP_ACCELERATION )
    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(PARTICLE_MECHANICS_APPLICATION, MP_VOLUME_ACCELERATION )
    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(PARTICLE_MECHANICS_APPLICATION, AUX_MOMENTA)
    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(PARTICLE_MECHANICS_APPLICATION, AUX_INERTIA)
    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(PARTICLE_MECHANICS_APPLICATION, AUX_RESIDUAL)

    // Stress Measures
    KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, Vector, MP_CAUCHY_STRESS_VECTOR )
    KRATOS_DEFINE_APPLICATION_VARIABLE(PARTICLE_MECHANICS_APPLICATION, double, MP_EQUIVALENT_STRESS)
    KRATOS_DEFINE_APPLICATION_VARIABLE(PARTICLE_MECHANICS_APPLICATION, double, MP_HARDENING_RATIO)

    // Strain Measures
    KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, Vector, MP_ALMANSI_STRAIN_VECTOR )
    KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, double, MP_DELTA_PLASTIC_STRAIN )
    KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, double, MP_DELTA_PLASTIC_VOLUMETRIC_STRAIN )
    KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, double, MP_DELTA_PLASTIC_DEVIATORIC_STRAIN )
    KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, double, MP_EQUIVALENT_PLASTIC_STRAIN )
    KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, double, MP_EQUIVALENT_PLASTIC_STRAIN_RATE)
    KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, double, MP_ACCUMULATED_PLASTIC_VOLUMETRIC_STRAIN )
    KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, double, MP_ACCUMULATED_PLASTIC_DEVIATORIC_STRAIN )

    // Constitutive law
    KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, ConstitutiveLaw::Pointer, CONSTITUTIVE_LAW_POINTER )
    // CL: Solid
    KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, double, RAYLEIGH_ALPHA )
    KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, double, RAYLEIGH_BETA )
    // CL: Mohr Coulomb
    KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, double, COHESION )
    KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, double, INTERNAL_DILATANCY_ANGLE )
    // CL: Mohr Coulomb Strain Softening
    KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, double, INTERNAL_FRICTION_ANGLE_RESIDUAL )
    KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, double, COHESION_RESIDUAL )
    KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, double, INTERNAL_DILATANCY_ANGLE_RESIDUAL )
    KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, double, SHAPE_FUNCTION_BETA )
    // CL: Johnson Cook
    KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, double, REFERENCE_STRAIN_RATE)
    KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, double, TAYLOR_QUINNEY_COEFFICIENT)
    KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, double, FAILURE_PLASTIC_STRAIN)
    // CL: d+/d- Damage
    KRATOS_DEFINE_APPLICATION_VARIABLE(PARTICLE_MECHANICS_APPLICATION, double, DAMAGE_TENSION)
    KRATOS_DEFINE_APPLICATION_VARIABLE(PARTICLE_MECHANICS_APPLICATION, double, UNIAXIAL_STRESS_TENSION)
    KRATOS_DEFINE_APPLICATION_VARIABLE(PARTICLE_MECHANICS_APPLICATION, double, UNIAXIAL_STRAIN_TENSION)
    KRATOS_DEFINE_APPLICATION_VARIABLE(PARTICLE_MECHANICS_APPLICATION, double, UNIAXIAL_DAMAGED_STRESS_TENSION)
    KRATOS_DEFINE_APPLICATION_VARIABLE(PARTICLE_MECHANICS_APPLICATION, double, THRESHOLD_TENSION)
    KRATOS_DEFINE_APPLICATION_VARIABLE(PARTICLE_MECHANICS_APPLICATION, double, DAMAGE_COMPRESSION)
    KRATOS_DEFINE_APPLICATION_VARIABLE(PARTICLE_MECHANICS_APPLICATION, double, UNIAXIAL_STRESS_COMPRESSION)
    KRATOS_DEFINE_APPLICATION_VARIABLE(PARTICLE_MECHANICS_APPLICATION, double, UNIAXIAL_STRAIN_COMPRESSION)
    KRATOS_DEFINE_APPLICATION_VARIABLE(PARTICLE_MECHANICS_APPLICATION, double, UNIAXIAL_DAMAGED_STRESS_COMPRESSION)
    KRATOS_DEFINE_APPLICATION_VARIABLE(PARTICLE_MECHANICS_APPLICATION, double, THRESHOLD_COMPRESSION)
    KRATOS_DEFINE_APPLICATION_VARIABLE(PARTICLE_MECHANICS_APPLICATION, double, YIELD_STRESS_TENSION)
    KRATOS_DEFINE_APPLICATION_VARIABLE(PARTICLE_MECHANICS_APPLICATION, double, FRACTURE_ENERGY_TENSION)
    KRATOS_DEFINE_APPLICATION_VARIABLE(PARTICLE_MECHANICS_APPLICATION, double, DAMAGE_ONSET_STRESS_COMPRESSION)
    KRATOS_DEFINE_APPLICATION_VARIABLE(PARTICLE_MECHANICS_APPLICATION, double, YIELD_STRESS_COMPRESSION)
    KRATOS_DEFINE_APPLICATION_VARIABLE(PARTICLE_MECHANICS_APPLICATION, double, RESIDUAL_STRESS_COMPRESSION)
    KRATOS_DEFINE_APPLICATION_VARIABLE(PARTICLE_MECHANICS_APPLICATION, double, YIELD_STRAIN_COMPRESSION)
    KRATOS_DEFINE_APPLICATION_VARIABLE(PARTICLE_MECHANICS_APPLICATION, double, FRACTURE_ENERGY_COMPRESSION)
    KRATOS_DEFINE_APPLICATION_VARIABLE(PARTICLE_MECHANICS_APPLICATION, double, BIAXIAL_COMPRESSION_MULTIPLIER)
    KRATOS_DEFINE_APPLICATION_VARIABLE(PARTICLE_MECHANICS_APPLICATION, double, INTEGRATION_IMPLEX)
    KRATOS_DEFINE_APPLICATION_VARIABLE(PARTICLE_MECHANICS_APPLICATION, double, BEZIER_CONTROLLER_S1)
    KRATOS_DEFINE_APPLICATION_VARIABLE(PARTICLE_MECHANICS_APPLICATION, double, BEZIER_CONTROLLER_EP1)
    KRATOS_DEFINE_APPLICATION_VARIABLE(PARTICLE_MECHANICS_APPLICATION, double, BEZIER_CONTROLLER_EP2)
    KRATOS_DEFINE_APPLICATION_VARIABLE(PARTICLE_MECHANICS_APPLICATION, double, BEZIER_CONTROLLER_EP3)
    KRATOS_DEFINE_APPLICATION_VARIABLE(PARTICLE_MECHANICS_APPLICATION, double, BEZIER_CONTROLLER_EP4)
    KRATOS_DEFINE_APPLICATION_VARIABLE(PARTICLE_MECHANICS_APPLICATION, double, SHEAR_COMPRESSION_REDUCTOR)
    KRATOS_DEFINE_APPLICATION_VARIABLE(PARTICLE_MECHANICS_APPLICATION, double, STRAIN_RATE_FACTOR_C1_TENSION)
    KRATOS_DEFINE_APPLICATION_VARIABLE(PARTICLE_MECHANICS_APPLICATION, double, STRAIN_RATE_FACTOR_C2_TENSION)
    KRATOS_DEFINE_APPLICATION_VARIABLE(PARTICLE_MECHANICS_APPLICATION, double, STRAIN_RATE_FACTOR_C1_COMPRESSION)
    KRATOS_DEFINE_APPLICATION_VARIABLE(PARTICLE_MECHANICS_APPLICATION, double, STRAIN_RATE_FACTOR_C2_COMPRESSION)
    KRATOS_DEFINE_APPLICATION_VARIABLE(PARTICLE_MECHANICS_APPLICATION, double, STRAIN_RATE_FACTOR_C1_YOUNGS_MOD)
    KRATOS_DEFINE_APPLICATION_VARIABLE(PARTICLE_MECHANICS_APPLICATION, double, STRAIN_RATE_FACTOR_C2_YOUNGS_MOD)
    KRATOS_DEFINE_APPLICATION_VARIABLE(PARTICLE_MECHANICS_APPLICATION, int, TENSION_YIELD_MODEL)

    // Mesh variables
    KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, std::vector<typename Geometry<Node<3>>::Pointer>, GEOMETRY_NEIGHBOURS)

    /* NODAL VARIABLES */
    // Conditions
    // Particle Conditions
    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(PARTICLE_MECHANICS_APPLICATION, MPC_COORD )
    KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, int, MPC_CONDITION_ID )
    KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, bool, MPC_IS_NEUMANN )
    KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, double, MPC_AREA )
    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(PARTICLE_MECHANICS_APPLICATION, MPC_NORMAL )
    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(PARTICLE_MECHANICS_APPLICATION, MPC_DISPLACEMENT )
    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(PARTICLE_MECHANICS_APPLICATION, MPC_IMPOSED_DISPLACEMENT )
    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(PARTICLE_MECHANICS_APPLICATION, MPC_VELOCITY )
    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(PARTICLE_MECHANICS_APPLICATION, MPC_IMPOSED_VELOCITY )
    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(PARTICLE_MECHANICS_APPLICATION, MPC_ACCELERATION )
    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(PARTICLE_MECHANICS_APPLICATION, MPC_IMPOSED_ACCELERATION )
    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(PARTICLE_MECHANICS_APPLICATION, MPC_CONTACT_FORCE )
    KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, int, PARTICLES_PER_CONDITION )

    // Essential Boundary Conditions
    KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, double, PENALTY_FACTOR )
    KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, int, MPC_BOUNDARY_CONDITION_TYPE )

    // Natural Boundary Conditions
    // Nodal load variables
    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( PARTICLE_MECHANICS_APPLICATION, POINT_LOAD )
    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( PARTICLE_MECHANICS_APPLICATION, LINE_LOAD )
    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( PARTICLE_MECHANICS_APPLICATION, SURFACE_LOAD )

    // Momentum
    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(PARTICLE_MECHANICS_APPLICATION, NODAL_MOMENTUM )
    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(PARTICLE_MECHANICS_APPLICATION, NODAL_INERTIA )
    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(PARTICLE_MECHANICS_APPLICATION, NODAL_INTERNAL_FORCE )

    // Solver related variables
    KRATOS_DEFINE_APPLICATION_VARIABLE(PARTICLE_MECHANICS_APPLICATION, bool, IGNORE_GEOMETRIC_STIFFNESS)
    KRATOS_DEFINE_APPLICATION_VARIABLE(PARTICLE_MECHANICS_APPLICATION, bool, IS_AXISYMMETRIC)
    KRATOS_DEFINE_APPLICATION_VARIABLE(PARTICLE_MECHANICS_APPLICATION, bool, IS_PQMPM)
    KRATOS_DEFINE_APPLICATION_VARIABLE(PARTICLE_MECHANICS_APPLICATION, bool, IS_MAKE_NORMAL_MP_IF_PQMPM_FAILS)
    KRATOS_DEFINE_APPLICATION_VARIABLE(PARTICLE_MECHANICS_APPLICATION, double, PQMPM_SUBPOINT_MIN_VOLUME_FRACTION)

    // Explicit time integration variables
    KRATOS_DEFINE_APPLICATION_VARIABLE(PARTICLE_MECHANICS_APPLICATION, bool, CALCULATE_MUSL_VELOCITY_FIELD)
    KRATOS_DEFINE_APPLICATION_VARIABLE(PARTICLE_MECHANICS_APPLICATION, bool, IS_EXPLICIT)
    KRATOS_DEFINE_APPLICATION_VARIABLE(PARTICLE_MECHANICS_APPLICATION, bool, IS_EXPLICIT_CENTRAL_DIFFERENCE)
    KRATOS_DEFINE_APPLICATION_VARIABLE(PARTICLE_MECHANICS_APPLICATION, int, EXPLICIT_STRESS_UPDATE_OPTION)
    KRATOS_DEFINE_APPLICATION_VARIABLE(PARTICLE_MECHANICS_APPLICATION, bool, CALCULATE_EXPLICIT_MP_STRESS)
    KRATOS_DEFINE_APPLICATION_VARIABLE(PARTICLE_MECHANICS_APPLICATION, bool, EXPLICIT_MAP_GRID_TO_MP)
    KRATOS_DEFINE_APPLICATION_VARIABLE(PARTICLE_MECHANICS_APPLICATION, bool, IS_FIX_EXPLICIT_MP_ON_GRID_EDGE)
    KRATOS_DEFINE_APPLICATION_VARIABLE(PARTICLE_MECHANICS_APPLICATION, bool, EXPLICIT_CONTACT_RELEASE)
    KRATOS_DEFINE_APPLICATION_VARIABLE(PARTICLE_MECHANICS_APPLICATION, std::string, EXPLICIT_CONTACT_RELEASE_MODEL_PART)

    // Cosim coupling variables
    KRATOS_DEFINE_APPLICATION_VARIABLE(PARTICLE_MECHANICS_APPLICATION, int, INTERFACE_EQUATION_ID)
    KRATOS_DEFINE_APPLICATION_VARIABLE(PARTICLE_MECHANICS_APPLICATION, bool, IS_COSIM_COUPLED)
}

#endif // KRATOS_PARTICLE_MECHANICS_APPLICATION_VARIABLES_H_INCLUDED  defined