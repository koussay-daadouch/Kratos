#include "compressible_potential_flow_application_variables.h"

namespace Kratos
{
// Degrees of freedom
KRATOS_CREATE_VARIABLE(double, VELOCITY_POTENTIAL)
KRATOS_CREATE_VARIABLE(double, AUXILIARY_VELOCITY_POTENTIAL)
KRATOS_CREATE_VARIABLE(double, PSI)

// Lagrange multipliers
KRATOS_CREATE_VARIABLE(double, LAGRANGE_MULTIPLIER_0)
KRATOS_CREATE_VARIABLE(double, LAGRANGE_MULTIPLIER_1)
KRATOS_CREATE_VARIABLE(double, LAGRANGE_MULTIPLIER_2)
KRATOS_CREATE_VARIABLE(double, LAGRANGE_MULTIPLIER_3)
KRATOS_CREATE_VARIABLE(double, LAGRANGE_MULTIPLIER_4)
KRATOS_CREATE_VARIABLE(double, LAGRANGE_MULTIPLIER_5)
KRATOS_CREATE_VARIABLE(double, LAGRANGE_MULTIPLIER_6)
KRATOS_CREATE_VARIABLE(double, LAGRANGE_MULTIPLIER_7)
KRATOS_CREATE_VARIABLE(double, LAGRANGE_MULTIPLIER_8)
KRATOS_CREATE_VARIABLE(double, LAGRANGE_MULTIPLIER_9)

//Embedded variables
KRATOS_CREATE_VARIABLE(double, GEOMETRY_DISTANCE)
KRATOS_CREATE_VARIABLE(double, ROTATION_ANGLE);

// Wake variables
KRATOS_CREATE_VARIABLE(double, WAKE_DISTANCE);
KRATOS_CREATE_VARIABLE(Vector, WAKE_ELEMENTAL_DISTANCES);
KRATOS_CREATE_VARIABLE(Vector, WAKE_ORIGIN);
KRATOS_CREATE_VARIABLE(int, TE_ELEMENT_COUNTER);

// Adjoint variables
KRATOS_CREATE_VARIABLE(double, ADJOINT_VELOCITY_POTENTIAL)
KRATOS_CREATE_VARIABLE(double, ADJOINT_AUXILIARY_VELOCITY_POTENTIAL)

// Flow field magnitudes
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(VELOCITY_LOWER)
KRATOS_CREATE_VARIABLE(double, PRESSURE_LOWER)
KRATOS_CREATE_VARIABLE(double, POTENTIAL_JUMP)
KRATOS_CREATE_VARIABLE(double, ENERGY_NORM_REFERENCE)
KRATOS_CREATE_VARIABLE(double, POTENTIAL_ENERGY_REFERENCE)
KRATOS_CREATE_VARIABLE(double, HEAT_CAPACITY_RATIO)

// Free stream magnitudes
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(FREE_STREAM_VELOCITY)
KRATOS_CREATE_VARIABLE(double, FREE_STREAM_DENSITY)
KRATOS_CREATE_VARIABLE(double, FREE_STREAM_MACH)

// Integral magnitudes
KRATOS_CREATE_VARIABLE(double, LIFT_COEFFICIENT)
KRATOS_CREATE_VARIABLE(double, DRAG_COEFFICIENT)
KRATOS_CREATE_VARIABLE(double, MOMENT_COEFFICIENT)
KRATOS_CREATE_VARIABLE(double, LIFT_COEFFICIENT_JUMP)
KRATOS_CREATE_VARIABLE(double, LIFT_COEFFICIENT_FAR_FIELD)
KRATOS_CREATE_VARIABLE(double, DRAG_COEFFICIENT_FAR_FIELD)

// Geometrical variables
KRATOS_CREATE_VARIABLE(double, REFERENCE_CHORD)
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(WAKE_NORMAL)
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(WING_SPAN_DIRECTION)

// Markers
KRATOS_CREATE_VARIABLE(int, WAKE)
KRATOS_CREATE_VARIABLE(int, KUTTA)
KRATOS_CREATE_VARIABLE(int, WING_TIP)
KRATOS_CREATE_VARIABLE(bool, TRAILING_EDGE)
KRATOS_CREATE_VARIABLE(bool, UPPER_SURFACE)
KRATOS_CREATE_VARIABLE(bool, LOWER_SURFACE)
KRATOS_CREATE_VARIABLE(bool, UPPER_WAKE)
KRATOS_CREATE_VARIABLE(bool, LOWER_WAKE)
KRATOS_CREATE_VARIABLE(int, AIRFOIL)

// To be removed
KRATOS_CREATE_VARIABLE(int, TRAILING_EDGE_ELEMENT)
KRATOS_CREATE_VARIABLE(int, DECOUPLED_TRAILING_EDGE_ELEMENT)
KRATOS_CREATE_VARIABLE(int, DEACTIVATED_WAKE)
KRATOS_CREATE_VARIABLE(int, ALL_TRAILING_EDGE)
KRATOS_CREATE_VARIABLE(int, ZERO_VELOCITY_CONDITION)
} // namespace Kratos
