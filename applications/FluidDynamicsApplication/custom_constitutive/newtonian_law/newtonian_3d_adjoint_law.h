//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

#if !defined (KRATOS_NEWTONIAN_ADJOINT_LAW_3D_H_INCLUDED)
#define  KRATOS_NEWTONIAN_ADJOINT_LAW_3D_H_INCLUDED

// System includes

// External includes

// Project includes

// Application includes
#include "custom_constitutive/fluid_adjoint_constitutive_law.h"


namespace Kratos
{

///@addtogroup FluidDynamicsApplication
///@{

///@name Kratos Classes
///@{

/// This class contains the common infrastructure for fluid constitutive laws.
class KRATOS_API(FLUID_DYNAMICS_APPLICATION) Newtonian3DAdjointLaw : public FluidAdjointConstitutiveLaw
{
public:

    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(Newtonian3DAdjointLaw);

    using BaseType = FluidAdjointConstitutiveLaw;

    using GeometryType = typename BaseType::GeometryType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    Newtonian3DAdjointLaw(ConstitutiveLaw& rConstitutiveLaw);

    /// Copy constructor.
    Newtonian3DAdjointLaw (const Newtonian3DAdjointLaw& rOther);

    /// Destructor
    virtual ~Newtonian3DAdjointLaw() override;

    ///@}
    ///@name Operations
    ///@{

    void CalculateMaterialResponseCauchyDerivative(
        ConstitutiveLaw::Parameters& rValuesDerivative,
        ConstitutiveLaw::Parameters& rValues,
        const IndexType NodeIndex,
        const Variable<double>& rDerivativeVariable,
        const double EffectiveViscosity,
        const double EffectiveViscosityDerivative) override;

    int Check(
        const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const ProcessInfo& rCurrentProcessInfo) override;

    double CalculateEffectiveViscosityDerivative(
        ConstitutiveLaw::Parameters& rValuesDerivative,
        ConstitutiveLaw::Parameters& rValues,
        const IndexType NodeIndex,
        const Variable<double>& rDerivativeVariable) override;

    ///@}
    ///@name Input and output
    ///@{

    /// @return A short string identifying this constitutive law instance.
    std::string Info() const override;

    /// Print basic information about this constitutive law instance.
    void PrintInfo(std::ostream& rOStream) const override;

    /// Print detailed information about this constitutive law instance and its managed data.
    void PrintData(std::ostream& rOStream) const override;

    ///@}

}; // Class Newtonian3DAdjointLaw

}  // namespace Kratos.

#endif // KRATOS_NEWTONIAN_ADJOINT_LAW_3D_H_INCLUDED  defined