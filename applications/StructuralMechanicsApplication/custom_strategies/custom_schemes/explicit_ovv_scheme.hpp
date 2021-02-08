//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license:
//kratos/license.txt
//
//  Main authors:    Ignasi de Pouplana
//
//

#if !defined(KRATOS_EXPLICIT_OVV_SCHEME_HPP_INCLUDED)
#define KRATOS_EXPLICIT_OVV_SCHEME_HPP_INCLUDED

/* System includes */
#include <fstream>

/* External includes */

/* Project includes */
#include "custom_strategies/custom_schemes/explicit_vv_scheme.hpp"
#include "utilities/variable_utils.h"
#include "custom_utilities/explicit_integration_utilities.h"

namespace Kratos {

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}

///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @class ExplicitOVVScheme
 * @ingroup StructuralMechanicsApplciation
 * @brief An explicit forward euler scheme with a split of the inertial term
 * @author Ignasi de Pouplana
 */
template <class TSparseSpace,
          class TDenseSpace //= DenseSpace<double>
          >
class ExplicitOVVScheme
    : public ExplicitVVScheme<TSparseSpace, TDenseSpace> {

public:
    ///@name Type Definitions
    ///@{

    /// The definition of the base type
    typedef ExplicitVVScheme<TSparseSpace, TDenseSpace> BaseType;

    /// Some definitions related with the base class
    typedef typename BaseType::DofsArrayType DofsArrayType;
    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;
    typedef typename BaseType::TSystemVectorType TSystemVectorType;
    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    /// The arrays of elements and nodes
    typedef ModelPart::ElementsContainerType ElementsArrayType;
    typedef ModelPart::NodesContainerType NodesArrayType;

    /// Definition of the size type
    typedef std::size_t SizeType;

    /// Definition of the index type
    typedef std::size_t IndexType;

    /// Definition fo the node iterator
    typedef typename ModelPart::NodeIterator NodeIterator;

    /// The definition of the numerical limit
    static constexpr double numerical_limit = std::numeric_limits<double>::epsilon();

    using BaseType::mDeltaTime;
    using BaseType::mAlpha;
    using BaseType::mBeta;
    using BaseType::mTheta1;

    /// Counted pointer of ExplicitOVVScheme
    KRATOS_CLASS_POINTER_DEFINITION(ExplicitOVVScheme);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor.
     * @details The ExplicitOVVScheme method
     */
    ExplicitOVVScheme()
        : ExplicitVVScheme<TSparseSpace, TDenseSpace>() {}

    /** Destructor.
    */
    virtual ~ExplicitOVVScheme() {}

    ///@}
    ///@name Operators
    ///@{

    /**
     * @brief This is the place to initialize the Scheme. This is intended to be called just once when the strategy is initialized
     * @param rModelPart The model of the problem to solve
     */
    void Initialize(ModelPart& rModelPart) override
    {
        KRATOS_TRY

        const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

        mDelta = r_current_process_info[LOAD_FACTOR];

        BaseType::Initialize(rModelPart);

        KRATOS_CATCH("")
    }

    /**
     * @brief This method updates the translation DoF
     * @param itCurrentNode The iterator of the current node
     * @param DisplacementPosition The position of the displacement dof on the database
     * @param DomainSize The current dimention of the problem
     */
    void PredictTranslationalDegreesOfFreedom(
        NodeIterator itCurrentNode,
        const IndexType DisplacementPosition,
        const SizeType DomainSize = 3
        ) override
    {
        const double nodal_mass = itCurrentNode->GetValue(NODAL_MASS);
        const array_1d<double, 3>& r_external_forces = itCurrentNode->FastGetSolutionStepValue(FORCE_RESIDUAL);
        const array_1d<double, 3>& r_previous_external_forces = itCurrentNode->FastGetSolutionStepValue(FORCE_RESIDUAL,1);
        const array_1d<double, 3>& r_current_internal_force = itCurrentNode->FastGetSolutionStepValue(NODAL_INERTIA);
        const array_1d<double, 3>& r_previous_internal_force = itCurrentNode->FastGetSolutionStepValue(NODAL_INERTIA,1);
        const array_1d<double, 3>& r_current_damping_force = itCurrentNode->FastGetSolutionStepValue(MIDDLE_VELOCITY);
        array_1d<double, 3>& r_current_displacement = itCurrentNode->FastGetSolutionStepValue(DISPLACEMENT);
        array_1d<double, 3>& r_current_velocity = itCurrentNode->FastGetSolutionStepValue(VELOCITY);

        std::array<bool, 3> fix_displacements = {false, false, false};
        fix_displacements[0] = (itCurrentNode->GetDof(DISPLACEMENT_X, DisplacementPosition).IsFixed());
        fix_displacements[1] = (itCurrentNode->GetDof(DISPLACEMENT_Y, DisplacementPosition + 1).IsFixed());
        if (DomainSize == 3)
            fix_displacements[2] = (itCurrentNode->GetDof(DISPLACEMENT_Z, DisplacementPosition + 2).IsFixed());

        // Solution of the explicit equation:
        if ((nodal_mass*mDelta*mDelta) > numerical_limit){
            for (IndexType j = 0; j < DomainSize; j++) {
                if (fix_displacements[j] == false) {
                    r_current_displacement[j] += r_current_velocity[j]*mDeltaTime + 0.5 * (mTheta1*r_external_forces[j]+(1.0-mTheta1)*r_previous_external_forces[j]
                                                                                           - (mTheta1*r_current_internal_force[j]+(1.0-mTheta1)*r_previous_internal_force[j])
                                                                                           - r_current_damping_force[j])/(nodal_mass*mDelta*mDelta) * mDeltaTime * mDeltaTime;
                    r_current_velocity[j] += 0.5 * mDeltaTime * (mTheta1*r_external_forces[j]+(1.0-mTheta1)*r_previous_external_forces[j]
                                                                 - (mTheta1*r_current_internal_force[j]+(1.0-mTheta1)*r_previous_internal_force[j])
                                                                 - r_current_damping_force[j])/(nodal_mass*mDelta*mDelta);
                }
            }
        }
        else {
            noalias(r_current_displacement) = ZeroVector(3);
            noalias(r_current_velocity) = ZeroVector(3);
        }
    }

    /**
     * @brief This method updates the translation DoF
     * @param itCurrentNode The iterator of the current node
     * @param DisplacementPosition The position of the displacement dof on the database
     * @param DomainSize The current dimention of the problem
     */
    void UpdateTranslationalDegreesOfFreedom(
        NodeIterator itCurrentNode,
        const IndexType DisplacementPosition,
        const SizeType DomainSize = 3
        ) override
    {
        const double nodal_mass = itCurrentNode->GetValue(NODAL_MASS);
        const array_1d<double, 3>& r_external_forces = itCurrentNode->FastGetSolutionStepValue(FORCE_RESIDUAL);
        const array_1d<double, 3>& r_previous_external_forces = itCurrentNode->FastGetSolutionStepValue(FORCE_RESIDUAL,1);
        const array_1d<double, 3>& r_current_internal_force = itCurrentNode->FastGetSolutionStepValue(NODAL_INERTIA);
        const array_1d<double, 3>& r_previous_internal_force = itCurrentNode->FastGetSolutionStepValue(NODAL_INERTIA,1);
        const array_1d<double, 3>& r_current_damping_force = itCurrentNode->FastGetSolutionStepValue(MIDDLE_VELOCITY);
        array_1d<double, 3>& r_current_velocity = itCurrentNode->FastGetSolutionStepValue(VELOCITY);

        std::array<bool, 3> fix_displacements = {false, false, false};
        fix_displacements[0] = (itCurrentNode->GetDof(DISPLACEMENT_X, DisplacementPosition).IsFixed());
        fix_displacements[1] = (itCurrentNode->GetDof(DISPLACEMENT_Y, DisplacementPosition + 1).IsFixed());
        if (DomainSize == 3)
            fix_displacements[2] = (itCurrentNode->GetDof(DISPLACEMENT_Z, DisplacementPosition + 2).IsFixed());

        // Solution of the explicit equation:
        if ((nodal_mass*mDelta*mDelta) > numerical_limit){
            for (IndexType j = 0; j < DomainSize; j++) {
                if (fix_displacements[j] == false) {
                    r_current_velocity[j] += 0.5 * mDeltaTime * (mTheta1*r_external_forces[j]+(1.0-mTheta1)*r_previous_external_forces[j]
                                                                 - (mTheta1*r_current_internal_force[j]+(1.0-mTheta1)*r_previous_internal_force[j])
                                                                 - r_current_damping_force[j])/(nodal_mass*mDelta*mDelta);
                }
            }
        }
        else {
            noalias(r_current_velocity) = ZeroVector(3);
        }

        const array_1d<double, 3>& r_previous_velocity = itCurrentNode->FastGetSolutionStepValue(VELOCITY,1);
        array_1d<double, 3>& r_current_acceleration = itCurrentNode->FastGetSolutionStepValue(ACCELERATION);

        noalias(r_current_acceleration) = (1.0/mDeltaTime) * (r_current_velocity - r_previous_velocity);
    }

    ///@}
    ///@name Operations
    ///@{

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:

    ///@}
    ///@name Protected Structs
    ///@{

    /**
     * @brief This struct contains the information related with the increment od time step
     */
    // struct DeltaTimeParameters {
    //     double Maximum;         // Maximum delta time
    //     double Fraction;        // Fraction of the delta time
    // };

    /**
     * @brief This struct contains the details of the time variables
     */
    // struct TimeVariables {
    //     double Current;        // n+1

    //     double Delta;          // Time step
    // };

    ///@name Protected static Member Variables
    ///@{

    // TimeVariables mTime;            /// This struct contains the details of the time variables
    // DeltaTimeParameters mDeltaTime; /// This struct contains the information related with the increment od time step

    ///@}
    ///@name Protected member Variables
    ///@{

    double mDelta;

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

}; /* Class ExplicitOVVScheme */

///@}

///@name Type Definitions
///@{

///@}

} /* namespace Kratos.*/

#endif /* KRATOS_EXPLICIT_OVV_SCHEME_HPP_INCLUDED  defined */
