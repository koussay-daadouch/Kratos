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

#if !defined(KRATOS_K_EPSILON_VMS_MONOLITHIC_K_BASED_EPSILON_K_BASED_WALL_CONDITION_DATA_H)
#define KRATOS_K_EPSILON_VMS_MONOLITHIC_K_BASED_EPSILON_K_BASED_WALL_CONDITION_DATA_H

// System includes

// External includes

// Project includes
#include "includes/node.h"
#include "includes/process_info.h"
#include "geometries/geometry.h"
#include "geometries/geometry_data.h"

// Application includes
#include "custom_conditions/vms_monolithic_k_based_wall_condition_derivatives.h"
#include "custom_conditions/vms_monolithic_k_based_wall_condition_derivative_utilities.h"

namespace Kratos
{
///@name Kratos Classes
///@{

namespace KEpsilonWallConditionData
{

template <unsigned int TDim, unsigned int TNumNodes>
class VMSMonolithicKBasedEpsilonKBasedWallConditionData
{
public:
    ///@name Type Definitions
    ///@{

    using IndexType = std::size_t;

    using NodeType = Node<3>;

    using GeometryType = Geometry<NodeType>;

    ///@}
    ///@name Static Operations
    ///@{

    static void Check(
        const Condition& rCondition,
        const ProcessInfo& rProcessInfo);

    static std::vector<const Variable<double>*> GetDofVariablesList();

    ///@}
    ///@name Classes
    ///@{

    class Fluid
    {
    public:
        ///@name Type Definitions
        ///@{

        using TResidualsDerivatives = VMSMonolithicKBasedWallConditionDerivatives<TDim, TNumNodes>;

        ///@}
        ///@name Classes
        ///@{

        class Primal
        {
        public:
            ///@name Type Definitions
            ///@{

            using Data = typename TResidualsDerivatives::Data;

            using ResidualsContributions = typename TResidualsDerivatives::ResidualsContributions;

            ///@}
        };

        class StateDerivatives
        {
        public:
            ///@name Classes
            ///@{

            class FirstDerivatives
            {
            public:
                ///@name Type Definitions
                ///@{

                using Data = typename TResidualsDerivatives::Data;

                using Velocity = typename TResidualsDerivatives::template VariableDerivatives<typename VMSMonolithicKBasedWallConditionDerivativeUtilities<TDim>::VelocityDerivative>;

                using TurbulenceModelVariable1 = typename TResidualsDerivatives::template VariableDerivatives<typename VMSMonolithicKBasedWallConditionDerivativeUtilities<TDim>::KDerivative>;

                using TurbulenceModelVariable2 = typename TResidualsDerivatives::template VariableDerivatives<typename VMSMonolithicKBasedWallConditionDerivativeUtilities<TDim>::NonRelatedDerivative>;

                ///@}
            };

            ///@}
        };

        // class SensitivityDerivatives
        // {
        // public:
        //     ///@name Type Definitions
        //     ///@{

        //     using Data = typename TResidualsDerivatives::Data;

        //     using Shape = typename TResidualsDerivatives::template VariableDerivatives<typename QSVMSDerivativeUtilities<TDim>::template ShapeDerivative<TNumNodes>>;

        //     ///@}
        // };

    ///@}
    };

    // class TurbulenceModelEquation1
    // {
    // public:
    //     ///@name Classes
    //     ///@{

    //     using EquationDataType = KEpsilonWallConditionData::KConditionDataDerivatives<TDim, TNumNodes>;

    //     using TResidualsDerivatives = ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<TDim, TNumNodes, typename EquationDataType::Data>;

    //     class Primal
    //     {
    //     public:
    //         ///@name Type Definitions
    //         ///@{

    //         using Data = typename TResidualsDerivatives::Data;

    //         using ResidualsContributions = typename TResidualsDerivatives::ResidualsContributions;

    //         ///@}
    //     };

    //     class StateDerivatives
    //     {
    //     public:
    //         ///@name Classes
    //         ///@{

    //         class FirstDerivatives
    //         {
    //         public:
    //             ///@name Type Definitions
    //             ///@{

    //             using Data = typename TResidualsDerivatives::Data;

    //             using Velocity = typename TResidualsDerivatives::template VariableDerivatives<typename EquationDataType::UDerivative>;

    //             using TurbulenceModelVariable1 = typename TResidualsDerivatives::template VariableDerivatives<typename EquationDataType::KDerivative>;

    //             using TurbulenceModelVariable2 = typename TResidualsDerivatives::template VariableDerivatives<typename EquationDataType::EpsilonDerivative>;

    //             ///@}
    //         };

    //         class SecondDerivatives
    //         {
    //         public:
    //             ///@name Type Definitions
    //             ///@{

    //             using Data = typename TResidualsDerivatives::Data;

    //             using TurbulenceModelVariableRate1 = typename TResidualsDerivatives::SecondDerivatives;

    //             ///@}

    //         };

    //         ///@}
    //     };

    //     class SensitivityDerivatives
    //     {
    //     public:
    //         ///@name Type Definitions
    //         ///@{

    //         using Data = typename TResidualsDerivatives::Data;

    //         using Shape = typename TResidualsDerivatives::template VariableDerivatives<typename EquationDataType::ShapeDerivative>;

    //         ///@}
    //     };

    //     ///@}
    // };

    // class TurbulenceModelEquation2
    // {
    // public:
    //     ///@name Classes
    //     ///@{

    //     using EquationDataType = KEpsilonWallConditionData::EpsilonConditionDataDerivatives<TDim, TNumNodes>;

    //     using TResidualsDerivatives = ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<TDim, TNumNodes, typename EquationDataType::Data>;

    //     class Primal
    //     {
    //     public:
    //         ///@name Type Definitions
    //         ///@{

    //         using Data = typename TResidualsDerivatives::Data;

    //         using ResidualsContributions = typename TResidualsDerivatives::ResidualsContributions;

    //         ///@}
    //     };

    //     class StateDerivatives
    //     {
    //     public:
    //         ///@name Classes
    //         ///@{

    //         class FirstDerivatives
    //         {
    //         public:
    //             ///@name Type Definitions
    //             ///@{

    //             using Data = typename TResidualsDerivatives::Data;

    //             using Velocity = typename TResidualsDerivatives::template VariableDerivatives<typename EquationDataType::UDerivative>;

    //             using TurbulenceModelVariable1 = typename TResidualsDerivatives::template VariableDerivatives<typename EquationDataType::KDerivative>;

    //             using TurbulenceModelVariable2 = typename TResidualsDerivatives::template VariableDerivatives<typename EquationDataType::EpsilonDerivative>;

    //             ///@}
    //         };

    //         class SecondDerivatives
    //         {
    //         public:
    //             ///@name Type Definitions
    //             ///@{

    //             using Data = typename TResidualsDerivatives::Data;

    //             using TurbulenceModelVariableRate2 = typename TResidualsDerivatives::SecondDerivatives;

    //             ///@}
    //         };

    //         ///@}
    //     };

    //     class SensitivityDerivatives
    //     {
    //     public:
    //         ///@name Type Definitions
    //         ///@{

    //         using Data = typename TResidualsDerivatives::Data;

    //         using Shape = typename TResidualsDerivatives::template VariableDerivatives<typename EquationDataType::ShapeDerivative>;

    //         ///@}
    //     };

    //     ///@}
    // };

    ///@}
};
} // namespace KEpsilonWallConditionData
} // namespace Kratos

#endif // KRATOS_K_EPSILON_VMS_MONOLITHIC_K_BASED_EPSILON_K_BASED_WALL_CONDITION_DATA_H