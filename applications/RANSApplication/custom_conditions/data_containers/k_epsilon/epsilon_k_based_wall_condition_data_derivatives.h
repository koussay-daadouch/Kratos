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

#if !defined(KRATOS_EPSILON_EPSILON_K_BASED_WALL_CONDITION_DATA_DERIVATIVES_H)
#define KRATOS_EPSILON_EPSILON_K_BASED_WALL_CONDITION_DATA_DERIVATIVES_H

// System includes

// External includes

// Project includes
#include "containers/variable.h"
#include "geometries/geometry.h"
#include "includes/node.h"
#include "includes/properties.h"
#include "includes/ublas_interface.h"

// Application includes
#include "custom_conditions/data_containers/k_epsilon/epsilon_k_based_wall_condition_data.h"

namespace Kratos
{
namespace KEpsilonWallConditionData
{
///@name Kratos Classes
///@{

template <unsigned int TDim>
class EpsilonKBasedWallConditionDataDerivatives
{
public:
    ///@name Type Definitions
    ///@{

    using NodeType = Node<3>;

    using GeometryType = Geometry<NodeType>;

    using IndexType = std::size_t;

    ///@}
    ///@name Classes
    ///@{

    // forward declaration of the data gauss point data holder
    class Data;

    class UDerivative
    {
    public:
        /// name@ Type Definitions
        ///@{

        using DataType = Data;

        static constexpr unsigned int TDerivativeDimension = TDim;

        static constexpr double TSelfWeight = 0.0;

        ///@}
        ///@name Life Cycle
        ///@{

        UDerivative(const DataType& rData) : mrData(rData) {}

        ///@}
        ///@name Operations
        ///@{

        double CalculateWallFluxDerivative(
            const IndexType ConditionNodeIndex,
            const IndexType DirectionIndex,
            const double W,
            const Vector& rN,
            const double WDerivative,
            const double DetJDerivative,
            const IndexType ParentElementNodeIndex,
            const double ParentElementW,
            const Vector& rParentElementN,
            const Matrix& rParentElementdNdX,
            const double ParentElementWDerivative,
            const double ParentElementDetJDerivative,
            const Matrix& ParentElementdNdXDerivative) const;

        double CalculateWallFluxDerivative(
            const IndexType DirectionIndex,
            const double W,
            const Vector& rN,
            const IndexType ParentElementNodeIndex,
            const double ParentElementW,
            const Vector& rParentElementN,
            const Matrix& rParentElementdNdX,
            const double ParentElementWDerivative,
            const double ParentElementDetJDerivative,
            const Matrix& ParentElementdNdXDerivative) const;

        ///@}

    private:
        ///@name Private Members
        ///@{

        const DataType& mrData;

        ///@}
    };

    class KDerivative
    {
    public:
        /// name@ Type Definitions
        ///@{

        using DataType = Data;

        static constexpr unsigned int TDerivativeDimension = 1;

        static constexpr double TSelfWeight = 0.0;

        ///@}
        ///@name Life Cycle
        ///@{

        KDerivative(const DataType& rData) : mrData(rData) {}

        ///@}
        ///@name Operations
        ///@{

        double CalculateWallFluxDerivative(
            const IndexType ConditionNodeIndex,
            const IndexType DirectionIndex,
            const double W,
            const Vector& rN,
            const double WDerivative,
            const double DetJDerivative,
            const IndexType ParentElementNodeIndex,
            const double ParentElementW,
            const Vector& rParentElementN,
            const Matrix& rParentElementdNdX,
            const double ParentElementWDerivative,
            const double ParentElementDetJDerivative,
            const Matrix& ParentElementdNdXDerivative) const;

        double CalculateWallFluxDerivative(
            const IndexType DirectionIndex,
            const double W,
            const Vector& rN,
            const IndexType ParentElementNodeIndex,
            const double ParentElementW,
            const Vector& rParentElementN,
            const Matrix& rParentElementdNdX,
            const double ParentElementWDerivative,
            const double ParentElementDetJDerivative,
            const Matrix& ParentElementdNdXDerivative) const;

        ///@}

    private:
        ///@name Private Members
        ///@{

        const DataType& mrData;

        ///@}
    };

    class EpsilonDerivative
    {
    public:
        /// name@ Type Definitions
        ///@{

        using DataType = Data;

        static constexpr unsigned int TDerivativeDimension = 1;

        static constexpr double TSelfWeight = 1.0;

        ///@}
        ///@name Life Cycle
        ///@{

        EpsilonDerivative(const DataType& rData) : mrData(rData) {}

        ///@}
        ///@name Operations
        ///@{

        double CalculateWallFluxDerivative(
            const IndexType ConditionNodeIndex,
            const IndexType DirectionIndex,
            const double W,
            const Vector& rN,
            const double WDerivative,
            const double DetJDerivative,
            const IndexType ParentElementNodeIndex,
            const double ParentElementW,
            const Vector& rParentElementN,
            const Matrix& rParentElementdNdX,
            const double ParentElementWDerivative,
            const double ParentElementDetJDerivative,
            const Matrix& ParentElementdNdXDerivative) const;

        double CalculateWallFluxDerivative(
            const IndexType DirectionIndex,
            const double W,
            const Vector& rN,
            const IndexType ParentElementNodeIndex,
            const double ParentElementW,
            const Vector& rParentElementN,
            const Matrix& rParentElementdNdX,
            const double ParentElementWDerivative,
            const double ParentElementDetJDerivative,
            const Matrix& ParentElementdNdXDerivative) const;

        ///@}

    private:
        ///@name Private Members
        ///@{

        const DataType& mrData;

        ///@}
    };

    class ShapeDerivative
    {
    public:
        /// name@ Type Definitions
        ///@{

        using DataType = Data;

        static constexpr unsigned int TDerivativeDimension = TDim;

        static constexpr double TSelfWeight = 0.0;

        ///@}
        ///@name Life Cycle
        ///@{

        ShapeDerivative(const DataType& rData) : mrData(rData) {}

        ///@}
        ///@name Operations
        ///@{

        double CalculateWallFluxDerivative(
            const IndexType ConditionNodeIndex,
            const IndexType DirectionIndex,
            const double W,
            const Vector& rN,
            const double WDerivative,
            const double DetJDerivative,
            const IndexType ParentElementNodeIndex,
            const double ParentElementW,
            const Vector& rParentElementN,
            const Matrix& rParentElementdNdX,
            const double ParentElementWDerivative,
            const double ParentElementDetJDerivative,
            const Matrix& ParentElementdNdXDerivative) const;

        double CalculateWallFluxDerivative(
            const IndexType DirectionIndex,
            const double W,
            const Vector& rN,
            const IndexType ParentElementNodeIndex,
            const double ParentElementW,
            const Vector& rParentElementN,
            const Matrix& rParentElementdNdX,
            const double ParentElementWDerivative,
            const double ParentElementDetJDerivative,
            const Matrix& ParentElementdNdXDerivative) const;

        ///@}

    private:
        ///@name Private Members
        ///@{

        const DataType& mrData;

        ///@}
    };



    class Data : public EpsilonKBasedWallConditionData
    {
    public:
        ///@name Type Definitions
        ///@{

        using BaseType = EpsilonKBasedWallConditionData;

        ///@}
        ///@name Life Cycle
        ///@{

        Data(
            const IndexType NumberOfGaussPoints,
            const GeometryType& rGeometry,
            const Properties& rProperties,
            const ProcessInfo& rProcessInfo);

        ///@}
        ///@name Static Operations
        ///@{

        using BaseType::GetScalarVariable;

        static const Variable<double>& GetAdjointScalarVariable();

        static void Check(
            const Condition& rCondition,
            const ProcessInfo& rCurrentProcessInfo);

        static const std::string GetName() { return "KEpsilonEpsilonKBasedWallConditionAdjointData"; }

        ///@}
        ///@name Operations
        ///@{

        void CalculateGaussPointData(
            const Vector& rN,
            const Vector& rParentElementN,
            const Matrix& rParentElementdNdX,
            const int Step = 0);

        double GetWallFlux() const { return mWallFlux; }

        ///@}

    protected:
        ///@name Protected Members
        ///@{

        const GeometryType& mrParentElementGeometry;
        const IndexType mNumberOfGaussPoints;

        using  BaseType::mEpsilonSigma;
        using  BaseType::mKappa;
        using  BaseType::mYPlus;
        using  BaseType::mCmu25;
        using  BaseType::mDensity;

        double mCmu;
        double mKinematicViscosity;
        double mUTau;
        double mTurbulentKineticEnergy;
        double mTurbulentKinematicViscosity;
        double mWallVelocityMagnitude;
        double mWallHeight;
        double mNormalMagnitude;
        double mBeta;
        double mYPlusLimit;
        double mWallFlux;

        array_1d<double, 3> mUnitNormal;
        array_1d<double, 3> mWallVelocity;

        ///@}
        ///@name Private Friends
        ///@{

        friend class UDerivative;
        friend class KDerivative;
        friend class EpsilonDerivative;
        friend class ShapeDerivative;

        ///@}
    };

    ///@}
};

///@}

} // namespace KEpsilonElementData

} // namespace Kratos

#endif // KRATOS_EPSILON_EPSILON_K_BASED_WALL_CONDITION_DATA_DERIVATIVES_H