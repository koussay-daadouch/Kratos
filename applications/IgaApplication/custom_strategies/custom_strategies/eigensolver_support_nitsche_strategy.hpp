//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//   License:        BSD License
//   Kratos default license: kratos/license.txt

#if !defined(KRATOS_EIGENSOLVER_SUPPORT_NITSCHE_STRATEGY )
#define  KRATOS_EIGENSOLVER_SUPPORT_NITSCHE_STRATEGY

// System includes

// External includes

// Project includes
#include "solving_strategies/strategies/solving_strategy.h"
#include "utilities/builtin_timer.h"
#include "spaces/ublas_space.h"

#include "containers/model.h"

// Application includes
#include "iga_application_variables.h"

namespace Kratos
{

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

/// Strategy for solving generalized eigenvalue problems.
template<class TSparseSpace,
         class TDenseSpace,
         class TLinearSolver
         >
class EigensolverSupportNitscheStrategy
    : public SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(EigensolverSupportNitscheStrategy);

    typedef SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

    typedef typename BaseType::TSchemeType::Pointer SchemePointerType;

    typedef typename BaseType::TBuilderAndSolverType::Pointer BuilderAndSolverPointerType;

    typedef typename TDenseSpace::VectorType DenseVectorType;

    typedef typename TDenseSpace::MatrixType DenseMatrixType;

    typedef TSparseSpace SparseSpaceType;

    typedef typename TSparseSpace::VectorPointerType SparseVectorPointerType;

    typedef typename TSparseSpace::MatrixPointerType SparseMatrixPointerType;

    typedef typename TSparseSpace::MatrixType SparseMatrixType;

    typedef typename TSparseSpace::VectorType SparseVectorType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    EigensolverSupportNitscheStrategy(
        ModelPart& rModelPart,
        SchemePointerType pScheme,
        BuilderAndSolverPointerType pBuilderAndSolver
        )
        : SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(rModelPart)
    {
        KRATOS_TRY

        mpScheme = pScheme;

        mpBuilderAndSolver = pBuilderAndSolver;

        // ensure initialization of system matrices in InitializeSolutionStep()
        mpBuilderAndSolver->SetDofSetIsInitializedFlag(false);

        // default echo level (mute)
        this->SetEchoLevel(0);

        // default rebuild level (build at each solution step)
        this->SetRebuildLevel(1);

        SparseMatrixType* AuxStiffnessMatrix = new SparseMatrixType;
        mpStiffnessMatrix = Kratos::shared_ptr<SparseMatrixType>(AuxStiffnessMatrix);
        SparseMatrixType* AuxQMatrix = new SparseMatrixType;
        mpQMatrix = Kratos::shared_ptr<SparseMatrixType>(AuxQMatrix);

        KRATOS_CATCH("")
    }

    /// Deleted copy constructor.
    EigensolverSupportNitscheStrategy(const EigensolverSupportNitscheStrategy& Other) = delete;

    /// Destructor.
    ~EigensolverSupportNitscheStrategy() override
    {
        // Clear() controls order of deallocation to avoid invalid memory access
        // in some special cases.
        // warning: BaseType::GetModelPart() may be invalid here.
        this->Clear();
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void SetIsInitialized(bool val)
    {
        mInitializeWasPerformed = val;
    }

    bool GetIsInitialized() const
    {
        return mInitializeWasPerformed;
    }

    void SetScheme(SchemePointerType pScheme)
    {
        mpScheme = pScheme;
    };

    SchemePointerType& pGetScheme()
    {
        return mpScheme;
    };

    void SetBuilderAndSolver(BuilderAndSolverPointerType pNewBuilderAndSolver)
    {
        mpBuilderAndSolver = pNewBuilderAndSolver;
    };

    BuilderAndSolverPointerType& pGetBuilderAndSolver()
    {
        return mpBuilderAndSolver;
    };

    SparseMatrixType& GetStiffnessMatrix()
    {
//         if (mpStiffnessMatrix == nullptr)
//         {
//             KRATOS_ERROR << "CHEK YOUR Stiffness MATRIX PLEASE IS INITIALIZED" << std::endl;
//         }
        return *mpStiffnessMatrix;
    }

    SparseMatrixType& GetQMatrix()
    {
//         if (mpQMatrix == nullptr)
//         {
//             KRATOS_ERROR << "CHEK YOUR Q MATRIX PLEASE IS INITIALIZED" << std::endl;
//         }
        return *mpQMatrix;
    }

    SparseMatrixPointerType& pGetStiffnessMatrix()
    {
//         if (mpStiffnessMatrix == nullptr)
//         {
//             KRATOS_ERROR << "CHEK YOUR Stiffness MATRIX PLEASE IS INITIALIZED" << std::endl;
//         }
        return mpStiffnessMatrix;
    }

    SparseMatrixPointerType& pGetQMatrix()
    {
//         if (mpQMatrix == nullptr)
//         {
//             KRATOS_ERROR << "CHEK YOUR Q MATRIX PLEASE IS INITIALIZED" << std::endl;
//         }
        return mpQMatrix;
    }

    void SetReformDofSetAtEachStepFlag(bool flag)
    {
        this->pGetBuilderAndSolver()->SetReshapeMatrixFlag(flag);
    }

    bool GetReformDofSetAtEachStepFlag() const
    {
        return this->pGetBuilderAndSolver()->GetReshapeMatrixFlag();
    }

    /// Set verbosity level of the solving strategy.
    /**
     * - 0 -> mute... no echo at all
     * - 1 -> print time and basic information
     * - 2 -> print linear solver data
     * - 3 -> print debug information
     */
    void SetEchoLevel(int Level) override
    {
        BaseType::SetEchoLevel(Level);
        this->pGetBuilderAndSolver()->SetEchoLevel(Level);
    }

    /**
     * Initialization to be performed once before using the strategy.
     */
    void Initialize() override
    {
        KRATOS_TRY

        ModelPart& rModelPart = BaseType::GetModelPart();
        const int rank = rModelPart.GetCommunicator().MyPID();

        KRATOS_INFO_IF("EigensolverSupportNitscheStrategy", BaseType::GetEchoLevel() > 2 && rank == 0)
            <<  "Entering Initialize" << std::endl;

        if (mInitializeWasPerformed == false)
        {
            SchemePointerType& pScheme = this->pGetScheme();

            if (pScheme->SchemeIsInitialized() == false)
                pScheme->Initialize(rModelPart);

            if (pScheme->ElementsAreInitialized() == false)
                pScheme->InitializeElements(rModelPart);

            if (pScheme->ConditionsAreInitialized() == false)
                pScheme->InitializeConditions(rModelPart);
        }

        KRATOS_INFO_IF("EigensolverSupportNitscheStrategy", BaseType::GetEchoLevel() > 2 && rank == 0)
            <<  "Exiting Initialize" << std::endl;

        KRATOS_CATCH("")
    }

    /**
     * Clears the internal storage
     */
    void Clear() override
    {
        KRATOS_TRY

        // if the preconditioner is saved between solves, it should be cleared here
        BuilderAndSolverPointerType& pBuilderAndSolver = this->pGetBuilderAndSolver();
        pBuilderAndSolver->GetLinearSystemSolver()->Clear();

        if (this->pGetStiffnessMatrix() != nullptr)
            this->pGetStiffnessMatrix() = nullptr;

        if (this->pGetQMatrix() != nullptr)
            this->pGetQMatrix() = nullptr;

        // Re-setting internal flag to ensure that the dof sets are recalculated
        pBuilderAndSolver->SetDofSetIsInitializedFlag(false);

        pBuilderAndSolver->Clear();

        this->pGetScheme()->Clear();

        mInitializeWasPerformed = false;

        KRATOS_CATCH("")
    }

    /**
     * Performs all the required operations that should be done (for each step)
     * before solving the solution step.
     * A member variable should be used as a flag to make sure this function is called only once per step.
     */
    void InitializeSolutionStep() override
    {
        KRATOS_TRY

        ModelPart& rModelPart = BaseType::GetModelPart();
        const int rank = rModelPart.GetCommunicator().MyPID();

        KRATOS_INFO_IF("EigensolverSupportNitscheStrategy", BaseType::GetEchoLevel() > 2 && rank == 0)
            <<  "Entering InitializeSolutionStep" << std::endl;

        BuilderAndSolverPointerType& pBuilderAndSolver = this->pGetBuilderAndSolver();
        SchemePointerType& pScheme = this->pGetScheme();
        SparseMatrixPointerType& pQMatrix = this->pGetQMatrix();
        SparseMatrixType& rQMatrix = this->GetQMatrix();

        // Initialize dummy vectors
        SparseVectorPointerType pDx = SparseSpaceType::CreateEmptyVectorPointer();
        SparseVectorPointerType pb = SparseSpaceType::CreateEmptyVectorPointer();
        auto& rDx = *pDx;
        auto& rb = *pb;

        // Reset solution dofs
        BuiltinTimer system_construction_time;
        if (pBuilderAndSolver->GetDofSetIsInitializedFlag() == false ||
            pBuilderAndSolver->GetReshapeMatrixFlag() == true)
        {
            // Set up list of dofs
            BuiltinTimer setup_dofs_time;
            pBuilderAndSolver->SetUpDofSet(pScheme, rModelPart);

            KRATOS_INFO_IF("Setup Dofs Time", BaseType::GetEchoLevel() > 0 && rank == 0)
                << setup_dofs_time.ElapsedSeconds() << std::endl;

            // Set global equation ids
            BuiltinTimer setup_system_time;
            pBuilderAndSolver->SetUpSystem(rModelPart);

            KRATOS_INFO_IF("Setup System Time", BaseType::GetEchoLevel() > 0 && rank == 0)
                << setup_system_time.ElapsedSeconds() << std::endl;

            // Resize and initialize system matrices
            BuiltinTimer system_matrix_resize_time;
            SparseMatrixPointerType& pStiffnessMatrix = this->pGetStiffnessMatrix();

            // Stiffness matrix
            pBuilderAndSolver->ResizeAndInitializeVectors(
                pScheme, pStiffnessMatrix, pDx, pb, rModelPart);

            // Q matrix
            pBuilderAndSolver->ResizeAndInitializeVectors(
                pScheme, pQMatrix, pDx, pb, rModelPart);

            KRATOS_INFO_IF("System Matrix Resize Time", BaseType::GetEchoLevel() > 0 && rank == 0)
                << system_matrix_resize_time.ElapsedSeconds() << std::endl;
        }
        else
        {
            SparseSpaceType::Resize(rb, SparseSpaceType::Size1(rQMatrix));
            SparseSpaceType::Set(rb, 0.0);
            SparseSpaceType::Resize(rDx, SparseSpaceType::Size1(rQMatrix));
            SparseSpaceType::Set(rDx, 0.0);
        }

        KRATOS_INFO_IF("System Construction Time", BaseType::GetEchoLevel() > 0 && rank == 0)
            << system_construction_time.ElapsedSeconds() << std::endl;

        // Initial operations ... things that are constant over the solution
        // step
        pBuilderAndSolver->InitializeSolutionStep(BaseType::GetModelPart(),
                                                  rQMatrix, rDx, rb);

        // Initial operations ... things that are constant over the solution
        // step
        pScheme->InitializeSolutionStep(BaseType::GetModelPart(), rQMatrix, rDx, rb);

        KRATOS_INFO_IF("EigensolverSupportNitscheStrategy", BaseType::GetEchoLevel() > 2 && rank == 0)
            <<  "Exiting InitializeSolutionStep" << std::endl;

        KRATOS_CATCH("")
    }

    bool SolveSolutionStep() override
    {
        KRATOS_TRY;

        ModelPart& rModelPart = BaseType::GetModelPart();

        SchemePointerType& pScheme = this->pGetScheme();
        SparseMatrixType& rStiffnessMatrix = this->GetStiffnessMatrix();
        SparseMatrixType& rQMatrix = this->GetQMatrix();

        // Initialize dummy rhs vector
        SparseVectorType b;
        SparseSpaceType::Resize(b,SparseSpaceType::Size1(rStiffnessMatrix));
        SparseSpaceType::Set(b,0.0);

        rModelPart.GetProcessInfo()[BUILD_LEVEL] = 1;
        TSparseSpace::SetToZero(rStiffnessMatrix);
        this->pGetBuilderAndSolver()->Build(pScheme,rModelPart,rStiffnessMatrix,b);

        if (BaseType::GetEchoLevel() == 4) {
            TSparseSpace::WriteMatrixMarketMatrix("StiffnessMatrix.mm", rStiffnessMatrix, false);
        }

        std::cout<<"here we go again"<<std::endl;
        
        rModelPart.GetProcessInfo()[BUILD_LEVEL] = 5;
        TSparseSpace::SetToZero(rQMatrix);
        this->pGetBuilderAndSolver()->Build(pScheme,rModelPart,rQMatrix,b);
        
        if (BaseType::GetEchoLevel() == 4) {
            TSparseSpace::WriteMatrixMarketMatrix("QMatrix.mm", rQMatrix, false);
        }

        // Find the DOFs on the current interface boundary
        Model new_model;               
        ModelPart& new_model_part = new_model.CreateModelPart("new_model"); 

        for (auto& r_cond : rModelPart.Conditions()) {
            
            auto& r_geom = r_cond.GetGeometry();
            auto& r_N = r_geom.ShapeFunctionsValues();

            for (IndexType i = 0; i<r_N.size2();++i)
            {
                if(r_N(0,i) > 1e-6)
                {
                    new_model_part.AddNode(r_geom.pGetPoint(i));
                }
            }
        }
        double number_of_nodes = new_model_part.NumberOfNodes();

        //create the result vector
        Vector rResult;
        rResult.resize(number_of_nodes*3);

        IndexType i = 0;
        for (auto& r_node : new_model_part.Nodes()) {

            const IndexType index = i * 3;
            
            rResult[index]     = r_node.GetDof(DISPLACEMENT_X).EquationId();
            rResult[index + 1] = r_node.GetDof(DISPLACEMENT_Y).EquationId();
            rResult[index + 2] = r_node.GetDof(DISPLACEMENT_Z).EquationId();

            i++;
        }

        SparseMatrixType reduced_Stiffness;
        reduced_Stiffness = ZeroMatrix(number_of_nodes*3,number_of_nodes*3);

        SparseMatrixType reduced_Q;
        reduced_Q = ZeroMatrix(number_of_nodes*3,number_of_nodes*3);

        for (IndexType i = 0; i < number_of_nodes*3; i++) 
        {
            for (IndexType j = 0; j <= i; j++) 
            {
                reduced_Stiffness(i,j) = rStiffnessMatrix(rResult(i),rResult(j));
                if (i != j)
                    reduced_Stiffness(j,i) = rStiffnessMatrix(rResult(i),rResult(j));
            }
        }

        for (IndexType i = 0; i < number_of_nodes*3; i++) 
        {
            for (IndexType j = 0; j <= i; j++) 
            {
                reduced_Q(i,j) = rQMatrix(rResult(i),rResult(j));
                if (i != j)
                    reduced_Q(j,i) = rQMatrix(rResult(i),rResult(j));
            }
        }

        // Eigenvector matrix and eigenvalue vector are initialized by the solver
        DenseVectorType Eigenvalues;
        DenseMatrixType Eigenvectors;

        // Solve for eigenvalues and eigenvectors
        BuiltinTimer system_solve_time;
        this->pGetBuilderAndSolver()->GetLinearSystemSolver()->Solve(
                reduced_Q,
                reduced_Stiffness,
                Eigenvalues,
                Eigenvectors);

        KRATOS_INFO_IF("System Solve Time", BaseType::GetEchoLevel() > 0)
                << system_solve_time.ElapsedSeconds() << std::endl;

        this->AssignVariables(Eigenvalues,Eigenvectors);

        return true;
        KRATOS_CATCH("")
    }

    void FinalizeSolutionStep() override
    {
        KRATOS_TRY;

        const int rank = BaseType::GetModelPart().GetCommunicator().MyPID();
        KRATOS_INFO_IF("EigensolverSupportNitscheStrategy", BaseType::GetEchoLevel() > 2 && rank == 0)
            <<  "Entering FinalizeSolutionStep" << std::endl;

        SparseMatrixType& rQMatrix = this->GetQMatrix();
        SparseVectorPointerType pDx = SparseSpaceType::CreateEmptyVectorPointer();
        SparseVectorPointerType pb = SparseSpaceType::CreateEmptyVectorPointer();
        pGetBuilderAndSolver()->FinalizeSolutionStep(
            BaseType::GetModelPart(), rQMatrix, *pDx, *pb);
        pGetScheme()->FinalizeSolutionStep(BaseType::GetModelPart(),
                                           rQMatrix, *pDx, *pb);
        KRATOS_INFO_IF("EigensolverSupportNitscheStrategy", BaseType::GetEchoLevel() > 2 && rank == 0)
            <<  "Exiting FinalizeSolutionStep" << std::endl;

        KRATOS_CATCH("");
    }

    /**
     * Function to perform expensive checks.
     * It is designed to be called ONCE to verify that the input is correct.
     */
    int Check() override
    {
        KRATOS_TRY

        ModelPart& rModelPart = BaseType::GetModelPart();
        const int rank = rModelPart.GetCommunicator().MyPID();

        KRATOS_INFO_IF("EigensolverSupportNitscheStrategy", BaseType::GetEchoLevel() > 2 && rank == 0)
            <<  "Entering Check" << std::endl;

        // check the model part
        BaseType::Check();

        // check the scheme
        this->pGetScheme()->Check(rModelPart);

        // check the builder and solver
        this->pGetBuilderAndSolver()->Check(rModelPart);

        KRATOS_INFO_IF("EigensolverSupportNitscheStrategy", BaseType::GetEchoLevel() > 2 && rank == 0)
            <<  "Exiting Check" << std::endl;

        return 0;

        KRATOS_CATCH("")
    }

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
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

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

    SchemePointerType mpScheme;

    BuilderAndSolverPointerType mpBuilderAndSolver;

    SparseMatrixPointerType mpStiffnessMatrix;

    SparseMatrixPointerType mpQMatrix;

    bool mInitializeWasPerformed = false;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{


    /// Assign eigenvalues and eigenvectors to kratos variables.
    void AssignVariables(DenseVectorType& rEigenvalues, DenseMatrixType& rEigenvectors)
    {
        ModelPart& rModelPart = BaseType::GetModelPart();
        const std::size_t NumEigenvalues = rEigenvalues.size();

        // store eigenvalues in process info
        rModelPart.GetProcessInfo()[EIGENVALUE_VECTOR] = rEigenvalues;

        const auto& r_dof_set = this->pGetBuilderAndSolver()->GetDofSet();

        for (ModelPart::NodeIterator itNode = rModelPart.NodesBegin(); itNode!= rModelPart.NodesEnd(); itNode++) {
            ModelPart::NodeType::DofsContainerType& NodeDofs = itNode->GetDofs();
            const std::size_t NumNodeDofs = NodeDofs.size();
            Matrix& rNodeEigenvectors = itNode->GetValue(EIGENVECTOR_MATRIX);
            if (rNodeEigenvectors.size1() != NumEigenvalues || rNodeEigenvectors.size2() != NumNodeDofs) {
                rNodeEigenvectors.resize(NumEigenvalues,NumNodeDofs,false);
            }

            // fill the EIGENVECTOR_MATRIX
            for (std::size_t i = 0; i < NumEigenvalues; i++) {
                for (std::size_t j = 0; j < NumNodeDofs; j++)
                {
                    const auto itDof = std::begin(NodeDofs) + j;
                    bool is_active = !(r_dof_set.find(**itDof) == r_dof_set.end());
                    if ((*itDof)->IsFree() && is_active) {
                       rNodeEigenvectors(i,j) = rEigenvectors(i,(*itDof)->EquationId());
                    }
                    else {
                       rNodeEigenvectors(i,j) = 0.0;
                    }
                }
            }
        }
    }


    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}

}; /* Class EigensolverSupportNitscheStrategy */

///@}

///@name Type Definitions
///@{


///@}

} /* namespace Kratos */

#endif /* KRATOS_EIGENSOLVER_SUPPORT_NITSCHE_STRATEGY  defined */
