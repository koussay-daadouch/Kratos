//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: Nelson Lafontaine $
//   Date:                $Date: 2008-10-13 08:56:42 $
//   Revision:            $Revision: 1.5 $
//
//

#if !defined(KRATOS_MPI_EXPLICIT_SOLVER_STRATEGY)
#define  KRATOS_MPI_EXPLICIT_SOLVER_STRATEGY

// /* External includes */
// #include "boost/smart_ptr.hpp"

// System includes

// Project includes
#include "utilities/timer.h"
#include "custom_utilities/neighbours_calculator.h"
#include "custom_utilities/mpi_neighbours_calculator.h"
#include "custom_utilities/create_and_destroy.h"

#include "custom_elements/spheric_particle.h" //M: le afegit jo.. no hi era. cal que hi sigui oi???
#include "includes/variables.h"

/* System includes */
#include <limits>
#include <iostream>
#include <iomanip>
#include <iostream>

/* External includes */
#ifdef _OPENMP
#include <omp.h>
#endif

#include "boost/smart_ptr.hpp"

/* Project includes */
#include "includes/define.h"
#include "utilities/openmp_utils.h"
#include "includes/model_part.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/schemes/scheme.h"

// #include "custom_utilities/neighbours_calculator.h"
#include "custom_utilities/mpi_neighbours_calculator.h"
#include "custom_strategies/schemes/integration_scheme.h"
#include "custom_strategies/strategies/explicit_solver_strategy.h"

namespace Kratos
{
  ///@addtogroup ApplicationNameApplication
  ///@{

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
    
  bool gcontanct(boost::shared_ptr<Element> it1, boost::shared_ptr<Element> it2) 
  { 
      int idp11 = it1->GetGeometry()(0)->Id();
      int idp12 = it1->GetGeometry()(1)->Id();
      int idp21 = it2->GetGeometry()(0)->Id();
      int idp22 = it2->GetGeometry()(1)->Id();
      
      int it1min = idp11 < idp12 ? idp11 : idp12;
      int it2min = idp21 < idp22 ? idp21 : idp22;
      
      if      ( idp11 <  idp21 ) return true;
      else if ( idp11 == idp21 )
      {
          if(idp12 <  idp22 ) return true;
          else return false;
      }
      else return false;
  }
  
  /// Short class definition.
  /** Detail class definition.
  */
  template<
  class TSparseSpace,
  class TDenseSpace, 
  class TLinearSolver> 
  class MpiExplicitSolverStrategy : public ExplicitSolverStrategy<TSparseSpace,TDenseSpace,TLinearSolver>
     {
      public:
      ///@name Type Definitions
      ///@{
 
      typedef SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver>   BaseType;
      typedef typename BaseType::TDataType                              TDataType;
      typedef typename BaseType::TBuilderAndSolverType                  TBuilderAndSolverType;
      typedef typename BaseType::TSchemeType                            TSchemeType;
      typedef typename BaseType::DofsArrayType                          DofsArrayType;
      typedef typename Element::DofsVectorType                          DofsVectorType;
      
      typedef ModelPart::NodesContainerType                             NodesArrayType;
      typedef ModelPart::ElementsContainerType                          ElementsArrayType;
      typedef ModelPart::ConditionsContainerType                        ConditionsArrayType;     
      
      typedef ModelPart::NodesContainerType::ContainerType              NodesContainerType;
      typedef ModelPart::ElementsContainerType::ContainerType           ElementsContainerType;
      typedef ModelPart::ConditionsContainerType::ContainerType         ConditionsContainerType;
      
      typedef DiscreteElement                                           ParticleType;
      
      typedef WeakPointerVector<Element> ParticleWeakVectorType; 
      typedef WeakPointerVector<Element >::iterator ParticleWeakIteratorType;
      typedef ParticleWeakVectorType::ptr_iterator ParticleWeakIteratorType_ptr;
      
      /// Pointer definition of ExplicitSolverStrategy
      KRATOS_CLASS_POINTER_DEFINITION(MpiExplicitSolverStrategy);
 
      
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      MpiExplicitSolverStrategy(){}
      
      MpiExplicitSolverStrategy(ModelPart& model_part,
                             ModelPart& contacts_model_part, 
                             const int dimension,
                             const double enlargement_factor,
                             const double damping_ratio,
                             const double fraction_delta_time,
                             const double max_delta_time,
                             const double n_step_search,
                             const double safety_factor,
                             const bool MoveMeshFlag,
                             const bool delta_option,
                             const bool continuum_simulating_option,
                             typename IntegrationScheme::Pointer pScheme
      ) : ExplicitSolverStrategy<TSparseSpace,TDenseSpace,TLinearSolver>(
                             model_part,
                             contacts_model_part,
                             dimension,
                             enlargement_factor,
                             damping_ratio,
                             fraction_delta_time,
                             max_delta_time,
                             n_step_search,
                             safety_factor,
                             MoveMeshFlag,
                             delta_option,
                             continuum_simulating_option,
                             pScheme )
      {
      }

      /// Destructor.
      virtual ~MpiExplicitSolverStrategy(){}
          
    protected:

    private:
 
    typename IntegrationScheme::Pointer mpScheme;
    
    virtual void Synchronize(ModelPart& r_model_part, ModelPart& r_contact_model_part)
    {   
        r_model_part.GetCommunicator().SynchronizeNodalSolutionStepsData();
        r_model_part.GetCommunicator().SynchronizeDofs();
        
        //ContactModelpart Element Synchronize       
        r_contact_model_part.GetCommunicator().SynchronizeElementalNonHistoricalVariable(CONTACT_FAILURE);    
        r_contact_model_part.GetCommunicator().SynchronizeElementalNonHistoricalVariable(CONTACT_SIGMA);
        r_contact_model_part.GetCommunicator().SynchronizeElementalNonHistoricalVariable(CONTACT_TAU);
        r_contact_model_part.GetCommunicator().SynchronizeElementalNonHistoricalVariable(FAILURE_CRITERION_STATE);
        r_contact_model_part.GetCommunicator().SynchronizeElementalNonHistoricalVariable(FAILURE_CRITERION_OPTION);
 
        r_contact_model_part.GetCommunicator().SynchronizeElementalNonHistoricalVariable(CONTACT_SIGMA_MAX);
        r_contact_model_part.GetCommunicator().SynchronizeElementalNonHistoricalVariable(CONTACT_SIGMA_MIN);
        r_contact_model_part.GetCommunicator().SynchronizeElementalNonHistoricalVariable(CONTACT_TAU_ZERO);
        r_contact_model_part.GetCommunicator().SynchronizeElementalNonHistoricalVariable(CONTACT_INTERNAL_FRICC);

        r_contact_model_part.GetCommunicator().SynchronizeElementalNonHistoricalVariable(LOCAL_CONTACT_FORCE_LOW);
        r_contact_model_part.GetCommunicator().SynchronizeElementalNonHistoricalVariable(LOCAL_CONTACT_FORCE_HIGH);
    }
    
    virtual void Repart(ModelPart& r_model_part)
    {
        typedef Mpi_Neighbours_Calculator<ParticleType> NeighboursCalculatorType;
        
        NeighboursCalculatorType::Parallel_partitioning(r_model_part,true);
    }
    
    virtual ElementsArrayType& GetElements(ModelPart& r_model_part)
    {
        return r_model_part.GetCommunicator().LocalMesh().Elements();
    }

    virtual void SearchIniNeighbours(ModelPart& r_model_part,bool extension_option)
    { 
        typedef Mpi_Neighbours_Calculator<ParticleType> NeighboursCalculatorType;
              
        NeighboursCalculatorType neighbourCalc;
        neighbourCalc.Search_Ini_Neighbours(r_model_part, extension_option);
    }//SearchIniNeighbours


    virtual void SearchNeighbours(ModelPart& r_model_part,bool extension_option)
    {
        typedef Mpi_Neighbours_Calculator<ParticleType> NeighboursCalculatorType;
              
        NeighboursCalculatorType neighbourCalc;
        neighbourCalc.Search_Neighbours(r_model_part, extension_option);
    }//SearchNeighbours
    
    virtual void PrepareContactModelPart(ModelPart& r_model_part, ModelPart& mcontacts_model_part)
    {  
        mcontacts_model_part.GetCommunicator().SetNumberOfColors(r_model_part.GetCommunicator().GetNumberOfColors());
        mcontacts_model_part.GetCommunicator().NeighbourIndices() = r_model_part.GetCommunicator().NeighbourIndices();
    }
    
    virtual bool ContactElementsParallelCondition(ElementsArrayType::ptr_iterator it, ParticleWeakIteratorType_ptr continuum_ini_neighbour_iterator)
    {
        return ((*it)->GetValue(PARTITION_INDEX) != (*continuum_ini_neighbour_iterator).lock()->GetValue(PARTITION_INDEX));
    }
        
    //En aquest cas m'afegeixo jo a mi mateix
    virtual void Add_As_Own(ModelPart& r_model_part, ModelPart& mcontacts_model_part, ParticleWeakIteratorType_ptr continuum_ini_neighbour_iterator, Element::Pointer p_contact_element)
    {
        KRATOS_TRY
        
        mcontacts_model_part.Elements().push_back(p_contact_element);
        mcontacts_model_part.GetCommunicator().LocalMesh().Elements().push_back(p_contact_element);
        
        KRATOS_CATCH("")
    }
    
    //En aquest cas m'afegeixo jo al local i a la interface local corresponent amb la particio del vei ghost
    virtual void Add_As_Local(ModelPart& r_model_part, ModelPart& mcontacts_model_part, ParticleWeakIteratorType_ptr continuum_ini_neighbour_iterator, Element::Pointer p_contact_element)
    {
        KRATOS_TRY
        
        mcontacts_model_part.Elements().push_back(p_contact_element);
        
        Communicator::NeighbourIndicesContainerType communicator_ranks = r_model_part.GetCommunicator().NeighbourIndices();
        
        int NumberOfRanks = r_model_part.GetCommunicator().GetNumberOfColors();
        int destination = -1;
      
        for(int i = 0; i < NumberOfRanks; i++)
            if((*continuum_ini_neighbour_iterator).lock()->GetGeometry()(0)->GetSolutionStepValue(PARTITION_INDEX) == communicator_ranks[i])
                destination = i;
            
        mcontacts_model_part.GetCommunicator().LocalMesh().Elements().push_back(p_contact_element);
                        
        if(destination > -1)
        {   
            mcontacts_model_part.GetCommunicator().LocalMesh(destination).Elements().push_back(p_contact_element);
        }
        
        KRATOS_CATCH("")
    }
    
    //I aqui m'afegeixio yo com a ghost de la particio del vei local
    virtual void Add_As_Ghost(ModelPart& r_model_part, ModelPart& mcontacts_model_part, ParticleWeakIteratorType_ptr continuum_ini_neighbour_iterator, Element::Pointer p_contact_element)
    {
        KRATOS_TRY
        
//         mcontacts_model_part.Elements().push_back(p_contact_element);
        
        Communicator::NeighbourIndicesContainerType communicator_ranks = r_model_part.GetCommunicator().NeighbourIndices();
        
        int NumberOfRanks = r_model_part.GetCommunicator().GetNumberOfColors();
        int destination = -1;
      
        for(int i = 0; i < NumberOfRanks; i++)
            if((*continuum_ini_neighbour_iterator).lock()->GetGeometry()(0)->GetSolutionStepValue(PARTITION_INDEX) == communicator_ranks[i])
                destination = i;
                        
        if(destination > -1)
        {   
            mcontacts_model_part.GetCommunicator().GhostMesh().Elements().push_back(p_contact_element);
            mcontacts_model_part.GetCommunicator().GhostMesh(destination).Elements().push_back(p_contact_element);
        }
        
        KRATOS_CATCH("")
    }
    
    virtual void Sort_Contact_Modelpart(ModelPart& mcontacts_model_part)
    {   
//         std::sort(mcontacts_model_part.Elements().begin(),mcontacts_model_part.Elements().end(),gcontanct);
//         
//         std::sort(mcontacts_model_part.GetCommunicator().LocalMesh().Elements().ptr_begin(),mcontacts_model_part.GetCommunicator().LocalMesh().Elements().ptr_end(),gcontanct);
//         std::sort(mcontacts_model_part.GetCommunicator().GhostMesh().Elements().ptr_begin(),mcontacts_model_part.GetCommunicator().GhostMesh().Elements().ptr_end(),gcontanct);
      
        for(int i = 0; i < mcontacts_model_part.GetCommunicator().GetNumberOfColors(); i++)
        {
            std::sort(mcontacts_model_part.GetCommunicator().LocalMesh(i).Elements().ptr_begin(),mcontacts_model_part.GetCommunicator().LocalMesh(i).Elements().ptr_end(),gcontanct);
            std::sort(mcontacts_model_part.GetCommunicator().GhostMesh(i).Elements().ptr_begin(),mcontacts_model_part.GetCommunicator().GhostMesh(i).Elements().ptr_end(),gcontanct);
        }
    }
    
    virtual void Reassign_Ids(ModelPart& mcontacts_model_part)
    {
        int contacts_model_part_size = mcontacts_model_part.GetCommunicator().LocalMesh().Elements().size();
        int iteratorId = -1;
       
        MpiDiscreteParticleConfigure<3>::ReduceIds(contacts_model_part_size,iteratorId);
        
        if(iteratorId == -1)
            std::cout << "Something went wrong :(" << std::endl;
        
        for (ElementsArrayType::ptr_iterator it = mcontacts_model_part.GetCommunicator().LocalMesh().Elements().ptr_begin(); it != mcontacts_model_part.GetCommunicator().LocalMesh().Elements().ptr_end(); ++it)
        {
            (*it)->SetId(iteratorId++);
        }
        
        mcontacts_model_part.GetCommunicator().SynchronizeElementalIds();
    }

  
  }; // Class MpiExplicitSolverStrategy  


        
 /*
  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream, 
                    MpiExplicitSolverStrategy& rThis){return rIStream;}

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream, 
                    const MpiExplicitSolverStrategy& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
    */
  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_FILENAME_H_INCLUDED  defined 




