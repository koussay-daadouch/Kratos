//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:   BSD License
//      Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "testing/testing.h"
#include "includes/model_part.h"
#include "utilities/auxiliar_model_part_utilities.h"

// Utilities
#include "utilities/cpp_tests_utilities.h"

namespace Kratos {
namespace Testing {

/**
* Checks the correct work of the Auxiliar model parts utility GetData for scalar data on Node_historicalDatalocation 
*/
KRATOS_TEST_CASE_IN_SUITE(AuxiliarModelPartUtilities_GetData_double_Node_historical, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& this_model_part = current_model.CreateModelPart("Main");
    this_model_part.AddNodalSolutionStepVariable(DISPLACEMENT_X);
    
    // Variable addition
    CppTestsUtilities::Create2DGeometry(this_model_part, "Element2D3N");

    //Assign random values to x_direction displacement to each node, also save it in test_values to compare later
    int i=0;
    double disp_test_value = 1.26;    
    std::vector<double> test_values;
    test_values.resize(this_model_part.NumberOfNodes());

    for (auto i_node = this_model_part.NodesBegin(); i_node != this_model_part.NodesEnd(); i_node++) 
    {
		i_node->FastGetSolutionStepValue(DISPLACEMENT_X) = static_cast<double>(disp_test_value * i);
        test_values[i] = (disp_test_value * i);
        i++;
    }

    std::cout << "Welcome to testing!" << std::endl;

    auto data = AuxiliarModelPartUtilities(this_model_part).GetVariableData(DISPLACEMENT_X, DataLocation::NodeHistorical);

    KRATOS_CHECK_VECTOR_NEAR(data, test_values, 0.004);

}

/**
* Checks the correct work of the Auxiliar model parts utility GetData for scalar data on NodeNonHistorical Datalocation
*/ 
KRATOS_TEST_CASE_IN_SUITE(AuxiliarModelPartUtilities_GetData_double_Node_Nonhistorical, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& this_model_part = current_model.CreateModelPart("Main");
    this_model_part.AddNodalSolutionStepVariable(DISPLACEMENT_X);
    
    // Variable addition
    CppTestsUtilities::Create2DGeometry(this_model_part, "Element2D3N");

    //Assign random values to x_direction displacement to each node, also save it in test_values to compare later
    int i=0;
    double disp_test_value = 1.26;    
    std::vector<double> test_values;
    test_values.resize(this_model_part.NumberOfNodes());

    for (auto i_node = this_model_part.NodesBegin(); i_node != this_model_part.NodesEnd(); i_node++) 
    {
		i_node->GetValue(DISPLACEMENT_X) = static_cast<double>(disp_test_value * i);
        test_values[i] = (disp_test_value * i);
        i++;
    }

    std::cout << "Welcome to testing!" << std::endl;

    auto data = AuxiliarModelPartUtilities(this_model_part).GetVariableData(DISPLACEMENT_X, DataLocation::NodeNonHistorical);

    KRATOS_CHECK_VECTOR_NEAR(data, test_values, 0.004);

}

/**
* Checks the correct work of the Auxiliar model parts utility GetData for scalar data on Element Datalocation
*/ 
KRATOS_TEST_CASE_IN_SUITE(AuxiliarModelPartUtilities_GetData_double_Element, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& this_model_part = current_model.CreateModelPart("Main");
    this_model_part.AddNodalSolutionStepVariable(DISPLACEMENT_X);
    
    // Variable addition
    CppTestsUtilities::Create2DGeometry(this_model_part, "Element2D3N", true, true);

    //Assign random values to x_direction displacement to each node, also save it in test_values to compare later
    int i=0;
    double disp_test_value = 2.55;    
    std::vector<double> test_values;
    test_values.resize(this_model_part.NumberOfNodes());

    for (auto i_elem = this_model_part.ElementsBegin(); i_elem != this_model_part.ElementsEnd(); i_elem++) 
    {
		i_elem->GetValue(DISPLACEMENT_X) = static_cast<double>(disp_test_value * i);
        test_values[i] = (disp_test_value * i);
        i++;
    }

    std::cout << "Welcome to testing!" << std::endl;

    auto data = AuxiliarModelPartUtilities(this_model_part).GetVariableData(DISPLACEMENT_X, DataLocation::Element);

    KRATOS_CHECK_VECTOR_NEAR(data, test_values, 0.004);

}


} // namespace Testing
} // namespace Kratos.
