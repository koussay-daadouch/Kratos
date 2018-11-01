# import Kratos
import KratosMultiphysics
import KratosMultiphysics.ParticleMechanicsApplication
import run_cpp_unit_tests

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import from Test Factories (with general analysis flows)
from particle_mechanics_test_factory import CLLinearElastic2DQuadTest as TCLLinearElastic2DQuadTest

# Import from Test Factories (with differet analysis flows)
# from test_patch_test_particles import TestLoadingConditionsPoint

def AssembleTestSuites():
    ''' Populates the test suites to run.

    Populates the test suites to run. At least, it should populate the suites:
    "small", "nighlty" and "all"

    Return
    ------

    suites: A dictionary of suites
        The set of suites with its test_cases added.
    '''

    suites = KratosUnittest.KratosSuites

    ### Create a test suit with the selected tests (Small tests):
    ### These tests have to be very fast!
    ### Execution time << 1 sec on a regular PC !!!

    ## These tests are executed by the continuous integration tool
    smallSuite = suites['small']
    # smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestLoadingConditionsPoint]))

    ## These tests are executed in the nightly build
    nightSuite = suites['nightly']
    nightSuite.addTests(smallSuite)
    nightSuite.addTest(TCLLinearElastic2DQuadTest('test_execution'))

    ### Adding Validation Tests
    ## For very long tests that should not be in nighly and you can use to validate
    validationSuite = suites['validation']

    ### Create a test suit that contains all the tests:
    allSuite = suites['all']
    allSuite.addTests(nightSuite) # already contains the smallSuite
    allSuite.addTests(validationSuite)

    return suites

if __name__ == '__main__':
    KratosMultiphysics.Logger.PrintInfo("Unittests", "\nRunning cpp unit tests ...")
    run_cpp_unit_tests.run()
    KratosMultiphysics.Logger.PrintInfo("Unittests", "Finished running cpp unit tests!")

    KratosMultiphysics.Logger.PrintInfo("Unittests", "\nRunning python tests ...")
    KratosUnittest.runTests(AssembleTestSuites())
    KratosMultiphysics.Logger.PrintInfo("Unittests", "Finished python tests!")
