# We import the libraries
import KratosMultiphysics
import KratosMultiphysics.mpi as KratosMPI
import KratosMultiphysics.MeshingApplication as MeshingApplication
from KratosMultiphysics.compare_two_files_check_process import CompareTwoFilesCheckProcess


import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.testing.utilities import ReadModelPart
import json
import os
import math

def GetFilePath(fileName):
    return os.path.dirname(os.path.realpath(__file__)) + "/" + fileName

class TestMPIParMmg(KratosUnittest.TestCase):

    def test_mpi_sphere(self):
        KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)

        # We create the model part
        current_model = KratosMultiphysics.Model()
        main_model_part = current_model.CreateModelPart("MainModelPart")
        main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 3)

        # We add the variables needed
        main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
        main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE_GRADIENT)


        # We import the model main_model_part
        file_path = GetFilePath("/parmmg_eulerian_test/background_mesh_sphere")
        ReadModelPart(file_path, main_model_part)

        communicator = main_model_part.GetCommunicator().GetDataCommunicator()

        for node in main_model_part.Nodes:
            distance = math.sqrt(node.X**2+node.Y**2+node.Z**2) - 1.0/2.0
            node.SetSolutionStepValue(KratosMultiphysics.DISTANCE,distance)

        ##COMPUTE DISTANCE GRADIENT AND NODAL_H
        local_gradient = KratosMultiphysics.ComputeNodalGradientProcess3D(main_model_part,
        KratosMultiphysics.DISTANCE,
        KratosMultiphysics.DISTANCE_GRADIENT,
        KratosMultiphysics.NODAL_AREA)

        local_gradient.Execute()
        find_nodal_h = KratosMultiphysics.FindNodalHNonHistoricalProcess(main_model_part)
        find_nodal_h.Execute()

        ##COMPUTE LEVEL SET METRIC
        metric_parameters = KratosMultiphysics.Parameters("""
        {
            "minimal_size"                             : 1.5,
            "sizing_parameters": {
                "reference_variable_name"               : "DISTANCE",
                "boundary_layer_max_distance"           : 2.0,
                "interpolation"                         : "constant"
            },
            "enforce_current"                      : false,
            "anisotropy_remeshing"                 : false
        }
        """)
        metric_process = KratosMultiphysics.MeshingApplication.ComputeLevelSetSolMetricProcess3D(
            main_model_part,
            KratosMultiphysics.DISTANCE_GRADIENT,
            metric_parameters)
        metric_process.Execute()

        ##PERFORM REMESHING
        pmmg_parameters = KratosMultiphysics.Parameters("""
        {
            "save_external_files"              : true,
            "save_colors_files"                : true,
            "initialize_entities"              : false,
            "preserve_flags"                   : false
        }
        """)
        pmmg_process = KratosMultiphysics.MeshingApplication.ParMmgProcess3D(main_model_part.GetRootModelPart(), pmmg_parameters)
        pmmg_process.Execute()

        KratosMultiphysics.ModelPartIO("remeshed_sphere_"+str(communicator.Rank()), KratosMultiphysics.IO.WRITE ).WriteModelPart(main_model_part)

        # We check the solution
        check_parameters = KratosMultiphysics.Parameters("""
        {
            "reference_file_name"   : "reference_file_name",
            "output_file_name"      : "output_file_name",
            "comparison_type"       : "deterministic"

        }
        """)
        check_parameters["reference_file_name"].SetString("parmmg_eulerian_test/parmmg_sphere_reference_mdpa_rank_"+str(communicator.Rank())+".mdpa")
        check_parameters["output_file_name"].SetString("remeshed_sphere_"+str(communicator.Rank())+".mdpa")
        check_files = CompareTwoFilesCheckProcess(check_parameters)

        check_files.ExecuteInitialize()
        check_files.ExecuteBeforeSolutionLoop()
        check_files.ExecuteInitializeSolutionStep()
        check_files.ExecuteFinalizeSolutionStep()
        check_files.ExecuteFinalize()



if __name__ == '__main__':
    KratosUnittest.main()