from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics as KM
from . import recoverer

class StandardGradientRecoverer(recoverer.GradientRecoverer):
    def __init__(self, project_parameters, model_part):
        recoverer.GradientRecoverer.__init__(self, project_parameters, model_part)
    def RecoverGradientOfScalar(self, scalar_variable, gradient_variable):
        self.cplusplus_recovery_tool.CalculateGradient(self.model_part, scalar_variable, gradient_variable)
    def RecoverGradientOfVector(self, vector_variable, gradient_variable_x, gradient_variable_y, gradient_variable_z):
        self.cplusplus_recovery_tool.CalculateGradient(self.model_part, vector_variable, gradient_variable_x, gradient_variable_y, gradient_variable_z)
    def RecoverGradientOfVelocity(self):
        self.RecoverGradientOfVector(KM.VELOCITY, KM.VELOCITY_X_GRADIENT, KM.VELOCITY_Y_GRADIENT, KM.VELOCITY_Z_GRADIENT)
    def RecoverPressureGradient(self):
        self.RecoverGradientOfScalar(KM.PRESSURE, KM.RECOVERED_PRESSURE_GRADIENT)

class StandardMaterialAccelerationRecoverer(recoverer.MaterialAccelerationRecoverer):
    def __init__(self, project_parameters, model_part):
        recoverer.MaterialAccelerationRecoverer.__init__(self, project_parameters, model_part)
    def RecoverMaterialAcceleration(self):
        self.cplusplus_recovery_tool.CalculateVectorMaterialDerivative(self.model_part, KM.VELOCITY, KM.ACCELERATION, KM.MATERIAL_ACCELERATION)

class StandardLaplacianRecoverer(recoverer.LaplacianRecoverer):
    def __init__(self, project_parameters, model_part):
        recoverer.LaplacianRecoverer.__init__(self, project_parameters, model_part)
    def RecoverVectorLaplacian(self, vector_variable, laplacian_variable):
        self.cplusplus_recovery_tool.CalculateVectorLaplacian(self.model_part, vector_variable, laplacian_variable)
