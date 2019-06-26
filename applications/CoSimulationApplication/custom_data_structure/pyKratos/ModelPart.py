from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# pyKratos imports
from .Node import Node
from .Variables import *

# Other imports
from collections import OrderedDict

class ModelPart(object):

    class PointerVectorSet(OrderedDict):
        def __iter__(self):
            self.vals_list = iter(list(self.values()))
            return self

        def __next__(self):
            return next(self.vals_list)

    def __init__(self, name="default", buffer_size=1, internal_construction=False):
        if not internal_construction:
            raise Exception("Creation of standalone ModelParts is not possible, please use Model.CreateModelPart()!")

        self.__parent_model_part = None
        self.__sub_model_parts = PointerVectorSet()
        self.__nodes           = PointerVectorSet()
        self.__elements        = PointerVectorSet()
        self.__properties      = PointerVectorSet()

        self.__buffer_size = buffer_size

        self.__hist_variables = set()

        if("." in name):
            RuntimeError("Name of the modelpart cannot contain a . (dot) Please rename ! ")
        if(name == ""):
            RuntimeError("No empty names for modelpart are allowed. Please rename ! ")

        self.Name = name
        self.ProcessInfo = {TIME: 0.0, DELTA_TIME: 0.0}  # empty dictionary

    ### Methods related to historical variables ###

    def AddNodalSolutionStepVariable(self, variable):
        if self.NumberOfNodes() > 0:
            # this is forbidden since it creates problems with the memory management of historical variables
            raise Exception("Variables can only be added before adding Nodes!")

        self.__hist_variables.add(variable)

    def CloneSolutionStep(self):
        for node in self.Nodes:
            node.CloneSolutionStep()

    def CloneTimeStep(self, time):
        self.CloneSolutionStep()

        # the following is the equvalent of "SetAsTimeStepInfo" on the ProcessInfo
        old_time = self.ProcessInfo[TIME]
        self.ProcessInfo[TIME] = time
        self.ProcessInfo[DELTA_TIME] = time-old_time

    ### Methods related to SubModelParts ###
    @property
    def SubModelParts(self):
        return self.__sub_model_parts

    def CreateSubModelPart(self, name_smp):
        if name_smp in self.__sub_model_parts:
            raise RuntimeError('There is an already existing sub model part with name "{}" in model part: "{}"'.format(name_smp, self.Name))
        smp = ModelPart(name_smp, self.__buffer_size, True)
        smp.__parent_model_part = self
        smp.ProcessInfo = self.ProcessInfo

        self.__sub_model_parts[name_smp] = smp
        return smp

    def HasSubModelPart(self, name_smp):
        return name_smp in self.__sub_model_parts

    def GetSubModelPart(self, smp_name):
        try:
            return self.__sub_model_parts[smp_name]
        except KeyError:
            raise RuntimeError('SubModelPart "{}" not found'.format(smp_name))

    def IsSubModelPart(self):
        return self.__parent_model_part != None

    ### Methods related to Nodes ###
    @property
    def Nodes(self):
        return self.__nodes

    def NumberOfNodes(self):
        return len(self.__nodes)

    def GetNode(self, node_id):
        try:
            return self.__nodes[node_id]
        except KeyError:
            raise RuntimeError('Node index not found: {}'.format(node_id))

    def CreateNewNode(self, node_id, coord_x, coord_y, coord_z):
        if self.IsSubModelPart():
            new_node = self.__parent_model_part.CreateNewNode(node_id, coord_x, coord_y, coord_z, self.__hist_variables, self.__buffer_size)
            self.__nodes[node_id] = new_node
            return new_node
        else:
            existing_node = self.__nodes.get(node_id)
            if existing_node:
                if self.__Distance(existing_node.Coordinates(), [coord_x, coord_y, coord_z]) > 1E-15:
                    err_msg  = 'A node with Id #' + str(node_id) + ' already exists in the root model part with different Coordinates!'
                    err_msg += '\nExisting Coords: ' + str(existing_node.Coordinates())
                    err_msg += '\nNew Coords: '      + str([coord_x, coord_y, coord_z])
                    raise RuntimeError(err_msg)

                return existing_node
            else:
                new_node = self.Node(node_id, coord_x, coord_y, coord_z)
                self.__nodes[node_id] = new_node
                return new_node


    ### Methods related to Elements ###
    @property
    def Elements(self):
        return self.__elements

    def NumberOfElements(self):
        return len(self.__elements)

    def GetElement(self, element_id):
        try:
            return self.__elements[element_id]
        except KeyError:
            raise RuntimeError('Element index not found: {}'.format(element_id))

    def CreateNewElement(self, element_name, element_id, node_ids, property_id):
        if self.IsSubModelPart():
            new_element = self.__parent_model_part.CreateNewElement(element_name, element_id, node_ids, property_id)
            self.__elements[element_id] = new_element
            self.AddElement(new_element)
            return new_element
        else:
            element_nodes = [self.GetNode(node_id) for node_id in node_ids]
            new_element = Element(element_name, element_id, element_nodes, property_id) # TODO pass property? Or at least check if this property exists ...
            if element_id in self.__elements:
                existing_element = self.__elements[element_id]
                if existing_element != new_element:
                    raise RuntimeError('A different element with the same Id exists already!') # TODO check what Kratos does here

                return existing_element
            else:
                self.__elements[element_id] = new_element
                return new_element

    # add properties
    def AddProperties(self, dict_of_properties):
        self.Properties.update(dict_of_properties)

    def CreateNewProperties(self, prop_id):
        return None

    def RemoveElement(self):
        pass

    def HasNodalSolutionStepVariable(self):
        pass

    def WriteMesh(self):
        pass

    def Check(self):
        pass

    def __str__(self):
        return "ModelPart:\n    Number of Nodes: {0}\n    Nunber of Elements: {1}".format(len(self.NodesMap), len(self.ElementsMap))

