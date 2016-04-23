# ##### BEGIN GPL LICENSE BLOCK #####
#
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation; either version 2
#  of the License, or (at your option) any later version.#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software Foundation,
#  Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
#
# ##### END GPL LICENSE BLOCK #####

# <pep8 compliant>

"""
This script stores the operators for HexBlender, which is what happens
after a button is pressed on the GUI
"""

# blender imports
import bpy
import bmesh
from bpy.app.handlers import persistent
from bl_operators.presets import AddPresetBase
from bpy_extras.io_utils import ExportHelper
from mathutils import Vector

# hexblender imports
from hexblender import hexblender_helpers 
from hexblender.hex_interp_subdiv import hex_interp_subdiv
from hexblender.hex_interp_subdiv_no_priorities import hex_interp_subdiv_no_priorities
from hexblender.regularize_elements import regularize_elements
from hexblender import find_hex
from hexblender.regularize_surface import regularize_surface

# python imports
import mathutils
import array
import glob
import os
import random
import re
import subprocess
import time
import shutil
import pickle
import json
import numpy as np
from math import pi
from scipy.linalg import logm

# Use per module class registration/unregistration
def register():
    bpy.utils.register_module(__name__)


def unregister():
    bpy.utils.unregister_module(__name__)


#HexBlender Operators:
class HEXBLENDER_OT_import_pickle(bpy.types.Operator):
    """Export vertices in Continuity format"""

    bl_idname = "hexblender.import_pickle"
    bl_label = "Import Pickle File"
    bl_description = ("Import a pickle file, generally created in Continuity.")
    bl_options = {'REGISTER'}

    filepath = bpy.props.StringProperty(subtype = "FILE_PATH")

    @classmethod
    def poll(cls, context):
        rd = context.scene.render
        return context.scene

    def execute(self, context):
        now = time.time()
        fid = open(self.filepath,'rb')
        try:
            data = pickle.load(fid)
        except Exception as msg:
            print("Problem importing pickle: %s" % msg)
            print("  Will try to import using latin1 encoding")
            data = pickle.load(fid, encoding='latin1')
        print("reading pickle file took: %f" % (time.time() - now))

        verts = data["verts"]
        
        try:
            faces = data["faces"]
            facesKeyPres = True
        except:
            facesKeyPres = False
        
        mats = data["mats"]
        elems = data["elems"]
        
        # two separate calls because if faces are present, the numbering of 
        # faces will be preserved. If they are not present, the  numbering of 
        # faces will NOT be preserved, but the vertex numbers will be, of course.
        if facesKeyPres:
            # This will return an empy new_cube_faces
            new_cube_faces = hexblender_helpers.add_new_data(np.array(verts), 
                                                             [], 
                                                             faces, 
                                                             [], 
                                                             newObj=True)

            # here, not using continuity ordering convention so change ordering
            # takes the clockwise order of verts, all 6 faces
            new_cube_faces = []                               
            for cube in data['elems_tr']: 
                new_cube_faces.append([int(value) for value in cube[:4]])
                new_cube_faces.append([int(value) for value in cube[4:]])
                for i in range(4):
                    new_cube_faces.append([int(cube[i]),
                                           int(cube[(i+1)%4]),
                                           int(cube[4+(i+1)%4]),
                                           int(cube[(i+4)])])
        else:
            # only read in elems and make elems_tr de novo so people importing 
            # from continuity don't have to deal with it
            elems_tr = []
            for hex in elems:
                elems_tr.append([hex[0],hex[1],hex[3],hex[2],hex[4],hex[5],hex[7],hex[6]])
            new_cube_faces = hexblender_helpers.add_new_data(np.array(verts), 
                                                             [], 
                                                             [], 
                                                             elems_tr, 
                                                             newObj=True)
        print("Called add_new_data: %f" % (time.time() - now))

        print("Created new_cube_faces: %f" % (time.time() - now))
        # add materials to faces
        hexblender_helpers.set_mat_data_hex(mats, 
                                            new_cube_faces,
                                            subdivided = False)
        print("Called set_mat_data: %f" % (time.time() - now))

        # cache the original element/material data from the Continuity pickle 
        # file which makes our lives easier when exporting
        hexblender_helpers.cache_data(elems, mats)
        print("Called cached_data: %f" % (time.time() - now))

        print("Successfully imported %s" % self.filepath)
        return {'FINISHED'}

    def invoke(self, context, event):
        context.window_manager.fileselect_add(self)
        return {'RUNNING_MODAL'}


class HEXBLENDER_OT_export_vertices(bpy.types.Operator, ExportHelper):
    """Export vertices in Continuity format"""

    bl_idname = "hexblender.export_vertices"
    bl_label = "Export Vertices"
    bl_description = ("Export vertices in Continuity format.")
    bl_options = {'REGISTER'}

    # ExportHelper mixin class uses this
    filename_ext = ".txt"

    filter_glob = bpy.props.StringProperty(
            default="*.txt",
            options={'HIDDEN'},
            )

    @classmethod
    def poll(cls, context):
        return context.object is not None

    def execute(self, context):
        return_code, vert_text = hexblender_helpers.write_out_vertices(self.filepath)

        return return_code

    def invoke(self, context, event):
        self.filepath = "verts.txt"
        context.window_manager.fileselect_add(self)
        return {'RUNNING_MODAL'}

class HEXBLENDER_OT_export_elements(bpy.types.Operator, ExportHelper):
    """Export elements in Continuity format"""

    bl_idname = "hexblender.export_elements"
    bl_label = "Export Elements"
    bl_description = ("Export elements in Continuity format.")
    bl_options = {'REGISTER'}

    # ExportHelper mixin class uses this
    filename_ext = ".txt"

    filter_glob = bpy.props.StringProperty(
            default="*.txt",
            options={'HIDDEN'},
            )

    @classmethod
    def poll(cls, context):
        return context.object is not None

    def execute(self, context):
        return_code, elem_string = hexblender_helpers.write_out_3d_elements(self.filepath)

        return return_code

    def invoke(self, context, event):
        self.filepath = "hex_elems.txt"
        context.window_manager.fileselect_add(self)
        return {'RUNNING_MODAL'}


########### MD 1/12/2016 #################################################
class HEXBLENDER_OT_export_vertexweights(bpy.types.Operator, ExportHelper):
    # created by Melody Dong 1/2016
    """Export scaled vertex weights to .txt file"""

    bl_idname = "hexblender.export_vertex_weights"
    bl_label = "Export Vertex Weights"
    bl_description = ("Export vertex weights to .txt file.")
    bl_options = {'REGISTER'}

    # ExportHelper mixin class uses this
    filename_ext = ".txt"

    filter_glob = bpy.props.StringProperty(
            default="*.txt",
            options={'HIDDEN'},
            )

    @classmethod
    def poll(cls, context):
        return context.object is not None

    def execute(self, context):
        return_code, vertweight_text = hexblender_helpers.compute_write_out_vertex_weights(context, self.filepath)
        return return_code

    def invoke(self, context, event):
        self.filepath = "vert_weights.txt"
        context.window_manager.fileselect_add(self)
        return {'RUNNING_MODAL'}     
############################################################################

class HEXBLENDER_OT_export_hermite_tricubic_derivs(bpy.types.Operator, ExportHelper):
    """Export Hermite Tricubic derivatives in Continuity format"""

    bl_idname = "hexblender.export_hermite_tricubic_derivs"
    bl_label = "Export Hermite Tricubic Derivatives"
    bl_description = ("Export Hermite Tricubic Derivatives in Continuity format.")
    bl_options = {'REGISTER'}

    # ExportHelper mixin class uses this
    filename_ext = ".txt"

    filter_glob = bpy.props.StringProperty(
            default="*.txt",
            options={'HIDDEN'},
            )

    @classmethod
    def poll(cls, context):
        return context.object is not None

    def execute(self, context):
        return_code, deriv_text = hexblender_helpers.compute_write_hermite_3d_derivs(context, self.filepath)

        return return_code

    def invoke(self, context, event):
        self.filepath = "tricubic_derivs.txt"
        context.window_manager.fileselect_add(self)
        return {'RUNNING_MODAL'}


class HEXBLENDER_OT_regularize_hexs(bpy.types.Operator):
    """Regularize hexahedrons"""

    bl_idname = "hexblender.regularize_hexs"
    bl_label = "Regularize Hexs"
    bl_description = ("Regularize hexahedrons using geometric flow method (Zhang et al, 2007, Surface Smoothing and Quality Improvement of Quadrilateral/Hexahedral Meshes with Geometric Flow)")
    bl_options = {'REGISTER'}

    @classmethod
    def poll(cls, context):
        return context.object is not None

    def execute(self, context):
        mesh, orig_mode = hexblender_helpers.get_active_mesh()

        b_mesh = bmesh.new()
        b_mesh.from_mesh(mesh)

        # verts = b_mesh.verts
        # faces = b_mesh.faces

        b_verts = b_mesh.verts
        b_faces = b_mesh.faces

        try:
            hex_scene = context.scene.hexblender
        except Exception as msg:
            print("ERROR: Unable to get hexblender scene!  %s" % msg)
            return {'FINISHED'}

        # input parameters
        preserveRidges=hex_scene.hexblender_properties.preserve_ridges
        immobExtraord=hex_scene.hexblender_properties.immobilize_extraord
        tangentMotion=hex_scene.hexblender_properties.tangent_motion
        normalMotion=hex_scene.hexblender_properties.normal_motion
        internalMotion=hex_scene.hexblender_properties.internal_motion

        # itt is the number of iterations set in the HexBlender panel GUI
        itt = hex_scene.hexblender_properties.regularize_hex_iters

        print("\n Hex Regularization:\t\t by Greg Sturgeon")
        print(" -------------------")
     
        # See if anything is selected.
        # XXX: I don't like this implementation, but I haven't come across
        # an easier/faster way to see if any verts are selected in 
        # Blender's API
        b_verts.ensure_lookup_table()
        verts = [v for v in b_verts if v.select]

        if len(verts) == 0:
            print("Error: Nothing is selected! Please enter Edit Mode and select vertices of hexes to regularize.")
            return {'FINISHED'}
        elif len(verts)==len(b_verts):
            print("Everything selected. Getting cached hex data...")
            cubes, hex_mat = hexblender_helpers.get_cached_data(mesh)
        else:
            # We want to only have an execution of the hex regularization script
            # if there are hexes selected. And we only want to regularize only these selected hexes.
            # If the entire mesh is selected (as with pressing "A") the output of FindHex(selectionOnly=True)
            # is the same as with selectionOnly=False
            cubes, hex_mat = hexblender_helpers.find_hex(selectionOnly=True, verts=verts)
            # CV: is it possible to cache hex_mat for a group of previously regularized vertices/hexes? 
            # this will save redundant find_hex calls when the user presses the Regularize Hex button repeatedly
            # ideally, we could look up the selected 'cubes' and 'hex_mat' from the cached data after a single Update Pickle Hexes
            # as many times as we want for any legal seletion without doing a find_hex each time
        
        if not cubes:
            print("Error: No cubes selected!")
            return {'FINISHED'}
        
        # 'priorities' is an internal switch that specifies whether the
        # regularization is broken up into pieces by topology region. 
        # The default line "[np.unique(hex_mat).tolist()]" returns a list 
        # of lists with one member list.
        # That member list contains all of the topology numbers that were found. 
        # Consequently, the default behavior is to consider all priority 
        # regions together, and the for-loop executes only once.
        # see the method HexInterpSubdivide() for more details
        priorities = [np.unique(hex_mat).tolist()]
        #priorities = [[0,1],[2,3]]
            
        # flip the priorities list so that the last list executes first, 
        # and the first list executes last
        priorities.reverse()
        for priority_group in priorities:
            currHexesInMatlGrp = np.array([ind for ind,val in enumerate(hex_mat) 
                if val in priority_group])

            cubes_group = np.array(cubes).copy()
            cubes_group = cubes_group[np.array(currHexesInMatlGrp),:].tolist()
            
            allVertInds = np.unique(np.array(cubes_group[:]).flatten())

            isVertIncluded = np.array([vert in list(allVertInds) 
                for vert in range(len(mesh.vertices))])

            vertSelectionMap = np.cumsum(isVertIncluded)-1

            # initialize verts_array of selected vertex coordinates
            verts_array = np.zeros([np.shape(allVertInds)[0],3])

            for ind,val in enumerate(verts_array):
                verts_array[ind,0] = mesh.vertices[allVertInds[ind]].co.x
                verts_array[ind,1] = mesh.vertices[allVertInds[ind]].co.y
                verts_array[ind,2] = mesh.vertices[allVertInds[ind]].co.z

            # cubes_reduced defines the hex connectivity by (presumably) Continuity convention
            cubes_reduced = vertSelectionMap[np.array(cubes_group)]

            # cubes_trans defines the hex connectivity by Blender convention
            cubes_trans = np.zeros(np.shape(cubes_reduced),'int32')
            for ind,cube in enumerate(cubes_reduced):
                cubes_trans[ind,0] = cube[0]
                cubes_trans[ind,1] = cube[1]
                cubes_trans[ind,2] = cube[3]
                cubes_trans[ind,3] = cube[2]
                cubes_trans[ind,4] = cube[4]
                cubes_trans[ind,5] = cube[5]
                cubes_trans[ind,6] = cube[7]
                cubes_trans[ind,7] = cube[6]
            nN = regularize_elements(cubes_trans,         # e in regularize_elements.py
                                     verts_array,         # n in regularize_elements.py
                                     itt,                 # itt in '''
                                     preserveRidges=preserveRidges,
                                     immobExtraord=immobExtraord,
                                     tangentMotion=tangentMotion,
                                     normalMotion=normalMotion,
                                     internalMotion=internalMotion)
            # import pdb; pdb.set_trace()
            # nN is exactly the same as verts_array...
            for ind,vert_num in enumerate(allVertInds):
                b_verts[vert_num].co.x = nN[ind,0]
                b_verts[vert_num].co.y = nN[ind,1]
                b_verts[vert_num].co.z = nN[ind,2]

        print("Regularization completed!")

        b_mesh.to_mesh(mesh)
        hexblender_helpers.reset_mode(orig_mode)

        return {'FINISHED'}


class HEXBLENDER_OT_regularize_sfc(bpy.types.Operator):
    """Regularize surface"""

    bl_idname = "hexblender.regularize_sfc"
    bl_label = "Regularize Surface"
    bl_description = ("Regularize Surfaces")
    bl_options = {'REGISTER'}

    @classmethod
    def poll(cls, context):
        return context.object is not None

    def execute(self, context):
        mesh, orig_mode = hexblender_helpers.get_active_mesh()

        b_mesh = bmesh.new()
        b_mesh.from_mesh(mesh)

        verts = b_mesh.verts
        faces = b_mesh.faces

        try:
            hex_scene = context.scene.hexblender
        except Exception as msg:
            print("ERROR: Unable to get hexblender scene!  %s" % msg)
            return {'FINISHED'}
        
        itt = hex_scene.hexblender_properties.regularize_sfc_iters
        
        verts.ensure_lookup_table()
        selected_verts_obj = [v for v in verts if v.select]
        if len(selected_verts_obj) == 0:
            print("ERROR: nothing selected!")
            return {'FINISHED'}

        matList = hexblender_helpers.get_material_data(b_mesh)
        tmp, f, n = hexblender_helpers.get_boundary_data([], mesh=mesh)
        f = np.array(f)
        n = np.array(n)
        
        # we need the id lists to more readily conform to existing code
        selected_faces_obj = [p for p in faces if p.select]
        selected_faces_id = [f.index for f in selected_faces_obj]
        selected_vert_id = [v.index for v in selected_verts_obj]

        # create mapping back to original verts.  should also check to 
        # make sure the domain is connected
        # first get rid of vertices that don't belong to a face
        fSelected = np.array(selected_faces_id)
        f_verts_Selected = f[fSelected,:].copy().flatten()
        nSelected = np.array([vert for vert in selected_vert_id if vert in f_verts_Selected])
        nboolSelected = np.array([vert in list(nSelected) for vert in range(np.shape(n)[0])])
        vertSel_map = np.cumsum(nboolSelected)-1

        n_reg = n[nSelected,:].copy()
        f_reg = vertSel_map[f[fSelected,:]]

        nF,nN = regularize_surface(f_reg,n_reg,itt)
        
        #change the verts that have been modified
        n[nSelected,:] = nN
        
        #modify the base object in Blender, "mesh.verts"
        for ind, vert_num in enumerate(nSelected):
            verts[vert_num].co.x = nN[ind,0]
            verts[vert_num].co.y = nN[ind,1]
            verts[vert_num].co.z = nN[ind,2]
        
        print("Regularization completed!")
        b_mesh.to_mesh(mesh)
        hexblender_helpers.reset_mode(orig_mode)

        return {'FINISHED'}


class HEXBLENDER_OT_edge_lengths(bpy.types.Operator):
    """Get length of selected edges"""

    bl_idname = "hexblender.edge_lengths"
    bl_label = "Edge Lengths"
    bl_description = ("Edge Lengths")
    bl_options = {'REGISTER'}

    @classmethod
    def poll(cls, context):
        return context.object is not None

    def execute(self, context):
        # We have to be in edit more for this to work as
        # the easiest way I have seen to get the edge length is
        # by using bmesh
        if context.object.mode != 'EDIT':
            print("NOTICE: Edge Lengths command needs to run in EDIT MODE!")
            print("NOTICE: Changing to EDIT mode...")
            bpy.ops.object.mode_set(mode='EDIT')

        mesh, orig_mode = hexblender_helpers.get_active_mesh()

        b_mesh = bmesh.new()
        b_mesh.from_mesh(mesh)

        edges = b_mesh.edges

        sum_edges = 0.
        sqrt_sum = 0.
        length_list = []

        selected_edges = [e for e in edges if e.select]
        print("\n Number of edges selected: %s" % len(selected_edges))
        if len(selected_edges) == 1:
            print("\n Length of Edge:")
            print(selected_edges[0].calc_length())

        elif len(selected_edges) > 2:

            for edge in selected_edges:
                length_list.append(edge.calc_length())
                sum_edges = sum_edges + length_list[-1]

            average_edges = sum_edges / (len(selected_edges))

            for index in range(len(selected_edges)):
                extension = length_list[index] - average_edges
                extension = extension * extension
                sqrt_sum = sqrt_sum + extension

            stddev_edges = np.sqrt(sqrt_sum / (len(selected_edges)-1))
            length_list.sort()

            print("\n Standard deviation of edges length:")
            print(" -----------------------------------")
            print(" StdDev: %0.2f +/- %0.2f" % (average_edges, stddev_edges))
            print("  5%%: %.2f \t 95%%: %.2f" % (average_edges - (2*stddev_edges), average_edges + (2*stddev_edges)))
            print(" MIN: %.2f \t MAX: %.2f" % (length_list[0], length_list[len(selected_edges)-1]))

        else:
            print("ERROR: No edges selected")

        b_mesh.to_mesh(mesh)
        hexblender_helpers.reset_mode(orig_mode)
        return {'FINISHED'}


class HEXBLENDER_OT_print_data(bpy.types.Operator):
    """Print information of currently selected data"""

    bl_idname = "hexblender.print_data"
    bl_label = "Print Selection-Data"
    bl_description = ("Print info of selected data")
    bl_options = {'REGISTER'}

    @classmethod
    def poll(cls, context):
        return context.object is not None

    def execute(self, context):
        # We have to be in edit more for this to work as I'm using bmesh

        mesh, orig_mode = hexblender_helpers.get_active_mesh()

        faces = mesh.polygons
        edges = mesh.edges
        verts = mesh.vertices

        selected_faces_obj = [f for f in faces if f.select]
        selected_faces_index = [f.index for f in faces if f.select]
        selected_edges_index = [e.index for e in edges if e.select]
        selected_verts_index = [v.index for v in verts if v.select]

        if len(selected_verts_index) == 0:
            print("Nothing selected!")
            return {'FINISHED'}
        
        num_verts_in_face = len(faces[0].vertices)
        
        #cheap way to not call this when its a triangular mesh
        if not num_verts_in_face == 3:

            cubesAll, MatListAll = hexblender_helpers.get_cached_data(mesh)

            print("finding selected cubes...")
            cubes,HexMat = hexblender_helpers.find_hex(selectionOnly=True)

        print("\n Print Selection Data:")
        print(" ---------------------")

        print("\nFaces selected:")
        print(selected_faces_index)

        print("\nVerts indexed from 1\n -----")
        print((np.array(selected_verts_index)+1).tolist())
        
        if not num_verts_in_face==3:
            numberStr = ""
            numberList = []
            cubesAllDict = dict()
            for ind,cube in enumerate(cubesAll):
                currTuple = tuple(sorted(cube))
                cubesAllDict[currTuple] = ind
            for cube in cubes:
                cubeTuple = tuple(sorted(cube))
                if cubeTuple in cubesAllDict:
                    parentCubeIndex = cubesAllDict[cubeTuple]
                numberList.append(parentCubeIndex)
            numberList = sorted(numberList)
            for cubeNumber in numberList:
                numberStr += str(cubeNumber+1) + ", "
            
            #find the one-neighborhood of the selected cubes. This is principally for debugging hex derivatives
            oneNeighElems = []
            for cubeNumber in numberList:
                currCube = cubesAll[cubeNumber]
                for nodeNum in currCube: #indexed from zero
                    adjElements = np.nonzero(np.array(cubesAll) == nodeNum)[0]
                    oneNeighElems.append(adjElements.tolist())
            oneNeighElems = sum(oneNeighElems,[])
            oneNeighElems = list(set(oneNeighElems))
            for cubeNumber in numberList:
                oneNeighElems.remove(cubeNumber)
            
            print("\nOne neighborhood of selected cubes (indexed from 0) = ")
            print(sorted(oneNeighElems))
            
            print("\nSelected cubes' vertices: (cubes indexed from 1, vertices indexed from 0)")
            for index in numberList:
                print("%-6d: " %(index+1) , cubesAll[index])
            
            print("\nCube numbers (indexed from 1) = ")
            print(numberStr)
            print("(%d cube number(s) printed)" %len(numberList))
            
            # Print element-face pair for selected faces
            faceCubeList = []
            faceCubeFaceList = []
            for face in selected_faces_index:
                cubeIndex = 0
                faceDir = -1
                for cube in cubesAll:
                    nodesFound = []
                    numNodesFound = 0
                    for ind in range(num_verts_in_face):
                        for cubeNodeIndex,cubeNode in enumerate(cube):
                            if cubeNode == faces[face].vertices[ind]:
                                nodesFound.append(cubeNodeIndex+1)
                                numNodesFound += 1
                    if (numNodesFound == 4):
                        nodesFound.sort()
                        if (nodesFound == [1,3,5,7]):
                            faceDir = 1
                        elif (nodesFound == [2,4,6,8]):
                            faceDir = 2
                        elif (nodesFound == [1,2,5,6]):
                            faceDir = 3
                        elif (nodesFound == [3,4,7,8]):
                            faceDir = 4
                        elif (nodesFound == [1,2,3,4]):
                            faceDir = 5
                        elif (nodesFound == [5,6,7,8]):
                            faceDir = 6
                        break                
                    cubeIndex += 1
                faceCubeList.append(cubeIndex+1)
                faceCubeFaceList.append(faceDir)
            print(faceCubeList)
            print(faceCubeFaceList)
        
        hexblender_helpers.reset_mode(orig_mode)
        return {'FINISHED'}


class HEXBLENDER_OT_adjust_nodes_linear(bpy.types.Operator):
    """Adjust nodes in linear direction"""

    bl_idname = "hexblender.adjust_nodes_linear"
    bl_label = "Adjust nodes linear"
    bl_description = ("Adjust selected vertex positions to lie on line defined by bounding vertices")
    bl_options = {'REGISTER'}

    @classmethod
    def poll(cls, context):
        return context.object is not None

    def execute(self, context):

        mesh, orig_mode = hexblender_helpers.get_active_mesh()

        b_mesh = bmesh.new()
        b_mesh.from_mesh(mesh)

        if b_mesh is None:
            print("Do not have a bmesh object!")
            return {'FINISHED'}

        verts = b_mesh.verts
        verts.ensure_lookup_table()
        selected_verts_index = [v.index for v in verts if v.select]
        numSelectedVerts = len(selected_verts_index)
        maxDistance = 0
        v1Index = -1
        v2Index = -1
        for index1 in selected_verts_index:
            for index2 in selected_verts_index:
                v1 = verts[index1]
                v2 = verts[index2]
                v1Coords = [v1.co.x, v1.co.y, v1.co.z]
                v2Coords = [v2.co.x, v2.co.y, v2.co.z]
                currentDistance =   np.sqrt((v2Coords[0] - v1Coords[0])*(v2Coords[0] - v1Coords[0]) + 
                                            (v2Coords[1] - v1Coords[1])*(v2Coords[1] - v1Coords[1]) +
                                            (v2Coords[2] - v1Coords[2])*(v2Coords[2] - v1Coords[2]))
                if (currentDistance > maxDistance):
                    maxDistance = currentDistance
                    v1Index = index1
                    v2Index = index2
                    
        v1 = verts[v1Index]
        v2 = verts[v2Index]
        v1Coords = [v1.co.x, v1.co.y, v1.co.z]
        v2Coords = [v2.co.x, v2.co.y, v2.co.z]
        v1v2Vector = [(v2Coords[0] - v1Coords[0])/(maxDistance),
                      (v2Coords[1] - v1Coords[1])/(maxDistance),
                      (v2Coords[2] - v1Coords[2])/(maxDistance)]
        distanceArray = []
        for index in selected_verts_index:
            v = verts[index]
            vCoords = [v.co.x, v.co.y, v.co.z]
            v1Distance =   np.sqrt((vCoords[0] - v1Coords[0])*(vCoords[0] - v1Coords[0]) + 
                                   (vCoords[1] - v1Coords[1])*(vCoords[1] - v1Coords[1]) +
                                   (vCoords[2] - v1Coords[2])*(vCoords[2] - v1Coords[2]))
            distanceArray.append(tuple([index, v1Distance]))
            
        distanceArray.sort(key = lambda distance : distance[1] )

        for i in range(1,numSelectedVerts):
            index, distance = distanceArray[i]
            currentDistance = i*maxDistance/(numSelectedVerts-1)
            v = verts[index]
            v.co.x = v1Coords[0] + v1v2Vector[0]*currentDistance
            v.co.y = v1Coords[1] + v1v2Vector[1]*currentDistance
            v.co.z = v1Coords[2] + v1v2Vector[2]*currentDistance
            
        b_mesh.to_mesh(mesh)
        hexblender_helpers.reset_mode(orig_mode)
        return {'FINISHED'}


class HEXBLENDER_OT_adjust_nodes_circ(bpy.types.Operator):
    """Adjust nodes in linear direction"""

    bl_idname = "hexblender.adjust_nodes_circ"
    bl_label = "Circular nodal adjustment"
    bl_description = ("Adjust selected vertex positions to lie on circle with prescribed radius")
    bl_options = {'REGISTER'}

    @classmethod
    def poll(cls, context):
        return context.object is not None

    def execute(self, context):

        try:
            hex_scene = context.scene.hexblender
        except Exception as msg:
            print("ERROR: Unable to get hexblender scene!  %s" % msg)
            return {'FINISHED'}

        mesh, orig_mode = hexblender_helpers.get_active_mesh()

        b_mesh = bmesh.new()
        b_mesh.from_mesh(mesh)

        if b_mesh is None:
            print("Do not have a bmesh object!")
            return {'FINISHED'}

        verts = b_mesh.verts
        verts.ensure_lookup_table()
        selected_verts_index = [v.index for v in verts if v.select]
        numCircularVerts = len(selected_verts_index)
        sumDistance = 0
        sumCoords = [0, 0, 0]
        for index in selected_verts_index:
            v = verts[index]
            sumCoords[0] += v.co.x
            sumCoords[1] += v.co.y
            sumCoords[2] += v.co.z
        cCoords = [sumCoords[0]/float(numCircularVerts), sumCoords[1]/float(numCircularVerts), sumCoords[2]/float(numCircularVerts)] 
        
        circRadius = hex_scene.hexblender_properties.radius_of_nodal_adjustment
        avgDistance = circRadius
        if (circRadius == 0):
            for index in selected_verts_index:
                v = verts[index]
                vCoords = [v.co.x, v.co.y, v.co.z]
                currentDistance =   np.sqrt((vCoords[0] - cCoords[0])*(vCoords[0] - cCoords[0]) + 
                                            (vCoords[1] - cCoords[1])*(vCoords[1] - cCoords[1]) +
                                            (vCoords[2] - cCoords[2])*(vCoords[2] - cCoords[2]))
                sumDistance += currentDistance
            avgDistance = sumDistance/float(numCircularVerts)
        
        for index in selected_verts_index:
            v = verts[index]
            vCoords = [v.co.x, v.co.y, v.co.z]
            currentDistance =   np.sqrt((vCoords[0] - cCoords[0])*(vCoords[0] - cCoords[0]) + 
                                        (vCoords[1] - cCoords[1])*(vCoords[1] - cCoords[1]) +
                                        (vCoords[2] - cCoords[2])*(vCoords[2] - cCoords[2]))
            currentVector = [(vCoords[0] - cCoords[0])*(avgDistance/currentDistance),
                             (vCoords[1] - cCoords[1])*(avgDistance/currentDistance),
                             (vCoords[2] - cCoords[2])*(avgDistance/currentDistance)]
            v.co.x = cCoords[0] + currentVector[0]
            v.co.y = cCoords[1] + currentVector[1]
            v.co.z = cCoords[2] + currentVector[2]
        
        b_mesh.to_mesh(mesh)
        hexblender_helpers.reset_mode(orig_mode)
        return {'FINISHED'}


class HEXBLENDER_OT_hex_interp_subdivide(bpy.types.Operator):
    """Divide each hex element into 8 elements"""

    bl_idname = "hexblender.hex_interp_subdivide"
    bl_label = "Hex Interp Subdivide"
    bl_description = ("Perform interpolatory subdivision (Li-Bao, 2005) on mesh surfaces")
    bl_options = {'REGISTER'}

    @classmethod
    def poll(cls, context):
        return context.object is not None

    def execute(self, context):
        try:
            hex_scene = context.scene.hexblender
        except Exception as msg:
            print("ERROR: Unable to get hexblender scene!  %s" % msg)
            return {'FINISHED'}

        mesh, orig_mode = hexblender_helpers.get_active_mesh()


        if mesh is None:
            print("Do not have a mesh object!")
            return {'FINISHED'}

        cubes, matList = hexblender_helpers.get_cached_data(mesh)

        e,tmp,n = hexblender_helpers.get_boundary_data(cubes)
        e = np.array(e)
        n = np.array(n)
        
        # some changes to hexInterpSubdiv - now you can do a hexahedral 
        # interpolation on piecewise parts of
        # the mesh, demarcated by topology number (right now, the 
        # topology numbers must be manually placed in the source code). 
        # With these keyword arguments absent, it proceeds treating the 
        # mesh as one large region
        
        # right now, the keyword arguments are "matList", which is a 
        # list of element priority numbers note right now, 
        # contiguousTopologies() is not called.
        
        # 'priorities' is a list of lists. Each constituent list is a 
        # grouping of topology numbers that have 'equal' priority. For 
        # regions that interface two topology regions, the area with 
        # higher 'priority' determines where the vertices are placed in 
        # that connecting region. So if   priorities = [[0,1],[2,3]]
        # materials 0 and 1 are "priority group 0", and materials 
        # 2 and 3 are "priority group 1". At the location where priority 
        # group 0 and priority group 1 articulate, HexInterpSubdiv() 
        # will calculate new vertex locations for each region, but 
        # priority group 0's calculated vertex locations will be used.

        interpPriorities = hexblender_helpers.get_priorities(
            hex_scene.hexblender_properties.interp_priorities)

        #interpPriorities = [[0,4,5,6,7,8,1,2,3,9,10]]
        #interpPriorities = [[0,1,2,3,4,5,6,7,8,9,12,13,14,15],[10,11],[16]]
        #interpPriorities = [[0,1]]

        thinPlateMapping = hex_scene.hexblender_properties.thin_plate_mapping

        splinePriorities = hexblender_helpers.get_priorities(
            hex_scene.hexblender_properties.spline_priorities)
        
        #splinePriorities = [[0,4,5,6,7,8,1,2,3,9,10]]
        #splinePriorities = [[0],[1],[2],[3],[4],[5],[6],[7],[8],[9],[10]]
        #splinePriorities = [[0],[1],[2],[3],[4],[5],[6],[7]]
        #splinePriorities = [[0],[1],[2],[3],[4],[5],[6],[7],[8],[9],[10],[11],[12],[13],[14],[15],[16]]
        #splinePriorities = [[0],[1]]
        
        #eN, nN = HexInterpSubdiv(e,n)
        #eN, nN = HexInterpSubdiv(e,n, LinearSubdivOnly = True)
        #eN, nN = HexInterpSubdiv(e,n,MatList=matList,priorities=[[0,1],[2,3]])
        subdiv_type = int(hex_scene.hexblender_properties.subdiv_types)
        if subdiv_type == 0:
            print("Subdivide by the following priorities:")
            print("InterpPriorities: %s" % interpPriorities)
            print("SplinePriorities: %s" % splinePriorities)
            eN, nN = hex_interp_subdiv(e,
                                   n,
                                   MatList=matList,
                                   priorities=interpPriorities, 
                                   thinPlateMapping=thinPlateMapping, 
                                   thinPlateRegions=splinePriorities)
        elif subdiv_type == 1:
            print("Subdivide by element")
            eN, nN = hex_interp_subdiv_no_priorities(e,
                                   n,
                                   MatList=matList,
                                   priorities=None, 
                                   thinPlateMapping=thinPlateMapping, 
                                   thinPlateRegions=splinePriorities)
        else:
            print("!! Invalid subdivision type: %s" % subdiv_type)
            return

            
        # get existing mesh material info
        existing_mats = hexblender_helpers.get_mats_used_unused_sorted(mesh)
        new_cube_faces = hexblender_helpers.add_new_data(nN, [], [], eN, newObj=True)

        # We generally don't want to create a NEW object, but I'm trying
        # to get that route working first before trying to modify the
        # existing object.
        # This command is currently not undoable, so having a new object
        # (if memory allows) isn't that bad of an idea
        #new_cube_faces = hexblender_helpers.add_new_data(nN, [], [], eN)

        hexblender_helpers.set_mat_data_hex(matList, 
                                            new_cube_faces, 
                                            existing_mats = existing_mats)

        # execute these lines to pickle the subdivded mesh so you dont have 
        # to find it again
        # this will NOT preserve the element numbers that FindHex() would report
        HexMatNew = []
        for val in matList:
            HexMatNew.extend([val, val, val, val, val, val, val, val])

        cubes_normal, newfaces, newverts = hexblender_helpers.get_boundary_data(eN)
        HexMatNew2 = hexblender_helpers.contiguous_regions(cubes_normal, HexMatNew)
        #SaveTmpFile(elems=cubes_normal,elems_tr=eN,faces=newfaces,verts=newverts,mats=HexMatNew2)
        #Blender.Window.FileSelector(FileCallbackPickle, "Export Pickle File", Blender.sys.makename(ext='.pickle'))
        
        # I am perplexed as to why, but I have to turn cubes_normal into an 
        # array (then I turn it back into a list) 
        # or else the IDProperty type chokes only here and nowhere else
        context.scene.objects.active = context.scene.objects.active 
        hexblender_helpers.cache_data(np.array(cubes_normal).tolist(),HexMatNew2)

        if hex_scene.hexblender_properties.delete_orig_mesh:
            hexblender_helpers.delete_all(mesh.name)
        
        hexblender_helpers.reset_mode(orig_mode)
        return {'FINISHED'}


class HEXBLENDER_OT_update_pickle_hexs(bpy.types.Operator):
    """Find/Update Hexs"""

    bl_idname = "hexblender.update_pickle_hexs"
    bl_label = "Find/Update Hexs"
    bl_description = ("Find/update hexahedral elements, harmonize topologies, and cache results")
    bl_options = {'REGISTER'}

    # filepath = bpy.props.StringProperty(subtype = "FILE_PATH")

    @classmethod
    def poll(cls, context):
        return context.object is not None

    def execute(self, context):

        try:
            hex_scene = context.scene.hexblender
        except Exception as msg:
            print("ERROR: Unable to get hexblender scene!  %s" % msg)
            return {'FINISHED'}

        # recacheHexs=hex_scene.hexblender_properties.recacheHexs

        # fid = open(self.filepath,'wb')

        print("\n Search Hexs: \n")
        cubes, hexMat = hexblender_helpers.find_hex()
        hexMatNew = hexblender_helpers.contiguous_regions(cubes,hexMat)
        cubes_tpose, faces, verts = hexblender_helpers.get_boundary_data(cubes)

        # if recacheHexs:
        hexblender_helpers.cache_data(cubes, hexMatNew)

        # hexblender_helpers.write_pickle(fid, elems = cubes,
        #               elems_tr = cubes_tpose,
        #               faces = faces,
        #               verts = verts,
        #               mats = hexMatNew)

        # fid.close()

        return {'FINISHED'}

    # def invoke(self, context, event):
    #     context.window_manager.fileselect_add(self)
    #     return {'RUNNING_MODAL'}


class HEXBLENDER_OT_export_face_keys(bpy.types.Operator, ExportHelper):
    """Export Face Keys"""

    bl_idname = "hexblender.export_face_keys"
    bl_label = "Export Face Keys"
    bl_description = ("Export Face Keys")
    bl_options = {'REGISTER'}

    # ExportHelper mixin class uses this
    filename_ext = ".txt"

    filter_glob = bpy.props.StringProperty(
            default="*.txt",
            options={'HIDDEN'},
            )

    @classmethod
    def poll(cls, context):
        return context.object is not None

    def execute(self, context):
        return_code, elem_string = hexblender_helpers.write_out_2d_elements(self.filepath)

        return return_code

    def invoke(self, context, event):
        self.filepath = "face_elems.txt"
        context.window_manager.fileselect_add(self)
        return {'RUNNING_MODAL'}


class HEXBLENDER_OT_std_faces(bpy.types.Operator):
    """Standard Deviation of Faces"""

    bl_idname = "hexblender.std_faces"
    bl_label = "Standard Deviation of Faces"
    bl_description = ("Standard Deviation of Faces")
    bl_options = {'REGISTER'}

    @classmethod
    def poll(cls, context):
        return context.object is not None

    def execute(self, context):
        mesh, orig_mode = hexblender_helpers.get_active_mesh()

        faces = mesh.polygons

        sum_faces = 0.
        sqrt_sum = 0.
        area_list = []

        selected_faces_obj = [p for p in faces if p.select]
        num_selected = len(selected_faces_obj)

        if num_selected == 0:
            print(" ERROR: no faces selected")
        elif num_selected == 1:
            print("\n Area of selected face:")
            print(selected_faces_obj[0].area)
        else:

            for face in selected_faces_obj:
                sum_faces = sum_faces + face.area
                area_list.append(face.area)

            average_faces = sum_faces / num_selected

            for face in selected_faces_obj:
                extension = face.area - average_faces
                extension = extension * extension
                sqrt_sum = sqrt_sum + extension

            stddev_faces = np.sqrt(sqrt_sum / (num_selected - 1))
            area_list.sort()

            print("\n Standard deviation of faces area:")
            print(" ---------------------------------")
            print(" StdDev: %0.2f +/- %0.2f" % (average_faces, stddev_faces))
            print("  5%%: %.2f \t 95%%: %.2f" % (average_faces - (2*stddev_faces), average_faces + (2*stddev_faces)))
            print(" MIN: %.2f \t MAX: %.2f" % (area_list[0], area_list[num_selected-1]))

        hexblender_helpers.reset_mode(orig_mode)
        return {'FINISHED'}       


class HEXBLENDER_OT_std_angles(bpy.types.Operator):
    """Standard Deviation of Angles"""

    bl_idname = "hexblender.std_angles"
    bl_label = "Standard Deviation of Angles"
    bl_description = ("Standard Deviation of Angles")
    bl_options = {'REGISTER'}

    @classmethod
    def poll(cls, context):
        return context.object is not None

    def execute(self, context):
        mesh, orig_mode = hexblender_helpers.get_active_mesh()

        edges = mesh.edges
        faces = mesh.polygons
        verts = mesh.vertices

        selected_faces_obj = [p for p in faces if p.select]
        num_selected = len(selected_faces_obj)

        if num_selected == 0:
            print(" ERROR: no faces selected")
        else:
            neighbor_list = []
            for i in range(len(faces)):
                for j in range(4):
                    current_couple = []
                    current_couple.append(faces[i].edge_keys[j])
                    current_couple.append(faces[i].edge_keys[(j+1)%4])
                    neighbor_list.append(current_couple)

            for i in range(len(neighbor_list)):
                edges_index = []
                for j in range(2):
                    edge = neighbor_list[i][j]
                    for k in range(len(edges)):
                        if edge == edges[k].key:
                            neighbor_list[i][j] = k
                            break

            edges_vec = []
            for edge in edges:
                edge_vec = [0,0,0]
                edge_vec[0] = verts[edge.vertices[1]].co.x -\
                              verts[edge.vertices[0]].co.x
                edge_vec[1] = verts[edge.vertices[1]].co.y -\
                              verts[edge.vertices[0]].co.y
                edge_vec[2] = verts[edge.vertices[1]].co.z -\
                              verts[edge.vertices[0]].co.z
                edges_vec.append(edge_vec)

            angles = []
            for i in range(len(neighbor_list)):
                vec1 = Vector(edges_vec[neighbor_list[i][0]])
                vec2 = Vector(edges_vec[neighbor_list[i][1]])
                # Convert to degrees
                angle = vec1.angle(vec2) * (180/pi)
                angles.append(angle)

            sum_angles = 0.
            sqrt_sum = 0.
            area_list = []

            for angle in angles:
                sum_angles = sum_angles + angle

            average_angles = sum_angles / (len(angles))

            for angle in angles:
                extension = angle - average_angles
                extension = extension * extension
                sqrt_sum = sqrt_sum + extension

            stddev_angles = np.sqrt(sqrt_sum / (len(angles)-1))
            angles.sort()

            print("\n Standard deviation of angles:")
            print(" -----------------------------")
            print(" StdDev: %0.2f +/- %0.2f" % (average_angles, stddev_angles))
            print("  5%%: %.2f \t 95%%: %.2f" % (average_angles - (2*stddev_angles), average_angles + (2*stddev_angles)))
            print(" MIN: %.2f \t MAX: %.2f" % (angles[0], angles[len(angles)-1]))

        hexblender_helpers.reset_mode(orig_mode)
        return {'FINISHED'}       

class HEXBLENDER_OT_select_vertex(bpy.types.Operator):
    """Select Vertex"""

    bl_idname = "hexblender.select_vertex"
    bl_label = "Select Vertex"
    bl_description = ("Select Vertex")
    bl_options = {'REGISTER'}

    @classmethod
    def poll(cls, context):
        return context.object is not None

    def execute(self, context):
        mesh, orig_mode = hexblender_helpers.get_active_mesh()

        verts = mesh.vertices

        try:
            hex_scene = context.scene.hexblender
        except Exception as msg:
            print("ERROR: Unable to get hexblender scene!  %s" % msg)
            return

        vert_index_string = hex_scene.hexblender_properties.vert_index_number
        vert_index = hexblender_helpers.parse_list(vert_index_string)

        if max(vert_index) > len(verts) or min(vert_index) <= 0:
            print("ERROR: Vertex index is out of range, valid range: 1 to %d" % len(verts))
        else:
            for v in vert_index:
                # We are showing vert indices one-based, but need to 
                # look them up using zero-based values
                verts[int(v-1)].select = True

            if len(vert_index) == 1:
                # only move the cursor if we are only selecting 1 vert
                bpy.context.scene.cursor_location = verts[vert_index[0]-1].co

        hexblender_helpers.reset_mode(orig_mode)
        return {'FINISHED'}


class HEXBLENDER_OT_select_edge(bpy.types.Operator):
    """Select Edge"""

    bl_idname = "hexblender.select_edge"
    bl_label = "Select Edge"
    bl_description = ("Select Edge")
    bl_options = {'REGISTER'}

    @classmethod
    def poll(cls, context):
        return context.object is not None

    def execute(self, context):
        mesh, orig_mode = hexblender_helpers.get_active_mesh()

        edges = mesh.edges
        verts = mesh.vertices

        try:
            hex_scene = context.scene.hexblender
        except Exception as msg:
            print("ERROR: Unable to get hexblender scene!  %s" % msg)
            return

        edge_index_string = hex_scene.hexblender_properties.edge_index_number
        edge_index = hexblender_helpers.parse_list(edge_index_string)

        if max(edge_index) > len(edges) or min(edge_index) <= 0:
            print("ERROR: Edge index is out of range, valid range 1 to %d" % len(edges))
        else:
            for e in edge_index:
                edges[e-1].select = True

                if len(edge_index) == 1:
                    edge_co1 = verts[edges[e-1].vertices[0]].co
                    edge_co2 = verts[edges[e-1].vertices[1]].co
                    edge_co = []
                    for i in range(3):
                        edge_co.append((edge_co1[i]+edge_co2[i])/2.)

                    bpy.context.scene.cursor_location = edge_co

        hexblender_helpers.reset_mode(orig_mode)
        return {'FINISHED'}


class HEXBLENDER_OT_select_face(bpy.types.Operator):
    """Select Face"""

    bl_idname = "hexblender.select_face"
    bl_label = "Select Face"
    bl_description = ("Select Face")
    bl_options = {'REGISTER'}

    @classmethod
    def poll(cls, context):
        return context.object is not None

    def execute(self, context):
        mesh, orig_mode = hexblender_helpers.get_active_mesh()

        faces = mesh.polygons

        try:
            hex_scene = context.scene.hexblender
        except Exception as msg:
            print("ERROR: Unable to get hexblender scene!  %s" % msg)
            return

        face_index_string = hex_scene.hexblender_properties.face_index_number
        face_index = hexblender_helpers.parse_list(face_index_string)

        if max(face_index) > len(faces) or min(face_index) <= 0:
            print("ERROR: Face index is out of range, valid range 1 to : %d" % len(faces))
        else:
            for f in face_index:
                faces[f-1].select = True

            if len(face_index) == 1:
                bpy.context.scene.cursor_location = faces[f-1].center

        hexblender_helpers.reset_mode(orig_mode)
        return {'FINISHED'}              

class HEXBLENDER_OT_select_hex(bpy.types.Operator):
    """Select Hex Element"""

    bl_idname = "hexblender.select_hex"
    bl_label = "Select Hex"
    bl_description = ("Select Hex Element")
    bl_options = {'REGISTER'}

    @classmethod
    def poll(cls, context):
        return context.object is not None

    def execute(self, context):
        mesh, orig_mode = hexblender_helpers.get_active_mesh()
        verts = mesh.vertices

        try:
            hex_scene = context.scene.hexblender
        except Exception as msg:
            print("ERROR: Unable to get hexblender scene!  %s" % msg)
            return

        hex_index_string = hex_scene.hexblender_properties.hex_index_number
        hex_index = hexblender_helpers.parse_list(hex_index_string)
        bpy.types.Mesh.cached_data_json = bpy.props.StringProperty()

        # see if we have already cached the data
        if len(mesh.cached_data_json) > 0:
            try:
                cubes, mat_list = hexblender_helpers.reload_hex(json.loads(mesh.cached_data_json))
                print("Reloaded stored hex data!")
            except Exception as msg:
                print("Error reloading hex data: %s" % msg)
        else:
            print("finding cubes...")
            cubes, mat_list = hexblender_helpers.find_hex(selectionOnly=False)

        # Lets make sure we're in vertex select mode
        # Vertex, Edge, Face
        context.tool_settings.mesh_select_mode = (True, False, False)

        cubes = np.array(cubes)
        max_elem = len(cubes)
        if max(hex_index) > (len(cubes)) or min(hex_index) <= 0:
            print("ERROR: hex_index is out of range. Valid range 1 to %d"  % len(cubes))
        else:
            for h in hex_index:
                cached_data = json.loads(mesh.cached_data_json)
                cubeVerts = cached_data[str(h-1)]
                for vertInd in cubeVerts:
                    verts[vertInd].select = True

        hexblender_helpers.reset_mode(orig_mode)
        return {'FINISHED'}   


class HEXBLENDER_OT_rotate_hexs(bpy.types.Operator):
    """Rotate Hexs"""

    bl_idname = "hexblender.rotate_hexs"
    bl_label = "Rotate Hexs"
    bl_description = ("Rotate orientation of local coordinate frame (i.e. xi1, xi2, xi3 directions defined by order of element vertices) for all hexes in a topology/material region")
    bl_options = {'REGISTER'}

    filepath = bpy.props.StringProperty(subtype = "FILE_PATH")

    @classmethod
    def poll(cls, context):
        return context.object is not None

    def execute(self, context):
        mesh, orig_mode = hexblender_helpers.get_active_mesh()

        try:
            hex_scene = context.scene.hexblender
        except Exception as msg:
            print("ERROR: Unable to get hexblender scene!  %s" % msg)
            return

        rotate_type = int(hex_scene.hexblender_properties.rotate_types)
        rotate_matl_index = hex_scene.hexblender_properties.rotate_matl_index

        cubes, mat_list = hexblender_helpers.get_cached_data(mesh)

        current_hex_numbers = np.where(np.array(mat_list) == rotate_matl_index)[0]


        #DEFINING TYPES OF ROTATIONS
        rotationIndices =       [[0,1,2,3,4,5,6,7], #  0: do nothing
                                 [0,4,1,5,2,6,3,7], #  1: keep origin at nodal position 1,           
                                                    #     turn xi3 into xi1,                
                                                    #     xi1 into xi2,                     
                                                    #     xi2 into xi3 
                                 [1,3,0,2,5,7,4,6], #  2: move origin to nodal position 2,           
                                                    #     turn xi2 into xi1,                
                                                    #     turn xi1 into xi2 and flip,       
                                                    #     preserve xi3
                                 [1,5,3,7,0,4,2,6], #  3: move origin to nodal position 2,           
                                                    #     turn xi3 into xi1,                
                                                    #     preserve xi2,                     
                                                    #     turn xi1 and xi3 and flip
                                 [2,0,3,1,6,4,7,5], #  4: move origin to nodal position 3,           
                                                    #     turn xi2 into xi1 and flip,       
                                                    #     turn xi1 into xi2,                
                                                    #     preserve xi3
                                 [2,6,0,4,3,7,1,5], #  5: move origin to nodal position 3,           
                                                    #     turn xi3 into xi1,                
                                                    #     flip xi2,                         
                                                    #     turn xi2 into xi3 
                                 [3,2,1,0,7,6,5,4], #  6: move origin to nodal position 4,           
                                                    #     flip xi1,                         
                                                    #     flip xi2,                         
                                                    #     preserve xi3
                                 [3,7,2,6,1,5,0,4], #  7: move origin to nodal position 4,           
                                                    #     turn xi3 into xi1,                
                                                    #     turn xi1 into xi2 and flip,       
                                                    #     turn xi2 into xi3 and flip
                                 [4,6,5,7,0,2,1,3], #  8: move origin to nodal position 5,           
                                                    #     turn xi2 into xi1,                
                                                    #     turn xi1 into xi2,                
                                                    #     flip xi3
                                 [4,0,6,2,5,1,7,3], #  9: move origin to nodal position 5,           
                                                    #     turn xi3 into xi1 and flip,       
                                                    #     preserve xi2,                     
                                                    #     turn xi1 into xi3
                                 [5,4,7,6,1,0,3,2], #  10: move origin to nodal position 6,          
                                                    #      flip xi1,                         
                                                    #      keep xi2,                         
                                                    #      flip xi3
                                 [5,1,4,0,7,3,6,2], #  11: move origin to nodal position 6,          
                                                    #      turn xi3 into xi1 and flip,       
                                                    #      turn xi1 into xi2 and flip,       
                                                    #      turn xi2 into xi3
                                 [6,7,4,5,2,3,0,1], #  12: move origin to nodal position 7,          
                                                    #      preserve xi1,                     
                                                    #      flip xi2,                         
                                                    #      flip xi3
                                 [6,2,7,3,4,0,5,1], #  13: move origin to nodal position 7,          
                                                    #      turn xi3 into xi1 and flip,       
                                                    #      turn xi1 into xi2,                
                                                    #      turn xi2 into xi3 and flip
                                 [7,5,6,4,3,1,2,0], #  14: move origin to nodal position 8,          
                                                    #      turn xi2 into xi1 and flip,       
                                                    #      turn xi1 into xi2 and flip,       
                                                    #      flip xi3
                                 [7,3,5,1,6,2,4,0]] #  15: move origin to nodal position 8,          
                                                    #      turn xi3 into xi1 and flip,       
                                                    #      flip xi2,                         
                                                    #      turn xi1 into xi3 and flip
                                 
        currRotationIndices = np.array(rotationIndices[rotate_type])
        
        cached_data = json.loads(mesh.cached_data_json)
        for cubeNum in current_hex_numbers:
            #currHex = np.array(mesh.properties.pop(str(cubeNum)))
            currHex = np.array(cached_data.pop(str(cubeNum)))
            cached_data[str(cubeNum)] = currHex[currRotationIndices].tolist()

        mesh.cached_data_json = json.dumps(cached_data)
        print("\nDone rotating hexes")
        return {'FINISHED'}


class HEXBLENDER_OT_export_hermite_bicubic_derivs(bpy.types.Operator, ExportHelper):
    """Export hermite bicubic derivs"""

    bl_idname = "hexblender.export_hermite_bicubic_derivs"
    bl_label = "Export Hermite Bicubic Derivatives"
    bl_description = ("Export Hermite Bicubic Derivatives in Continuity format.")
    bl_options = {'REGISTER'}

    # ExportHelper mixin class uses this
    filename_ext = ".txt"

    filter_glob = bpy.props.StringProperty(
            default="*.txt",
            options={'HIDDEN'},
            )

    @classmethod
    def poll(cls, context):
        return context.object is not None

    def execute(self, context):
        return_code, derivs = hexblender_helpers.compute_write_hermite_2d_derivs(context, self.filepath)
        return return_code

    def invoke(self, context, event):
        self.filepath = "bicubic_derivs.txt"
        context.window_manager.fileselect_add(self)
        return {'RUNNING_MODAL'}


class HEXBLENDER_OT_export_fit_data_points(bpy.types.Operator, ExportHelper):
    """Export Fit Data Points"""

    bl_idname = "hexblender.export_fit_data_points"
    bl_label = "Export Fit Data Points"
    bl_description = ("Export Fit Data Points")
    bl_options = {'REGISTER'}

    # ExportHelper mixin class uses this
    filename_ext = ".txt"

    filter_glob = bpy.props.StringProperty(
            default="*.txt",
            options={'HIDDEN'},
            )

    @classmethod
    def poll(cls, context):
        return context.object is not None

    def execute(self, context):
        f = open(self.filepath, 'w')
        try:
            mesh, orig_mode = hexblender_helpers.get_active_mesh()
            verts = mesh.vertices
        except Exception as msg:
            print("Error getting vertices: %" % msg)
            return {'CANCELLED'}
        
        # get face materials for vertex materials. I am hoping when a 
        # vertex is shared by more than 1 material, it goes to the 
        # lowest one, but never checked...
        # the line that should be doing that is the "break" because it 
        # cycles through low numbers first
        
        elemMaterials = []
        currFaces = [[] for val in range(len(mesh.polygons))]
        for ind,face in enumerate(mesh.polygons):
            elemMaterials.append(face.material_index)
            for ind2,vert in enumerate(face.vertices):
                currFaces[ind].append(face.vertices[ind2])
        #~ print "currFaces = \n", currFaces
        #~ print "curr elem materials =     \n", elemMaterials
        # I THINK NEEDS MATERIAL NUMBERS INDEXED FROM ZERO FOR NOW, 
        # HAVENT TESTED OTHERWISE BUT TRIED TO GENERALIZE
        num_materials = len(np.unique(elemMaterials))
        mat_groupings = [[] for val in np.unique(elemMaterials)]
        
        for ind, mat_num in enumerate(np.unique(elemMaterials)):
             face_indices = np.nonzero(mat_num==elemMaterials)[0]
            #~ print "face  indices = ", face_indices
             if type(face_indices) == int: #catch single int returned, put in list
                 face_indices = [face_indices]
             for face_index in face_indices:
                mat_groupings[ind].extend(currFaces[face_index])
                mat_groupings[ind] = list(set(mat_groupings[ind]))
        #~ print "unique face mats = \n", np.unique(elemMaterials)
        #~ print "mat-groupings = \n", mat_groupings
        elemMaterialsList = np.unique(elemMaterials).tolist()
        #~ print "elem mats list = ", elemMaterialsList
        vertMaterials = []
        for vert in range(len(verts)):
            for ind, material in enumerate(elemMaterialsList):
                if vert in mat_groupings[ind]:
                    vertMaterials.append('group'+str(material+1))
                    break

        # === Header ===
        f.write('coord_1_val\tcoord_1_weight\tcoord_2_val\tcoord_2_weight\tcoord_3_val\tcoord_3_weight\tlabel_val\tData\n')

        # === Vertex List ===
        for i, v in enumerate(verts):
            f.write('%.6f\t1\t%.6f\t1\t%.6f\t1\t' % tuple(v.co))
            f.write('%s\t' % vertMaterials[i])
            f.write('%d\n' %(i+1))

        f.close()
        print("Successfully exported %s" % self.filepath)
        hexblender_helpers.reset_mode(orig_mode)
        return {'FINISHED'}

    def invoke(self, context, event):
        self.filepath = "fitting_data.txt"
        context.window_manager.fileselect_add(self)
        return {'RUNNING_MODAL'}

class HEXBLENDER_OT_export_bundled_3d_mesh(bpy.types.Operator, ExportHelper):
    """Export verts, elements and derivs for a tricubic hermite mesh"""

    bl_idname = "hexblender.export_bundled_3d_mesh"
    bl_label = "Export 3D verts, elems and derivatives"
    bl_description = ("Will create a pickle file of the vertices, elements and after subdividing, the derivatives for a tricubic hermite mesh.")
    bl_options = {'REGISTER'}

    # ExportHelper mixin class uses this
    filename_ext = ".pickle"

    filter_glob = bpy.props.StringProperty(
            default="*.pickle",
            options={'HIDDEN'},
            )

    @classmethod
    def poll(cls, context):
        return context.object is not None

    def execute(self, context):
        return_code, vert_string = hexblender_helpers.write_out_vertices()

        if 'FINISHED' in return_code:
            return_code, elem_string = hexblender_helpers.write_out_3d_elements()

        if 'FINISHED' in return_code:
            return_code, derivs = hexblender_helpers.compute_write_hermite_3d_derivs(context)

            # write data to pickle file
            data = [[vert_string],[elem_string],[derivs]]
            pickle.dump( data, open(self.filepath, "wb"), protocol = 2)

            print("Finished exporting %s" % self.filepath) 
            return {'FINISHED'}       
        else:
            return return_code

    def invoke(self, context, event):
        self.filepath = "bundled_3d_mesh.pickle"
        context.window_manager.fileselect_add(self)
        return {'RUNNING_MODAL'}


class HEXBLENDER_OT_export_bundled_2d_mesh(bpy.types.Operator, ExportHelper):
    """Export verts, elements and derivs for a bicubic hermite mesh"""

    bl_idname = "hexblender.export_bundled_2d_mesh"
    bl_label = "Export All components of a bicubic hermite mesh"
    bl_description = ("Will export the vertices, elements and derivatives for a bicubic hermite mesh into a file that Continuity can import.")
    bl_options = {'REGISTER'}

    # ExportHelper mixin class uses this
    filename_ext = ".pickle"

    filter_glob = bpy.props.StringProperty(
            default="*.pickle",
            options={'HIDDEN'},
            )

    @classmethod
    def poll(cls, context):
        return context.object is not None

    def execute(self, context):
        return_code, vert_string = hexblender_helpers.write_out_vertices()

        if 'FINISHED' in return_code:
            return_code, elem_string = hexblender_helpers.write_out_2d_elements()

        if 'FINISHED' in return_code:
            return_code, derivs = hexblender_helpers.compute_write_hermite_2d_derivs(context)

            # write data to pickle file
            data = [[vert_string],[elem_string],[derivs]]
            pickle.dump( data, open(self.filepath, "wb"), protocol = 2)

            print("Finished exporting %s" % self.filepath) 
            return {'FINISHED'}       
        else:
            return return_code

    def invoke(self, context, event):
        self.filepath = "bundled_2d_mesh.pickle"
        context.window_manager.fileselect_add(self)
        return {'RUNNING_MODAL'}

class HEXBLENDER_OT_set_tensor_output_file(bpy.types.Operator):
    '''
    Sets output filelocation
    '''
    bl_idname = "hexblender.set_tensor_output_file"
    bl_label = "Set Tensor Data Output Location"
    bl_description = ("Set Tensor Data Output File Location")
    bl_options = {'REGISTER'}

    filepath = bpy.props.StringProperty(subtype='FILE_PATH', default="")

    def execute(self, context):
        try:
            hex_scene = context.scene.hexblender
        except Exception as msg:
            print("ERROR: Unable to get hexblender scene!  %s" % msg)
            return {'FINISHED'}

        hex_scene.hexblender_properties.tensor_output_file = self.filepath
        return {'FINISHED'}

    def invoke(self, context, event):
        context.window_manager.fileselect_add(self)
        return {'RUNNING_MODAL'}

class HEXBLENDER_OT_get_tensor_input_file(bpy.types.Operator):
    '''
    Sets input filelocation
    '''
    bl_idname = "hexblender.get_tensor_input_file"
    bl_label = "Set Tensor Data Input Location"
    bl_description = ("Set Tensor Data Input File Location")
    bl_options = {'REGISTER'}

    filepath = bpy.props.StringProperty(subtype='FILE_PATH', default="")

    def execute(self, context):
        try:
            hex_scene = context.scene.hexblender
        except Exception as msg:
            print("ERROR: Unable to get hexblender scene!  %s" % msg)
            return {'FINISHED'}

        hex_scene.hexblender_properties.tensor_input_file = self.filepath
        return {'FINISHED'}

    def invoke(self, context, event):
        context.window_manager.fileselect_add(self)
        return {'RUNNING_MODAL'}


class HEXBLENDER_OT_compute_tensor_xform(bpy.types.Operator):
    '''
    This will compute, apply and save the transformed tensor data
    '''
    bl_idname = "hexblender.compute_tensor_xform"
    bl_label = "Compute and apply the transformation matrix"
    bl_description = ("Comptue and apply the transformation matrix")
    bl_options = {'REGISTER'}

    filepath = bpy.props.StringProperty(subtype='FILE_PATH', default="")

    def execute(self, context):
        try:
            hex_scene = context.scene.hexblender
        except Exception as msg:
            print("ERROR: Unable to get hexblender scene!  %s" % msg)
            return {'FINISHED'}

        aligned_obj_name = hex_scene.hexblender_properties.aligned_object
        unaligned_obj_name = hex_scene.hexblender_properties.unaligned_object

        aligned_data = hexblender_helpers.get_object_data(aligned_obj_name)
        unaligned_data = hexblender_helpers.get_object_data(unaligned_obj_name)

        if aligned_data == None or unaligned_data == None:
            return{'FINISHED'}

        print("Reading input data...")
        dt, coords, es, vs = hexblender_helpers.read_tensor_data_file(
                              hex_scene.hexblender_properties.tensor_input_file)
        print("Got it!")

        # Defines the coordinates to the transformed data that fits to
        # Continuity standards
        coords_seg = unaligned_data

        # number of data point
        num_data = len(aligned_data)

        # compute 4x4 transformation
        transformation1 = np.dot(unaligned_data.conj().T, unaligned_data)
        transformation2 = np.dot(unaligned_data.conj().T, aligned_data)
        transformation = np.linalg.solve(transformation1, transformation2)
        print("Computed transformation:")
        print(transformation)

        C = transformation[0:3, 0:3]

        # Apply transformation to scanner coordinates
        coords_temp = coords_seg
        coords_temp[:,3] = 1
        coords_r = np.dot(coords_temp, transformation)
        coords_r = np.delete(coords_r, 3, axis = 1)

        print('Applied transformation')

        ####################################################################
        # rotated diffusion tensors
        dt_r = np.zeros(dt.shape)
        # rotated diffusion tensors, log transforms
        dt_rl = np.zeros(dt.shape)
        dt_rscale = np.zeros(dt.shape)
        # roated diffusion tensors eigenvectors
        vpr = np.zeros(dt.shape)

        indsrem = []
        k = 0

        # Apply transformation to DT data
        for index in range(num_data):
            if index % 10000 == 0:
                print("%d tensors processed" % index)

            # DT rotation
            a = np.dot(C.conj().T, dt[:,:,index])
            dt_r[:,:,index] = np.dot(a, C)

            # finds the eigenvectors (v) and eigenvalues (e) of the DT data
            v, e = hexblender_helpers.eig_s(dt_r[:,:,index])

            # Eigenvectors
            vpr[:,:,index] = v

            # Eigenvalues
            # don't need epr(j,:)  = es_scale(j,:);
            # don't need DT_rscale(:,:,j) = vpr(:,:,j)*diag(epr(j,:))*vpr(:,:,j)';
            dt_rscale[:,:,index] = np.dot(vpr[:,:,index], np.dot(np.diag(es[index,:]), vpr[:,:,index].conj().T))

            # Compute matrix logarithm of DT
            #dt_rl[:,:,index] = logm(dt_r[:,:,index])
            dt_rl[:,:,index] = logm(dt_rscale[:,:,index])
            if np.isnan(dt_rl).any() or np.isinf(dt_rl).any():
                # mark data indices with nonreal tensors
                indsrem.append(index);


        #################################################################### 
        # Finds the index of the original coordinates at a certain slice in 
        # the model.  Because of signifcant figures, must first round the 
        # coordinates before finding slice at which to render 

        # find median image slice of DT data in scanner reference frametensors.
        median_image = np.median(coords[:,1])
        ind_slice = np.where(coords[:,1] == median_image)

        # initialize matrix to track DT data indices after removing bad voxels
        inds_slice2 = np.zeros((num_data)) 

        # DT data indices belonging to median slice get set to 1
        inds_slice2[ind_slice] = 1 

        # DT data indices belonging to bad voxels get reset back to 0
        inds_slice2[indsrem] = 0  

        # find DT data indices that are in median slice AND are not bad voxels
        ind_slice3 = np.where(inds_slice2==1)

        inds_dat = np.ones((num_data))
        inds_dat[indsrem] = 0
        inds_dat2 = np.where(inds_dat==1)

        #################################################################### 
        # save all DT data to Continuity data form
        hexblender_helpers.data_dt_form(hex_scene.hexblender_properties.tensor_output_file,
                                        coords_r[inds_dat2,:],
                                        dt_rl[:,:,inds_dat2])

        return {'FINISHED'}


class HEXBLENDER_OT_debug_hex_derivs(bpy.types.Operator):
    """Debug Hex Derivatives"""

    bl_idname = "hexblender.debug_hex_derivs"
    bl_label = "Debug Hex Derivatives"
    bl_description = ("Debug Hex Derivatives")
    bl_options = {'REGISTER'}

    @classmethod
    def poll(cls, context):
        return context.object is not None

    def execute(self, context):
        print("Not implemented yet")
        return {'FINISHED'}
        # I'll add this to the GUI
        print("\nDebug derivatives - ensure you specify an element number ")
        print("(or list of numbers of contiguous elements) indexed from zero ")
        print("in the original base mesh that need to be debugged in the variable hexesToDebug")
        
        # USER MUST SET ! Indexed from 0
        hexesToDebug = [11]
        #debugHexesOneNeigh = [1, 7, 8, 9, 13, 19, 20, 24]
        #debugHexesOneNeigh = [8, 10, 12, 12, 12, 12]
        #debugHexesOneNeigh = [7, 9, 11, 11, 11, 11]
        debugHexesOneNeigh = [6, 1, 3, 4, 6, 1]
        #debugHexesOneNeigh = []
        
        try:
            mesh, orig_mode = hexblender_helpers.get_active_mesh()
        except Exception as msg:
            print("Error getting mesh: %" % msg)
            return {'CANCELLED'}

        cubes, matList = hexblender_helpers.get_cached_data(mesh)

        baseHexes = sorted(hexesToDebug+debugHexesOneNeigh)
        
        hexesToDebugChildren = []
        adjHexesChildren = []
        for hexNum in hexesToDebug:
            hexesToDebugChildren.append(list(range(64*hexNum, 64*hexNum+64)))
        for hexNum in debugHexesOneNeigh:
            adjHexesChildren.append(list(range(64*hexNum, 64*hexNum+64)))
        #sum(x,[]) flattens a list of lists
        hexesToDebugChildren = list(set(sum(hexesToDebugChildren,[])))
        adjHexesChildren = list(set(sum(adjHexesChildren,[])))
        
        currHexes = np.concatenate((hexesToDebugChildren, adjHexesChildren))
        currHexes.sort()
        currDebugHexesSet = set(sorted(hexesToDebugChildren))
        hexMatNew = [1 if hexNum in currDebugHexesSet else 0 for hexNum in currHexes]
        cubes_group = np.array(cubes).copy()
        cubes_group = cubes_group[np.array(currHexes),:].tolist()
        
        # allVertInds has the ORIGINAL vertex numbers, before they are 
        # reduced from the mapping
        # vertSelectionMap tells you how to map the ORIGINAL vertex 
        # numbers to the NEW (REDUCED) vertex numbers
        allVertInds = np.unique(np.array(cubes_group[:]).flatten())
        allVertIndsSet = set(allVertInds)
        isVertIncluded = np.array([vert in allVertIndsSet for vert in range(len(mesh.verts))])
        vertSelectionMap = np.cumsum(isVertIncluded)-1
        
        # pluck out only the coordinates that we want and store 
        # them in verts_array
        verts_array = np.zeros([np.shape(allVertInds)[0],3])
        for ind,val in enumerate(verts_array):
            verts_array[ind,0] = mesh.verts[allVertInds[ind]].co.x
            verts_array[ind,1] = mesh.verts[allVertInds[ind]].co.y
            verts_array[ind,2] = mesh.verts[allVertInds[ind]].co.z
        
        # cubes of interest have been renumbered to have the NEW 
        # (REDUCED) vertex numbers
        cubes_reduced = vertSelectionMap[np.array(cubes_group)]
        cubes_trans = np.zeros([np.shape(cubes_reduced)[0],8])
        for ind,hex in enumerate(cubes_reduced):
            cubes_trans[ind,:] = [hex[0],hex[1],hex[3],hex[2],hex[4],hex[5],hex[7],hex[6]]
        
        hexblender_helpers.delete_all(mesh)
        newCubeFaces = hexblender_helpers.add_new_data(verts_array, [], [], cubes_trans)

        # I make the "element of interest" have material 1, 
        # "surrounding elements" material 0, because border faces end up 
        # belonging to material 1
        existing_mats = hexblender_helpers.get_mats_used_unused_sorted(mesh)
        hexblender_helpers.set_mat_data_hex(hexMatNew, 
                                            newCubeFaces,
                                            subdivided = False,
                                            existing_mats = existing_mats)
        
        # These are the suggested fixes
        print(baseHexes)
        print(cubes_trans)

        return {'FINISHED'}

class HEXBLENDER_OT_visualize_indices(bpy.types.Operator):
    bl_idname = "hexblender.visualize_indices"
    bl_label = "Index Visualizer"
    bl_description = "Toggle the visualization of indices"

    _handle = None

    @classmethod
    def poll(cls, context):
        return context.mode == "EDIT_MESH"

    def modal(self, context, event):
        if context.area:
            context.area.tag_redraw()

        try:
            hex_scene = context.scene.hexblender
        except Exception as msg:
            print("ERROR: Unable to get hexblender scene!  %s" % msg)
            return {'FINISHED'}

        # removal of callbacks when operator is called again
        if hex_scene.hexblender_properties.display_indices == -1:
            bpy.types.SpaceView3D.draw_handler_remove(self._handle, 'WINDOW')
            hex_scene.hexblender_properties.display_indices = 0
            return {"CANCELLED"}

        #hexblender_helpers.reset_mode("EDIT")
        return {"PASS_THROUGH"}

    def invoke(self, context, event):
        point_dict = {}
        try:
            hex_scene = context.scene.hexblender
        except Exception as msg:
            print("ERROR: Unable to get hexblender scene!  %s" % msg)
            return {'FINISHED'}

        mesh, orig_mode = hexblender_helpers.get_active_mesh()

        #if context.area.type == "VIEW_3D":
        if context.area.type == "PROPERTIES":
            if hex_scene.hexblender_properties.display_indices < 1:
                # operator is called for the first time, start everything
                hex_scene.hexblender_properties.display_indices = 1
                self._handle = bpy.types.SpaceView3D.draw_handler_add(
                    hexblender_helpers.draw_callback_px, (self, context, point_dict, mesh), 'WINDOW', 'POST_PIXEL')
                context.window_manager.modal_handler_add(self)
                hexblender_helpers.reset_mode(orig_mode)
                return {"RUNNING_MODAL"}
            else:
                # operator is called again, stop displaying
                hex_scene.hexblender_properties.display_indices = -1
                hexblender_helpers.reset_mode(orig_mode)
                return {'RUNNING_MODAL'}
        else:
            self.report({"WARNING"}, "View3D not found, can't run operator")
            hexblender_helpers.reset_mode(orig_mode)
            return {"CANCELLED"}



############## These are all just stubs for now!!!! ###################
class HEXBLENDER_OT_harmonize_topology(bpy.types.Operator):
    """Harmonize topology"""

    bl_idname = "hexblender.harmonize_topology"
    bl_label = "Harmonize Topology"
    bl_description = ("Harmonize element connectivity in selected topology regions`")
    bl_options = {'REGISTER'}

    @classmethod
    def poll(cls, context):
        return context.object is not None

    def execute(self, context):
        try:
            mesh, orig_mode = hexblender_helpers.get_active_mesh()
        except Exception as msg:
            print("Error getting mesh: %" % msg)
            return {'CANCELLED'}

        cubes, matList = hexblender_helpers.get_cached_data(mesh)

        #find out whether there are any non-intesecting regions that belong to the same labeled topology region
        #and separate them
        HexMatNew = contiguous_regions(list(ordered_hex_vertices), HexMat)
        
        print("Harmonizing topology regions...")
        
        #find verts and put them into an array
        verts_array = np.zeros([len(mesh.vertices),3])
        for ind, vert in enumerate(mesh.vertices):
            verts_array[ind,:] = vert.co
        
        #harmonize hexes
        ordered_hex_vertices = harmonize_topo(ordered_hex_vertices,verts_array,HexMatNew)

        print("Not Yet Implemented")
        return {'FINISHED'}


class HEXBLENDER_OT_linear_interp_subdivide(bpy.types.Operator):
    """Linear Sub"""

    bl_idname = "hexblender.linear_interp_subdivide"
    bl_label = "Linear Subdivide"
    bl_description = ("Perform linear subdivision")
    bl_options = {'REGISTER'}

    @classmethod
    def poll(cls, context):
        return context.object is not None

    def execute(self, context):
        print("Not Yet Implemented")
        return {'FINISHED'}


class HEXBLENDER_OT_fix_hex_keys(bpy.types.Operator):
    """Fix Hex Keys"""

    bl_idname = "hexblender.fix_hex_keys"
    bl_label = "Fix Hex Keys"
    bl_description = ("Fix Hex Keys")
    bl_options = {'REGISTER'}

    @classmethod
    def poll(cls, context):
        return context.object is not None

    def execute(self, context):
        print("Not Yet Implemented")
        return {'FINISHED'}


class HEXBLENDER_OT_recalc_normals(bpy.types.Operator):
    """Recalc/Flip Normals"""

    bl_idname = "hexblender.recalc_normals"
    bl_label = "Recalc/Flip Normals"
    bl_description = ("Recalc/Flip Normals")
    bl_options = {'REGISTER'}

    @classmethod
    def poll(cls, context):
        return context.object is not None

    def execute(self, context):
        print("Not Yet Implemented")
        return {'FINISHED'}

class HEXBLENDER_OT_show_fverts(bpy.types.Operator):
    """Show FVerts"""

    bl_idname = "hexblender.show_fverts"
    bl_label = "Show FVerts"
    bl_description = ("Show FVerts")
    bl_options = {'REGISTER'}

    @classmethod
    def poll(cls, context):
        return context.object is not None

    def execute(self, context):
        print("Not Yet Implemented")
        return {'FINISHED'}

class HEXBLENDER_OT_interp_subdivide(bpy.types.Operator):
    """Interp Subdivide (1 to 4)"""

    bl_idname = "hexblender.interp_subdivide"
    bl_label = "Interp Subdivide (1 to 4)"
    bl_description = ("Interp Subdivide (1 to 4)")
    bl_options = {'REGISTER'}

    @classmethod
    def poll(cls, context):
        return context.object is not None

    def execute(self, context):
        print("Not Yet Implemented")
        return {'FINISHED'}

class HEXBLENDER_OT_linear_subdivide(bpy.types.Operator):
    """Linear Subdivide (1 to 9)"""

    bl_idname = "hexblender.linear_subdivide"
    bl_label = "Linear Subdivide (1 to 9)"
    bl_description = ("Linear Subdivide (1 to 9)")
    bl_options = {'REGISTER'}

    @classmethod
    def poll(cls, context):
        return context.object is not None

    def execute(self, context):
        print("Not Yet Implemented")
        return {'FINISHED'}

class HEXBLENDER_OT_extrude(bpy.types.Operator):
    """Extrude"""

    bl_idname = "hexblender.extrude"
    bl_label = "Extrude"
    bl_description = ("Extrude")
    bl_options = {'REGISTER'}

    @classmethod
    def poll(cls, context):
        return context.object is not None

    def execute(self, context):
        print("Not Yet Implemented")
        return {'FINISHED'}

class HEXBLENDER_OT_change_offset(bpy.types.Operator):
    """Extrude"""

    bl_idname = "hexblender.change_offset"
    bl_label = "Change Offset"
    bl_description = ("Change Offset")
    bl_options = {'REGISTER'}

    @classmethod
    def poll(cls, context):
        return context.object is not None

    def execute(self, context):
        print("Not Yet Implemented")
        return {'FINISHED'}


class HEXBLENDER_OT_export_global_to_local_map(bpy.types.Operator):
    """Export Global to Local Map"""

    bl_idname = "hexblender.export_global_to_local_map"
    bl_label = "Export Global to Local Map"
    bl_description = ("Export Global to Local Map")
    bl_options = {'REGISTER'}

    @classmethod
    def poll(cls, context):
        return context.object is not None

    def execute(self, context):
        print("Not Yet Implemented")
        return {'FINISHED'}


class HEXBLENDER_OT_export_cont_fitting(bpy.types.Operator):
    """Export Cont Fitting"""

    bl_idname = "hexblender.export_cont_fitting"
    bl_label = "Export Cont Fitting"
    bl_description = ("Export Cont Fitting")
    bl_options = {'REGISTER'}

    @classmethod
    def poll(cls, context):
        return context.object is not None

    def execute(self, context):
        print("Not Yet Implemented")
        return {'FINISHED'}       


class HEXBLENDER_OT_get_vsoln_input_file(bpy.types.Operator):
    '''
    Sets input filelocation
    '''
    bl_idname = "hexblender.get_vsoln_file"
    bl_label = "Select voltage solution file (*.npy)"
    bl_description = ("Get voltage solution from file")
    bl_options = {'REGISTER'}

    filepath = bpy.props.StringProperty(subtype='FILE_PATH', default="")

    def execute(self, context):
        try:
            hex_scene = context.scene.hexblender
        except Exception as msg:
            print("ERROR: Unable to get hexblender scene!  %s" % msg)
            return {'FINISHED'}

        hex_scene.hexblender_properties.render_electrophysiology_vsoln_input_file = self.filepath
        return {'FINISHED'}

    def invoke(self, context, event):
        context.window_manager.fileselect_add(self)
        return {'RUNNING_MODAL'}


class HEXBLENDER_OT_render_electrophysiology_load(bpy.types.Operator):
    """Load electrophysiology solution from Vsoln"""

    bl_idname = "hexblender.render_electrophysiology_load"
    bl_label = "Load electrophysiology solution from Vsoln"
    bl_description = ("Load electrophysiology solution from Vsoln")
    bl_options = {'REGISTER'}

    @classmethod
    def poll(cls, context):
        return context.object is not None

    def execute(self, context):

        # function to set vertex colors based on initial frame
        def init_vcols(frame):
            m
            color_layer = mesh.vertex_colors[0]
            # create array to store direct mapping from vsoln_colormap to color layer
            # to avoid the two 'for' loops during subsequent updates
            reducedMap = np.zeros(len(color_layer.data),dtype='int')
            j=0
            for poly in mesh.polygons: # loop over all polygon faces in mesh
                for idx in poly.vertices: # loop over vertices of each polygon
                    color_layer.data[j].color = vsoln_colormap[frame,idx][0:3]
                    reducedMap[j] = idx
                    j = j+1
            # print(vsoln_colormap.shape)
            # print(vsoln_colormap[0,3][0:3])
            bpy.ops.object.mode_set(mode='VERTEX_PAINT')
            # print(reducedMap)
            return reducedMap

        # function to set vertex colors based on current frame
        def set_vcols(frame):
            obj = bpy.context.active_object
            mesh = obj.data
            color_layer = mesh.vertex_colors[0]
            j=0
            # use reducedMap to update color in one 'for' loop
            for polyvert in color_layer.data:
                polyvert.color = vsoln_colormap[frame,reducedMap[j]][0:3]
                j+=1
            # print(vsoln_colormap.shape)
            # print(vsoln_colormap[0,3][0:3])
            bpy.ops.object.mode_set(mode='VERTEX_PAINT')
            return {'FINISHED'}

        # handle to update vertex colors in viewport when frame is changed
        def handle_update_vsoln_colormap(scene):
            frame = scene.frame_current
            color_layer = mesh.vertex_colors[0]
            set_vcols(frame)
            bpy.ops.object.mode_set(mode='VERTEX_PAINT')
            return {'FINISHED'}

        try:
            hex_scene = context.scene.hexblender
        except Exception as msg:
            print("ERROR: Unable to get hexblender scene!  %s" % msg)
            return {'FINISHED'}

        mesh, orig_mode = hexblender_helpers.get_active_mesh()
        # get path to input voltage solution file (already converted to colormap)
        vsoln_input_file = hex_scene.hexblender_properties.render_electrophysiology_vsoln_input_file

        if vsoln_input_file == None:
            return{'FINISHED'}

        print("Reading voltage solution colormap data...")
        vsoln_colormap = np.load(vsoln_input_file)
        print("Got it!")
        #print(vsoln_colormap.shape)
        #print(vsoln_colormap[0,3][0:3])

        # create a new vertex color layer if one does not already exist
        if not mesh.vertex_colors:
            mesh.vertex_colors.new()

        # set vertex colors to voltage solution of first frame (t=0 ms)
        # return direct mapping from voltage solution to vertex color
        reducedMap = init_vcols(0)

        # set number of frames to length of voltage solution
        scn = bpy.context.scene
        scn.frame_start = 0
        scn.frame_end = vsoln_colormap.shape[0]
        # move current frame to first frame
        scn.frame_set(frame =0)

        # register the handle to run whenever frame is changed in the GUI
        bpy.app.handlers.frame_change_pre.append(handle_update_vsoln_colormap)
        return {'FINISHED'}