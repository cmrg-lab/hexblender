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
This script contains the custom properties used in HexBlender.

"""
# blender imports
import bpy
from hexblender.hexblender_helpers import get_material_data
from bpy.props import BoolProperty, CollectionProperty, EnumProperty, \
    FloatProperty, FloatVectorProperty, IntProperty, IntVectorProperty,\
    PointerProperty, StringProperty


# we use per module class registration/unregistration
def register():
    bpy.utils.register_module(__name__)


def unregister():
    bpy.utils.unregister_module(__name__)

def get_mat_list():
    mat_list = []

    if bpy.context.active_object is not None and\
       bpy.context.active_object.type == "MESH":

        for index in range(len(bpy.context.active_object.data.materials)):
            mat_list.append(index)
    else:
        mat_list = ''

    return mat_list

def get_interp_priorities(self):
    return self.get('interp_priorities', '%s' % get_mat_list())

def set_interp_priorities(self, value):
    self["interp_priorities"] = value

def get_spline_priorities(self):
    return self.get('spline_priorities', '%s' % get_mat_list())

def set_spline_priorities(self, value):
    self["spline_priorities"] = value

def get_reg_priorities(self):
    return self.get('reg_priorities', '%s' % get_mat_list())

def set_reg_priorities(self, value):
    self["reg_priorities"] = value


class HexBlenderObjectPanelProperty(bpy.types.PropertyGroup):
    tri_cubic_iters = IntProperty(
        name="Iterations", default=0, min=0, 
        description="Number of iterations to perform")

    interp_priorities = StringProperty(
        name="Interp Priorities", 
        description="Defines priority groupings (by topology/material region(s)) for simple interpolation. Ex: [0,1] or [0][1]", get=get_interp_priorities, set=set_interp_priorities)

    spline_priorities = StringProperty(
        name="Spline Priorities", 
        description="Defines priority groupings (by topology/material region(s)) for spline interpolation. Ex: [0,1] or [0],[1]",
        get=get_spline_priorities, set=set_spline_priorities)

    reg_priorities = StringProperty(
        name="Reg Priorities", 
        description="Defines the regularization priorities of the material surfaces. Ex: [0,1] or [0],[1]",
        get=get_reg_priorities, set=set_reg_priorities)

    tensor_input_file = StringProperty(
        name="Tensor input file", 
        description="Defines the file that contains tensor data that will be transformed.")

    tensor_output_file = StringProperty(
        name="Tensor output file", 
        description="Defines the file that will contain the transformed tensor data.")

    unaligned_object = StringProperty(
        name="The unaligned object", 
        description="Defines the unaligned coordinates for the transformation.")

    aligned_object = StringProperty(
        name="The aligned object", 
        description="Defines the aligned coordinates for the transformation.")

    regularize_hex_iters = IntProperty(
        name="Iterations", default=5, min=1, 
        description="Number of iterations to perform")

    regularize_sfc_iters = IntProperty(
        name="Iterations", default=10, min=1, 
        description="Number of iterations to perform")

    linear_sub_iters = IntProperty(
        name="Iterations", default=10, min=1, 
        description="Number of iterations to perform")

    material_num = IntProperty(
        name="Material Num", default=1, min=1, 
        description="Material Property Number")

    offset_ratio = IntProperty(
        name="Ratio (%)", default=10, min=1, 
        description="Offset Ratio")

    vert_index_number = StringProperty(
        name="Vetex Num(s)",
        description="Vertices wanted, one-based. Ex: 1-5,8,10")

    edge_index_number = StringProperty(
        name="Edge Num(s)",
        description="Edges wanted, one-based. Ex: 1-5,8,10")

    face_index_number = StringProperty(
        name="Face Num(s)",
        description="Faces wanted, one-based. Ex: 1-5,8,10")

    hex_index_number = StringProperty(
        name="Hex Num(s)",
        description="Hexes wanted, one-based. Ex: 1-5,8,10")

    radius_of_nodal_adjustment = FloatProperty(
        name="Radius", default=0.10, min=0.0, 
        description="Radius of nodal adjustment")

    extrude_offset = FloatProperty(
        name="Radius", default=0.10, min=0.0, 
        description="Extrude Offset")

    rotate_matl_index = IntProperty(
        name="Material", default=0, min=0, 
        description="The index (from 0) of the topology/material region to rotate")

    delete_orig_mesh = BoolProperty(
        name="Delete the original mesh", default=False,
        description="Delete the original mesh after subdividing")

    display_indices = IntProperty(
        name="Display indices",
        default=0)

    display_sel_only = BoolProperty(
        name="Selected only",
        description="Only display indices of selected vertices/edges/faces",
        default=True)

    display_vert_index = BoolProperty(
        name="Vertices",
        description="Display vertex indices", default=True)

    display_edge_index = BoolProperty(
        name="Edges",
        description="Display edge indices")

    display_face_index = BoolProperty(
        name="Faces",
        description="Display face indices")

    display_hex_index = BoolProperty(
        name="Hexes",
        description="Display hex indices")

    live_mode = BoolProperty(
        name="Live",
        description="Toggle live update of the selection, can be slow",
        default=False)

    subdiv_types_enum = [("0","Subdivide by priority region","Subdivide region-wise according to priority",0),
                         ("1","Subdivide by element (Experimental)","Subdivide element-wise",1)]

    subdiv_types = EnumProperty(name = "Subdivision method", 
                                description = "Select the method",
                                default = "0", 
                                items = subdiv_types_enum)

    export_subdiv_types_enum = [("0","Use existing subdivided mesh","Compute derivatives using existing mesh",0),
                         ("1","Subdivide by priority region","Subdivide region-wise according to priority and compute derivatives",1),
                         ("2","Subdivide by element (Experimental)","Subdivide element-wise, and compute derivatives",2)]

    export_subdiv_types = EnumProperty(name = "Subdivision method", 
                                description = "Select the method",
                                default = "0", 
                                items = export_subdiv_types_enum)

    thin_plate_mapping = BoolProperty(
        name="Use thin-plate spline mapping",
        description="Use thin-plate splines to compute interpolated subdivision surface (recommended for coarse linear meshes)",
        default=True)

    preserve_ridges = BoolProperty(
        name="Preserve ridges", default=True,
        description="")

    immobilize_extraord = BoolProperty(
        name="Immobilize extraordinary vertices", default=False,
        description="Disable tangent motion of extraordinary vertices")

    tangent_motion = BoolProperty(
        name="Tangent motion", default=True,
        description="Enable or disable tangent motion of regularized vertices")

    normal_motion = BoolProperty(
        name="Normal motion", default=True,
        description="Enable or disable normal motion of regularized vertices")

    internal_motion = BoolProperty(
        name="Internal motion", default=True,
        description="Enable or disable internal motion of regularized vertices")

    recacheHexs = BoolProperty(
        name="Update hex cache", default=True,
        description="Update hex cache with newly found hexs")

    render_electrophysiology_vsoln_input_file = StringProperty(
        name="", 
        description="Path to electrophysiology voltage solution (Vsoln*.npy).")

    update_vsoln_with_frame = BoolProperty(
        name="Update voltage solution with current timeline frame", 
        description="Check to view voltage solution as vertex paint layer updated by current timeline frame")

    render_biomechanics = StringProperty(
        name="Biomechanics", 
        description="Render biomechanics solution.")

# These are the original descriptions, keeping for the time being for
# reference
#  1: keep origin at nodal position 1, 
#     turn xi3 into xi1, xi1 into xi2, xi2 into xi3 
#  2: move origin to nodal position 2, 
#     turn xi2 into xi1, turn xi1 into xi2 and flip, preserve xi3
#  3: move origin to nodal position 2, 
#     turn xi3 into xi1, preserve xi2, turn xi1 and xi3 and flip
#  4: move origin to nodal position 3, 
#     turn xi2 into xi1 and flip, turn xi1 into xi2, preserve xi3
#  5: move origin to nodal position 3, 
#     turn xi3 into xi1,flip xi2, turn xi2 into xi3 
#  6: move origin to nodal position 4, 
#     flip xi1, flip xi2, preserve xi3
#  7: move origin to nodal position 4, 
#     turn xi3 into xi1, turn xi1 into xi2 and flip, turn xi2 into xi3 and flip
#  8: move origin to nodal position 5, 
#     turn xi2 into xi1,turn xi1 into xi2, flip xi3
#  9: move origin to nodal position 5, 
#     turn xi3 into xi1 and flip, preserve xi2, turn xi1 into xi3
#  10: move origin to nodal position 6, 
#      flip xi1, keep xi2, flip xi3
#  11: move origin to nodal position 6, 
#      turn xi3 into xi1 and flip, turn xi1 into xi2 and flip, turn xi2 into xi3
#  12: move origin to nodal position 7, 
#      preserve xi1, flip xi2, flip xi3
#  13: move origin to nodal position 7, 
#      turn xi3 into xi1 and flip, turn xi1 into xi2, turn xi2 into xi3 and flip
#  14: move origin to nodal position 8, 
#      turn xi2 into xi1 and flip, turn xi1 into xi2 and flip, flip xi3
#  15: move origin to nodal position 8, 
#      turn xi3 into xi1 and flip,flip xi2, turn xi1 into xi3 and flip

    rot_types_enum = [("0","Do Nothing","Do Nothing",0),
                      ("1",
                       "Origin at position 1...",
                       "Origin unchanged, xi1 -> xi2, xi2 -> xi3, xi3 -> xi1",1),
                      ("2",
                       "Move origin to pos 2, same xi3...",
                       "Move origin to nodal pos 2, xi1 -> xi2 and flip, xi2 -> xi1, preserve xi3",2),
                      ("3",
                       "Move origin to pos 2, same xi2...",
                       "Move origin to nodal pos 2, xi1 -> xi3 and flip, xi3 -> xi1, preserve xi2",3),
                      ("4",
                       "Move origin to pos 3, same xi3...",
                       "Move origin to nodal pos 3, xi1 -> xi2, xi2 -> xi1 and flip, preserve xi3",4),
                      ("5",
                       "Move origin to pos 3...",
                       "Move origin to nodal pos 3, xi3 -> xi1, flip xi2, xi2 -> xi3",5),
                      ("6",
                       "Move origin to pos 4, same xi3...",
                       "Move origin to nodal pos 4, flip xi1, flip xi2, preserve xi3",6),
                      ("7",
                       "Move origin to pos 4...",
                       "Move origin to nodal pos 3, xi1 -> xi2 and flip, xi2 -> xi3 and flip, xi3 -> xi1",7),
                      ("8",
                       "Move origin to pos 5...",
                       "Move origin to nodal pos 5, xi1 -> xi2, xi2 -> xi1, flip xi3",8),
                      ("9",
                       "Move origin to pos 5, same xi2...",
                       "Move origin to nodal pos 5, xi1 -> xi3, xi3 -> xi1 and flip, preserve xi2",9),
                      ("10",
                       "Move origin to pos 6, same xi2...",
                       "Move origin to nodal pos 6, flip xi1, preserve xi2, flip xi3",10),
                      ("11",
                       "Move origin to pos 6, same xi2...",
                       "Move origin to nodal pos 6, xi1 -> xi2 and flip, xi2 -> xi3, xi3 -> xi1 and flip",11),
                      ("12",
                       "Move origin to pos 7, same xi1...",
                       "Move origin to nodal pos 7, preserve xi1, flip xi2, flip xi3",12),
                      ("13",
                       "Move origin to pos 7...",
                       "Move origin to nodal pos 7, xi1 -> xi2, xi2 -> xi3 and flip, xi3 -> xi1 and flip",13),
                      ("14",
                       "Move origin to pos 8, flip xi3...",
                       "Move origin to nodal pos 8, xi1 -> xi2 and flip, xi2 -> xi1 and flip, flip xi3",14),
                      ("15",
                       "Move origin to pos 8, flip xi2...",
                       "Move origin to nodal pos 8, xi1 -> xi3 and flip, flip xi2, xi3 -> xi1",15)]

    rotate_types = EnumProperty(name = "Type", 
                                description = "Select the type of rotation",
                                default = "0", 
                                items = rot_types_enum)



# Main HexBlender Properties Class:
class HexBlenderPropertyGroup(bpy.types.PropertyGroup):
    hexblender_properties = PointerProperty(
        type=HexBlenderObjectPanelProperty,
        name="HexBlender Properties")
