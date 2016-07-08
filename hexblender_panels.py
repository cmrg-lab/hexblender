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
This script draws the panels and other UI elements for HexBlender.

"""

# blender imports
import bpy
#from bpy.types import Menu


# python imports
import re
import os

# hexblender imports
import hexblender

# we use per module class registration/unregistration
def register():
    bpy.utils.register_module(__name__)


def unregister():
    bpy.utils.unregister_module(__name__)


#HexBlendereGUI Panels:
class HEXBLENDER_PT_exports(bpy.types.Panel):
    bl_label = "Hex - Exports"
    bl_space_type = "VIEW_3D"
    bl_region_type = "TOOLS"
    bl_category = "Hex Render/Export"
    bl_options = {'DEFAULT_CLOSED'}

    def draw(self, context):
        layout = self.layout
        hex_scene = context.scene.hexblender

        row = layout.row(align=True)
        row.operator("hexblender.export_vertices", text="Cont. Vertices")
        row.prop(hex_scene.hexblender_properties, "export_vert_weights")

        row = layout.row(align=True)
        if int(hex_scene.hexblender_properties.export_vert_weights) == 1:
            row.prop(hex_scene.hexblender_properties, "group_name")
            row = layout.row(align=True)
            col = row.column(align=True)
            # CV 02/29/16: Need to change behavior of export_vertices to use export_vert_weights; for now, duplicate button
            #col.prop(hex_scene.hexblender_properties, "group_weight_add", icon='ZOOMIN', text="")
            #col.prop(hex_scene.hexblender_properties, "group_weight_remove", icon='ZOOMOUT', text="")
            row = layout.row(align=True)
        ###############################################

        row = layout.row(align=False)
        row.operator("hexblender.export_elements", text="Hex Elements")
        row.operator("hexblender.export_face_keys", text="Face Elements")

        #### MD 1/19/2016 - VERTEX WEIGHTS - SCALED ####
        layout.separator()
        row = layout.row(align=True)
        row.operator("hexblender.export_vertex_weights", text="Vertex Weights")
        row.prop(hex_scene.hexblender_properties, "vert_weight_scalar")

        row = layout.row(align=True)
        ####### PROOF OF CONCEPT - FUNCTIONALITY ########
        # TRICUBIC DERIVS
        layout.separator()
        row = layout.row(align=True)
        row.operator("hexblender.export_hermite_tricubic_derivs", text="Hermite Tricubic Derivs")
        row.prop(hex_scene.hexblender_properties, "tri_cubic_iters")
        row = layout.row(align=True)

        row.prop(hex_scene.hexblender_properties, "export_subdiv_types")
        row = layout.row(align=True)
        if int(hex_scene.hexblender_properties.export_subdiv_types) > 0:
            row.prop(hex_scene.hexblender_properties, "delete_orig_mesh")
            row = layout.row(align=True)

        if int(hex_scene.hexblender_properties.export_subdiv_types) == 1:
            row.prop(hex_scene.hexblender_properties, "interp_priorities")
            row = layout.row(align=True)
            row.prop(hex_scene.hexblender_properties, "spline_priorities")
            row = layout.row(align=True)
            row.prop(hex_scene.hexblender_properties, "reg_priorities")

        layout.separator()
        row = layout.row(align=True)
        row.operator("hexblender.export_bundled_3d_mesh", text="Export Bundled Tricubic Mesh")
        row.prop(hex_scene.hexblender_properties, "tri_cubic_iters")
        row = layout.row(align=True)
        row.prop(hex_scene.hexblender_properties, "interp_priorities")
        row = layout.row(align=True)
        row.prop(hex_scene.hexblender_properties, "spline_priorities")
        row = layout.row(align=True)
        row.prop(hex_scene.hexblender_properties, "reg_priorities")

        row = layout.row(align=False)
        row.operator("hexblender.export_hermite_bicubic_derivs", text="Hermite Bicubic Derivs")
        row.operator("hexblender.export_bundled_2d_mesh", text="Bundled Bicubic Mesh")

        row = layout.row(align=False)
        row.operator("hexblender.export_global_to_local_map", text="Global to Local Map")
        row.operator("hexblender.export_fit_data_points", text="Fit Data Points")

        #row = layout.row(align=True)
        #row.operator("hexblender.export_cont_fitting", text="Export Cont Fitting")


class HEXBLENDER_PT_imports(bpy.types.Panel):
    bl_label = "Hex - Imports"
    bl_space_type = "VIEW_3D"
    bl_region_type = "TOOLS"
    bl_category = "HexBlender Tools"

    def draw(self, context):
        layout = self.layout
        #scn = context.scene

        row = layout.row(align=True)
        row.operator("hexblender.import_pickle", text="Import Pickle File")


class HEXBLENDER_PT_hex_tools(bpy.types.Panel):
    bl_label = "Hex - Tools"
    bl_space_type = "VIEW_3D"
    bl_region_type = "TOOLS"
    bl_category = "HexBlender Tools"
    bl_options = {'DEFAULT_CLOSED'}

    def draw(self, context):
        layout = self.layout
        hex_scene = context.scene.hexblender

        row = layout.row(align=True)
        row.operator("hexblender.update_pickle_hexs", text="Find/update hexs")
        row = layout.row(align=True)
        row.operator("hexblender.rotate_hexs", text="Rotate Hexs")
        row.prop(hex_scene.hexblender_properties, "rotate_matl_index")
        row.prop(hex_scene.hexblender_properties, "rotate_types")

        layout.separator()


        row = layout.row(align=True)
        row.operator("hexblender.hex_interp_subdivide", text="Interpolatory subdivision (1-8)")
        row = layout.row(align=True)
        # layout.separator()
        row.prop(hex_scene.hexblender_properties, "subdiv_types")
        row = layout.row(align=True)
        row.prop(hex_scene.hexblender_properties, "thin_plate_mapping")
        if int(hex_scene.hexblender_properties.subdiv_types) == 0:
            row = layout.row(align=True)
            row.prop(hex_scene.hexblender_properties, "interp_priorities")
            row = layout.row(align=True)
            if hex_scene.hexblender_properties.thin_plate_mapping:
                row = layout.row(align=True)
                row.prop(hex_scene.hexblender_properties, "spline_priorities")
        row = layout.row(align=True)
        row.prop(hex_scene.hexblender_properties, "delete_orig_mesh")
        #row.operator("hexblender.linear_interp_subdivide", text="Linear 1-to-8 subdivision")
        layout.separator()

        row = layout.row(align=True)
        row.operator("hexblender.regularize_hexs", text="Regularize Hexs")
        row.prop(hex_scene.hexblender_properties, "regularize_hex_iters")
        row = layout.row(align=True)
        row.prop(hex_scene.hexblender_properties, "preserve_ridges")
        row = layout.row(align=True)
        row.prop(hex_scene.hexblender_properties, "immobilize_extraord")
        row = layout.row(align=True)
        row.prop(hex_scene.hexblender_properties, "tangent_motion")
        row = layout.row(align=True)
        row.prop(hex_scene.hexblender_properties, "normal_motion")
        row = layout.row(align=True)
        row.prop(hex_scene.hexblender_properties, "internal_motion")
        # row = layout.row(align=True)        
        # row.prop(hex_scene.hexblender_properties, "recacheHexs")

        # row = layout.row(align=True)
        # row.operator("hexblender.harmonize_topology", text="Harmonize topology region")

        # row = layout.row(align=True)
        # row.operator("hexblender.debug_hex_derivs", text="Debug Hex Derivs")

        # row = layout.row(align=True)
        # row.operator("hexblender.fix_hex_keys", text="Fix Hex Keys")

class HEXBLENDER_PT_sfc_tools(bpy.types.Panel):
    bl_label = "Hex- Quad Tools"
    bl_space_type = "VIEW_3D"
    bl_region_type = "TOOLS"
    bl_category = "HexBlender Tools"
    bl_options = {'DEFAULT_CLOSED'}

    def draw(self, context):
        layout = self.layout
        hex_scene = context.scene.hexblender

        row = layout.row(align=True)
        row.operator("hexblender.regularize_sfc", text="Regularize Surface")
        row.prop(hex_scene.hexblender_properties, "regularize_sfc_iters")

        row = layout.row(align=True)
        row.operator("hexblender.recalc_normals", text="Recalc/Flip Normals")

        row = layout.row(align=True)
        row.operator("hexblender.show_fverts", text="Show FVerts")

        row = layout.row(align=True)
        row.operator("hexblender.interp_subdivide", text="Interp Sub (1 to 4)")

        row = layout.row(align=True)
        row.operator("hexblender.linear_subdivide", text="Linear Subdivide (1 to 9)")
        row.prop(hex_scene.hexblender_properties, "linear_sub_iters")

        row = layout.row(align=True)
        row.operator("hexblender.extrude", text="Extrude")
        row.prop(hex_scene.hexblender_properties, "extrude_offset")

        row = layout.row(align=True)
        row.operator("hexblender.change_offset", text="Change Offset")
        row.prop(hex_scene.hexblender_properties, "extrude_offset")


class HEXBLENDER_PT_selection_tools(bpy.types.Panel):
    bl_label = "Hex - Selection Tools"
    bl_space_type = "VIEW_3D"
    bl_region_type = "TOOLS"
    bl_category = "HexBlender Tools"
    bl_options = {'DEFAULT_CLOSED'}

    def draw(self, context):
        layout = self.layout
        hex_scene = context.scene.hexblender

        row = layout.row(align=True)
        row.operator("hexblender.adjust_nodes_circ", text="Adjust vertex - circular")
        row.prop(hex_scene.hexblender_properties, "radius_of_nodal_adjustment")
        row = layout.row(align=True)
        row.operator("hexblender.adjust_nodes_linear", text="Adjust vertex - linear")
        layout.separator()

        row = layout.row(align=True)
        row.operator("hexblender.edge_lengths", text="Std Dev of Edge Lengths")
        row = layout.row(align=True)
        row.operator("hexblender.std_faces", text="Std Dev of Faces")
        row = layout.row(align=True)
        row.operator("hexblender.std_angles", text="Std Dev of Angles")
        row = layout.row(align=True)
        row.operator("hexblender.print_data", text="Print Selection Data")

        layout.separator()

        row = layout.row(align=True)
        row.prop(hex_scene.hexblender_properties, "vert_index_number")
        row = layout.row(align=True)
        row.operator("hexblender.select_vertex", text="Select Vert(s)")
        layout.separator()

        row = layout.row(align=True)
        row.prop(hex_scene.hexblender_properties, "edge_index_number")
        row = layout.row(align=True)
        row.operator("hexblender.select_edge", text="Select Edge(s)")
        layout.separator()

        row = layout.row(align=True)
        row.prop(hex_scene.hexblender_properties, "face_index_number")
        row = layout.row(align=True)
        row.operator("hexblender.select_face", text="Select Face(s)")
        layout.separator()

        row = layout.row(align=True)
        row.prop(hex_scene.hexblender_properties, "hex_index_number")
        row = layout.row(align=True)
        row.operator("hexblender.select_hex", text="Select Hex(es)")
        layout.separator()

        col = layout.column(align=True)
        col.operator("hexblender.visualize_indices", text="Visualize indices")
        row = col.row(align=True)
        row.active = (context.mode == "EDIT_MESH" and 
            hex_scene.hexblender_properties.display_indices == 1)

        row.prop(hex_scene.hexblender_properties, "display_vert_index", toggle=True)
        row.prop(hex_scene.hexblender_properties, "display_edge_index", toggle=True)
        row.prop(hex_scene.hexblender_properties, "display_face_index", toggle=True)
        row.prop(hex_scene.hexblender_properties, "display_hex_index", toggle=True)
        row = col.row(align=True)
        row.active = context.mode == "EDIT_MESH" and hex_scene.hexblender_properties.display_indices == 1
        row.prop(hex_scene.hexblender_properties, "display_sel_only")


class HEXBLENDER_PT_tensor_tools(bpy.types.Panel):
    bl_label = "Hex - Tensor Tools"
    bl_space_type = "VIEW_3D"
    bl_region_type = "TOOLS"
    bl_category = "HexBlender Tools"
    bl_options = {'DEFAULT_CLOSED'}

    def draw(self, context):
        layout = self.layout
        hex_scene = context.scene.hexblender

        row = layout.row(align=True)
        row.prop(hex_scene.hexblender_properties, "unaligned_object")
        row = layout.row(align=True)
        row.prop(hex_scene.hexblender_properties, "aligned_object")

        #row.operator("hexblender.tensor_apply_xform", text="Apply the transformation to the selected data file")
        row = layout.row(align=True)
        row.operator("hexblender.get_tensor_input_file", text="Set tensor input file")
        row = layout.row(align=True)
        row.prop(hex_scene.hexblender_properties, "tensor_input_file")

        row = layout.row(align=True)
        row.operator("hexblender.get_tensor_output_file", text="Set tensor output file")
        row = layout.row(align=True)
        row.prop(hex_scene.hexblender_properties, "tensor_output_file")

        row = layout.row(align=True)
        row.operator("hexblender.compute_tensor_xform", 
            text="Compute the transformation")


#HexBlendereGUI Panels:
class HEXBLENDER_PT_render(bpy.types.Panel):
    bl_label = "Hex - Render"
    bl_space_type = "VIEW_3D"
    bl_region_type = "TOOLS"
    bl_category = "Hex Render/Export"
    bl_options = {'DEFAULT_CLOSED'}

    def draw(self, context):
        layout = self.layout
        hex_scene = context.scene.hexblender

        row = layout.row(align=True)
        row.operator("hexblender.get_vsoln_file", text="Select voltage solution file")
        # row = layout.row(align=True)
        row.prop(hex_scene.hexblender_properties, "render_electrophysiology_vsoln_input_file")
        row = layout.row(align=True)
        row.operator("hexblender.render_electrophysiology_load", text="Load voltage solution")
        row = layout.row(align=True)
        row.prop(hex_scene.hexblender_properties, "update_vsoln_with_frame")

        layout.separator()

        row = layout.row(align=True)
        row.prop(hex_scene.hexblender_properties, "render_biomechanics")
