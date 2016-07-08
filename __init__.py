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

import os

bl_info = {
    "name": "HexBlender",
    "author": "Cardiac Mechanics Research Group, UCSD Bioengineering Dept.",
    "version": (1,7),
    "blender": (2, 70, 5),
    "location": "Properties > Object > HexBlender Panel and Properties -> Scene > HexBlender Panel",
    "description": "HexBlender Mesh Tools for Continuity",
    "warning": "",
    "wiki_url": "https://github.com/cmrglab/hexblender",
    "category": "Mesh Tools"
}

hexblender_info = {
    "supported_version_list": [(2, 70, 5)],
    "hexblender_source_list": [ "__init__.py",
                                "hexblender_panels.py",
                                "hexblender_operators.py",
                                "hexblender_helpers.py",
                                "hex_interp_subdiv.py"],
    "hexblender_addon_path": "",
}


#See if Removing this still works!!!
'''
if __name__ == '__main__':
    print("Came1")
    print("HexBlender is running as __main__")
    identify_source_version("")
'''

# To support reload properly, try to access a package var.
# If it's there, reload everything
if "bpy" in locals():
    print("Reloading HexBlender")
    import imp
    imp.reload(hexblender_panels)
    imp.reload(hexblender_operators)
    imp.reload(hexblender_helpers)
    imp.reload(hexblender_properties)
else:
    print("Importing HexBlender")
    from . import hexblender_panels
    from . import hexblender_operators
    from . import hexblender_helpers
    from . import hexblender_properties

import bpy
import sys


# We use per module class registration/unregistration
def register():

    bpy.utils.register_module(__name__)

    # Unregister and re-register panels to display them in order
    bpy.utils.unregister_class(hexblender_panels.HEXBLENDER_PT_imports)
    bpy.utils.unregister_class(hexblender_panels.HEXBLENDER_PT_exports)
    bpy.utils.unregister_class(hexblender_panels.HEXBLENDER_PT_hex_tools)
    bpy.utils.unregister_class(hexblender_panels.HEXBLENDER_PT_sfc_tools)
    bpy.utils.unregister_class(hexblender_panels.HEXBLENDER_PT_selection_tools)
    bpy.utils.unregister_class(hexblender_panels.HEXBLENDER_PT_tensor_tools)
    bpy.utils.unregister_class(hexblender_panels.HEXBLENDER_PT_render)


    bpy.utils.register_class(hexblender_panels.HEXBLENDER_PT_imports)
    bpy.utils.register_class(hexblender_panels.HEXBLENDER_PT_exports)

    bpy.utils.register_class(hexblender_panels.HEXBLENDER_PT_hex_tools)
    bpy.utils.register_class(hexblender_panels.HEXBLENDER_PT_sfc_tools)
    bpy.utils.register_class(hexblender_panels.HEXBLENDER_PT_selection_tools)
    bpy.utils.register_class(hexblender_panels.HEXBLENDER_PT_tensor_tools)

    bpy.utils.register_class(hexblender_panels.HEXBLENDER_PT_render)

    # Properties
    bpy.types.Scene.hexblender = bpy.props.PointerProperty(
        type=hexblender_properties.HexBlenderPropertyGroup)
    #bpy.types.Object.hexblender = bpy.props.PointerProperty(
    #    type=hexblender_properties.HexBlenderPropertyGroup)

    print("HexBlender registered")
    print("HexBlender Addon found: ", __file__)

    hexblender_info["hexblender_addon_path"] = os.path.dirname(__file__)
    print("HexBlender Addon Path is ",
          hexblender_info["hexblender_addon_path"])
    addon_path = os.path.dirname(__file__)


def unregister():

    bpy.utils.unregister_module(__name__)
    del bpy.types.Scene.hexblender

    print("HexBlender unregistered")


# for testing
if __name__ == '__main__':
    print("Came6")
    register()
