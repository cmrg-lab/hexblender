bl_info = {
    "name": "export mesh data to TXT",
    "author": "PATRICK BOELEN",
    "version": (1, 0),
    "blender": (2, 7, 0),
    "location": "View3D > Tool Shelf",
    "description": "EXPORT MESH DATA TO TXT",
    "warning": "",
    "wiki_url": "",
    "category": "Import-Export"}


import bpy

class exportToTXT(bpy.types.Operator):
    bl_idname = "export.export_to_txt"
    bl_label =  "Export To TXT"
    
    filepath = bpy.props.StringProperty(subtype="FILE_PATH")
    
    def execute(self, context):
        obverts = bpy.context.active_object.data.vertices
        obfaces = bpy.context.active_object.data.polygons

        verts = []
        faces = []
        

        
        for vertex in obverts:
            verts.append(tuple(vertex.co))
            
            
        for face in obfaces:
            faces.append(tuple(face.vertices))
            
            
        file = open(self.filepath, 'w')
        file.write(str("datos de Vertices:   "))
        file.write(str("\n\n\n"))
        file.write(str(verts))
        file.write(str("\n\n\n\n\n\n\n\n\n\n\n\n"))
        file.write(str("datos de Caras:   "))
        file.write(str("\n\n\n"))
        file.write(str(faces))

        
        
            
        return {'FINISHED'}
    
    
    def invoke(self, context, event):
        context.window_manager.fileselect_add(self)
        return {'RUNNING_MODAL' }
    
   
class exportToTXTPanel(bpy.types.Panel):
    bl_idname = "Export_To_TXT"
    bl_label =  "Export To TXT"
    bl_space_type = "VIEW_3D"
    bl_region_type = "TOOLS"
    bñ_context = "objectmode"
    
    def draw(self, context):
        layout = self.layout
        layout.operator("export.export_to_txt", text="Export to TXT")
    
    
            
def register():
    bpy.utils.register_module(__name__)

def unregister():
    bpy.utils.unregister_module(__name__)


if __name__ == "__main__":
    register()