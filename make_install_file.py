""" Simple script the will just create a zip file that can be used to
    easily install the HexBlender addon into Blender """
import os
import shutil
import zipfile

dirname = "hexblender"
zip_filename = "hexblender.zip"

# make a new directory
try:
    os.mkdir("hexblender")
except Exception as msg:
    print("ERROR: Unable to make directory: %s.  \n%s" % (dirname, msg))

# copy all python files into new directory
try:
    all_files = os.listdir(".")
    for f in all_files:
        if f.endswith(".py") and "make_install" not in f:
            shutil.copy(f, dirname)
        if "REVISION.txt" in f:
            shutil.copy(f, dirname)
except Exception as msg:
    print("ERROR: Unable to copy files.  \n%s" % msg)

# make the zip file
zipf = zipfile.ZipFile(zip_filename, 'w')
for root, dirs, files in os.walk(dirname):
    for f in files:
        zipf.write(os.path.join(root, f))
zipf.close()

# remove the directory
try:
    shutil.rmtree(dirname)
except Exception as msg:
    print("ERROR: Unable to remove directory: %s. \n%s" % (dirname, msg))
