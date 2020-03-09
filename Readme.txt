================================================================
Some usage notes on the HexBlender addon
================================================================

The HexBlender scripts are being developed using Blender 2.77.

=============================
Installing HexBlender addon
=============================
To install HexBlender as a Blender addon:
1) All the Hexblender code can be found within Continuity at src/blender/blender_2.7.  Run the make_install_file.py which will create a zip file named hexblender.zip, or simply use the existing hexblender.zip file if it is up to date.  Note, you can create that zip file with Python 2.7, it doesn't need to be Python 3 even though Blender 2.77 uses Python 3.5.1.

2) You will need to install scipy within Blender's Python before using Hexblender:

3) hexblender.zip includes all the files you need as a Blender add-on.

########### Mac #############
On Mac systems we have been successful in simply copying the scipy directory from Anaconda (http://continuum.io/downloads -- be sure to get the Python 3.5 version), into the appropriate place within the Blender installation.  Something like the following worked for us:
Mac: cp -r /Anaconda/lib/python3.5/site-packages/scipy /Applications/Blender_2.77a/blender.app/Contents/Resources/2.77/python/lib/python3.5/site-packages/

Alternatively, you can simply download Anaconda's scipy directory from here:
https://www.dropbox.com/s/w38yb331e64k2sl/scipy.zip?dl=0
and copy it into the site-packages directory within the Blender installation.

########### Windows #############
On Windows 64bit systems we include the site-packages directory from Anaconda (http://continuum.io/downloads) into our PYTHONPATH definition;
Instructions:
1. Download and install Anaconda
2. Open up Anaconda Prompt as an administrator
3. Type "conda create --name blender python=3.5" and press enter
4. Type "conda activate blender" and press enter
5. Type "pip install numpy" and press enter
6. Type "pip install scipy" and press enter
7. Set your PYTHONPATH environment variable to equal "<base folder path for anaconda>/envs/blender"

To do this, you'll need to add PYTHONPATH as a new system environmental variable:
For Windows 8 and 10, search for "system variables" in the Control Panel.
For Windows 7, Start -> Right Click on Computer -> Properties -> Advanced -> Environmental Variables.

########### Linux #############
On Linux 64bit systems we have been successful in putting the site-packages directory from Anaconda (http://continuum.io/downloads -- be sure to get the Python 3 64bit version) into our PYTHONPATH definition.
Anaconda can do this during the installation or you can do it yourself afterwards. Something like the following worked for us:
export PYTHONPATH=$HOME/anaconda3/lib/python3.4/site-packages:$PYTHONPATH
After doing that Hexblender did load and work successfully.

Alternatively, we have provided a prebuilt scipy module for people to install. This was built on Ubuntu 12.04 LTS:
http://cmrg.ucsd.edu/wiki/widget/data/scipy_0_14_0_ubuntu_12_04.tar.gz

3) From within Blender, go File -> User Preferences -> Addons, and select the button at the bottom: Install from File and select the zip file made in step one. 

4) Type HexBlender in the little search window on the Addon windows and you should then see the HexBlender addon.  Click on the box next to the dancing guy to activate it.  You will likely want to then hit the button at the botton "Save User Settings", so HexBlender is always activated when starting Blender.


================
Using HexBlender
================
The current thinking is that most of the HexBlender tools will be located in the Properties window, and more specifically in the Scene and Object panels of the Properties window.

Working Functions:                         Where located:
Import Pickle File                       on the Scene panel
Export Nodes                             on the Object panel under Export
Export Elements                          on the Obj panel under Export
Export Hermite Tricubic Derivs           on the Obj panel under Export 
Export Face Keys                         on Obj panel under Export
Regularize Hexs                          on the Obj panel under Hex Tools
Regularize Surface                       on the Obj panel under Quad Face Tools
Linear Nodal Adjustment                  on the Obj panel under Hex Tools
Circlular Nodal Adjustment               on the Obj panel under Hex Tools
Std Dev of Edge Lengths                  on the Obj panel under Selection Tools
Print Selection Data                     on the Obj panel under Selection Tools
Update/Pickle Hexs                       on Obj panel under Hex Tools
Std Dev of Faces                         on Obj panel under Selection Tools
Std Dev of Angle                         on Obj panel under Selection Tools
Show Vertex                              on Obj panel under Selection Tools
Show Edge                                on Obj panel under Selection Tools
Show Face                                on Obj panel under Selection Tools
Rotate Hex                               on Obj panel under Hex Tools
Export Fit Data Points                   on Obj panel under Export
Compute tensor transformation            on Obj panel under Tensor Tools
Show Hex                                 on Obj panel under Selection Tools
Visualize Indices                        on Obj panel under Selection Tools


There are also additional buttons that are currently displayed in the GUI but are not functional.  They will report "Not yet implemented" when executed.


===============
Known Issues
===============
None at this time.

======================================
Reporting Issues or run into problems?
======================================
Email jvandorn at ucsd dot edu
