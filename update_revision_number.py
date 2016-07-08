'''
We need a way to easily know what revision we are at and to be able to
update the revision to the proper value.
We previously used SVN and just used an incremental checkin id to be
our revision number.  We will roughly still follow that, but as Git uses
hashes for each commit, we ideally need to track that information as well.
Blender only allows integers in it's bl_info:version setting.  So, we'll
create a new file that will contain the most recent hash.

The commit procedure should be:
  1. Commit all of the changes to the normal code
  2. Run this script
  3. Run make_install_file.py
  4. Commit __init__.py, REVISION.txt and hexblender.zip

Steps 2 and 3 can easily be combined should that become necessary.

For now, our revision number will simply be: MajorVersion.NumberOfCommits

Script Requirement: the git executable needs to be in your PATH
'''
import os
import shutil
from subprocess import check_output

TEMP_FILENAME = "asdofiasdj.txt"
INIT_FILENAME = "__init__.py"

MAJOR_VERSION = 1

# Get number of commits
output = check_output("git log --oneline", shell = True )
num_of_commits = len(output.splitlines())

# Increase the number to account for the next commit
# which will contain __init__.py and hexblender.zip
new_commit_value = num_of_commits + 1

# Get the latest hash and store in revision file
# NOTE: This will get hash for whatever branch you are on
last_commit_hash = str(output[0:7], 'utf-8')

full_revision_number = "%d.%s.r%s" % (MAJOR_VERSION, new_commit_value, last_commit_hash)
print("New revision is: %s" % full_revision_number)

# As the revision file currently only stores the revision information,
# just open and write to the file, replacing the entire contents
try:
    f = open("REVISION.txt", "wt")
    f.write(full_revision_number)
    f.close()
except Exception as msg:
    print("Error updating REVISION.txt file!!! : %s" % msg)
    exit(1)

# Now update the __init.py__ file
try:
    output = open (TEMP_FILENAME, 'wt')
    with open(INIT_FILENAME, 'r') as file:
        for line in file:
            if '"version"' in line:
                line = '    "version": (%d,%d),\n' %(MAJOR_VERSION, new_commit_value)
            output.write(line)

    output.close()

    # copy the file
    shutil.copy(TEMP_FILENAME, INIT_FILENAME)
    print("__init__.py has been updated.")

    #remove temp file
    os.remove(TEMP_FILENAME)

except Exception as msg:
    print("Error updating __init__.py file!!! : %s" % msg)
