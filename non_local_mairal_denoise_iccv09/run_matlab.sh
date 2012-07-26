
   #!/bin/sh
   # if you can not use the mex files, uncomment this line.
   export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:./libs/
   export KMP_DUPLICATE_LIB_OK=true
   matlab -nodesktop -nosplash 
