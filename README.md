SAGECAL on Daliuge
===========
#### python files:
  if source root directory of daliuge is $DFMS, first put all python files of python directory into "$DFMS/deploy/cnlab" directory.
  then reinstall Daliuge.

#### compile
    1, create a soft link "lib" to the directory ligsagecal.a existed(ln -s).
	2, update Makefile if dependent libs were installed in different directories (like casacore, openBLAS).
	3, make all

    
#### dist_sagecal.json
  It is the logical graph and includes self-customized paths(like BashShellDrop's cmd) which you need to modify according to your path.