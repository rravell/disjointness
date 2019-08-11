#!/bin/bash

DIRECTORY="$HOME/git/disjointness/src"
DIRTYSCRIPT="dirtyscript.lisp"
VECTORFILE="IneffFunctionals44Lisp.txt"
OUTPUTFILE="NPAvalue.txt"

cat << EOF > $DIRECTORY/$DIRTYSCRIPT
#!/usr/bin/sbcl --script
(load "$HOME/quicklisp/setup.lisp")
(ql:quickload :npa-hierarchy)
(in-package :npa-user)
EOF

while read -r vector; do
	echo "(print (solve-problem (maximise $vector) (level 2)(scenario (4 4 4 4) (4 4 4 4))))" >> $DIRECTORY/$DIRTYSCRIPT
done < $DIRECTORY/$VECTORFILE

rm .npa_sdpa_tmp.dat-s
chmod +x $DIRECTORY/$DIRTYSCRIPT > $DIRECTORY/$OUTPUTFILE
$DIRECTORY/$DIRTYSCRIPT > $DIRECTORY/$OUTPUTFILE