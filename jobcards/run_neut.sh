#!/bin/bash

#run neut
neutroot2 NEUT_6.1.2_CC_numu.toml samples/neut.ghep.root

## Only one NUISANCE step is required for NEUT
#echo "Creating NUISANCE flat trees"
nuisflat -f GenericVectors -i NEUT:samples/neut.ghep.root -o samples/neut.flat.root

rm hEnumu_cv_o.root
