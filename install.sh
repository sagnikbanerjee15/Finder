#!/bin/bash

# Install OLego
cd $PWD/src/olego
make
chmod -R a+x *
cd ../..

# Install PsiCLASS
cd $PWD/src/assemblies_psiclass_modified
make 
chmod -R a+x *
cd ../..

# Install AUGUSTUS
cd $PWD/src/Augustus
chmod -R a+x *
cd ../..

# Install BRAKER2
cd $PWD/src/BRAKER
chmod -R a+x *
cd ../..

# Install GUSHR
cd $PWD/src/GUSHR
chmod -R a+x *
cd ../..