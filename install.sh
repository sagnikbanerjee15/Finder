#!/bin/bash

# Install OLego
cd $PWD/dep/olego
make
chmod -R a+x *
cd ../..

# Install PsiCLASS
cd $PWD/dep/assemblies_psiclass_modified
make 
chmod -R a+x *
cd ../..

# Install AUGUSTUS
cd $PWD/dep/Augustus
chmod -R a+x *
cd ../..

# Install BRAKER2
cd $PWD/dep/BRAKER
chmod -R a+x *
cd ../..

# Install GUSHR
cd $PWD/dep/GUSHR
chmod -R a+x *
cd ../..