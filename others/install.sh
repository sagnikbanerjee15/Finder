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

cd dep
tar -xvzf gmst_linux_64.tar.gz
tar -xvzf gmes_linux_64.tar.gz
chmod -R a+x *
cd gmes_linux_64
perl change_path_in_perl_scripts.pl "/usr/bin/env perl"
cd ..
gunzip gm_key_64.gz

# Remove the downloaded tar.gz files
rm gm*.tar.gz

# The key needs to be in the home directory
mv gm_key_64 .gm_key
mv .gm_key ~
cd ..

# Make everything executable
chmod -R a+x *