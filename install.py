#! /usr/bin/env python

import os

"""
Check the bashrc file to check if finder was previously installed
"""
temp_file=os.path.expanduser('~')+"/"+ '.bashrc.temp'
actual_file=os.path.expanduser('~')+"/"+ '.bashrc'

finder_installation_locations=[]
fhr=open(os.path.expanduser('~')+"/"+ '.bashrc',"r")
for line in fhr:
    if "export" in line and "Finder" in line:
        finder_installation_locations.append(line.strip().split("$PATH:")[-1])
fhr.close()

#print(finder_installation_locations)
remove_these_indices=[]
for i,each_installation in enumerate(finder_installation_locations):
    check_this_file="software_identity"
    verify_contents="FINDER: An automated software package to annotate eukaryotic genes from RNA-Seq data and associated protein sequences - Banerjee et, al 2021"
    #print(each_installation+"/"+check_this_file)
    if os.path.exists(each_installation+"/"+check_this_file)==True:
        #print(open(each_installation+"/"+check_this_file).read().split("\n")[0])
        #print(verify_contents)
        if verify_contents != open(each_installation+"/"+check_this_file).read().split("\n")[0]:
            remove_these_indices.append(i)
        else:
            print("Match found")
    else:
        print("File not found")
        remove_these_indices.append(i)

#print(remove_these_indices)
for i in remove_these_indices[::-1]:
    finder_installation_locations.pop(i)
    

#print(finder_installation_locations)
if len(finder_installation_locations)>0:
    fhr=open(os.path.expanduser('~')+"/"+ '.bashrc',"r")
    fhw = open(os.path.expanduser('~')+"/"+ '.bashrc.temp',"w")
    for line in fhr:
        if "export" in line and "Finder" in line:
            if line.strip().split("$PATH:")[-1] in finder_installation_locations:
                loc=line.strip().split("$PATH:")[-1] 
                print(f"{loc} will be removed from $PATH")
                continue
            else:
                fhw.write(line)
        else:
            fhw.write(line)
    fhw.close()
    
    cmd=f"mv {temp_file} {actual_file}"
    os.system(cmd)
    
    cmd=f"rm {temp_file}"
    os.system(cmd)
    
pwd = os.getcwd()

# Install OLego
os.chdir(pwd+"/dep/olego")
os.system("make")
os.system("chmod -R a+x *")
os.chdir(pwd)

# Install PsiCLASS
os.chdir(pwd+"/dep/psiclass_terminal_exon_length_modified")
os.system("make")
os.system("chmod -R a+x *")
os.chdir(pwd)

# Install Augustus
os.chdir(pwd+"/dep/Augustus")
os.system("chmod -R a+x *")
os.chdir(pwd)

# Install BRAKER2
os.chdir(pwd+"/BRAKER")
os.system("chmod -R a+x *")
os.chdir(pwd) 

# Install GUSHR
os.chdir(pwd+"/GUSHR")
os.system("chmod -R a+x *")
os.chdir(pwd)

os.chdir(pwd+"/dep")
os.system("tar -xvzf gmst_linux_64.tar.gz")
os.system("tar -xvzf gmes_linux_64.tar.gz")
os.system("chmod -R a+x *")
os.chdir(pwd+"/dep/gmes_linux_64")
os.system("perl change_path_in_perl_scripts.pl \"/usr/bin/env perl\"")
os.chdir(pwd+"/dep")
os.system("gunzip gm_key_64.gz")

# Remove the downloaded tar.gz files
os.system("rm gm*.tar.gz")

# The key needs to be in the home directory
os.system("mv gm_key_64 .gm_key")
os.system("mv .gm_key ~")
os.chdir(pwd)

# Make everything executable
os.system("chmod -R a+x *")






