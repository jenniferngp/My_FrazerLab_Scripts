
# 1. Make sure conda is directed to yours!
which conda 

# 2. Make sure curl is directed to /usr/bin/curl! The default is your conda's curl which does not work for R installation.
which curl
export PATH=/usr/bin/:$PATH
which curl

# 3. Download an R.tar.gz file (https://cran.r-project.org/mirrors.html) - choose any mirror
wget http://lib.stat.cmu.edu/R/CRAN/src/base/R-4/R-4.0.1.tar.gz

# 4. Unzip 
tar -xf R-4.0.1.tar.gz
cd R-4.0.1

# 5. Configure
# If you get "--with-x=yes (default) and X11 headers/libs are not available" error, it could be that you're on a node and not on the head node
# If you get "--with-readline" error, use "--with-readline=no --with-x=no" to your configure command
./configure --with-pcre1 --prefix=/frazer01/home/jennifer/software/R-4.0.1

# 6. Build
make && make install

# 7. Export or add to your bash profile
export PATH=~/software/R-4.0.1/bin:$PATH

# 8. Open R to check that it works and the version is correct
R

# 9. If you are doing single-cell analyses, HIGHLY HIGHLY recommend installing Seurat and Signac IMMEDIATELY. 
# They have specific dependencies and often clash with other packages. 
# Installing them down the line will cause a lot of problems and conflicts with other packages. 
