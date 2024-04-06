
# The easiest thing to do to install Seurat AND Signac would be to, unfortunately, builk a new R and immediately install Seurat, Signac, and all of its dependencies. 
# Seurat requires a new version of GCC while other packages like igraph require older versions of gcc
# Ex: igraph needs libgfortran.so.4 which only the old version of gcc has
# Simply run this or add to your bash profile to add the new gcc to LD_Library_PATH
# To check whether the gcc version has the GLIB you're looking for: strings /software/gcc-11.2.0/lib64/libstdc++.so.6.0.29 | grep "GLIB"
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:"/frazer01/home/jennifer/anaconda3/lib":"/software/gcc-11.2.0/lib64":"/software/gcc-7.2.0/lib64":"/software/libpng-1.6.35/lib":"/export/software/gcc-7.2.0/lib64/libgfortran.so.4":"/usr/bin/g++":"/software/mpfr-2.4.2/lib":"/usr/lib64/libmpfr.so.4"

# if you run into any errors related to a "lib" or "GLIB" while installing Seurat, you can fix by:
# 1. Find the path of the library/software/tool
locate [software_name]

# 2. Then add to the LD_LIBRARY_PATH above like I did with the other libraries and try installing Seurat again

# If you are unable to install Rcpp, spatstat, reticulate (which are very problematic), try these:
install.packages(“https://cran.r-project.org/src/contrib/Archive/RcppTOML/RcppTOML_0.1.3.tar.gz”, repos=NULL, type=“source”)
install.packages(“https://cran.r-project.org/src/contrib/Archive/spatstat/spatstat_1.64-1.tar.gz”, repos=NULL, type=“source”)
install.packages("https://cran.r-project.org/src/contrib/Archive/reticulate/reticulate_1.4.tar.gz", repos=NULL, type="source")

