# install packages needed from CRAN
list.of.packages <- c("proto","stringr","XML","xtable","randtoolbox")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]

if (length(new.packages)) {
  print("installing : ")
  print(new.packages)
  install.packages(new.packages, repos = "http://cran.rstudio.com/")
}
