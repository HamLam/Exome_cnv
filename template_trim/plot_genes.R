options(warn=-1);
.libPaths("/home/msistaff/lamx0031/R/x86_64-unknown-linux-gnu-library/3.1")
library(RMySQL)
library(calibrate)
library(plotrix)
library(zoo)

m <- dbDriver("MySQL")
con <-dbConnect(m,username="root",dbname="cnv",host="localhost",unix.socket="socket_path");
dir_path = "sample_path";
setwd(dir_path);
source("scripts_location/R_function.R");
