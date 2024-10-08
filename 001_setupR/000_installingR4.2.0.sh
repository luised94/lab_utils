# Installing R 4.2.0
# Purpose: Create a working environment that is similar to linux cluster to minimize the chance of conflict.
# See Brain Mode short version for thread

sudo apt-get update
sudo apt install r-base
sudo apt install build-essential libcurl4-gnutls-dev libxml2-dev libssl-dev

R_VERSION="R-4.2.0"
DIR_TO_INSTALL=$HOME
curl --output "$HOME/${R_VERSION}.tar.gz" "https://cran.r-project.org/src/base/$(echo ${R_VERSION} | cut  -d. -f1)/${R_VERSION}.tar.gz"

# Install Java
sudo apt-get install default-jdk
export JAVA_HOME=/usr/lib/jvm/java-11-openjdk-amd64
#echo "export JAVA_HOME=/usr/lib/jvm/java-11-openjdk-amd64" >> ~/my_config/dotfiles/bashrc
export PATH=$JAVA_HOME/bin:$PATH
#sudo R CMD javareconf

tar -xzvf ${R_VERSION}.tar.gz
cd ${R_VERSION}


./configure --enable-R-shlib --with-blas --with-lapack --without-x
make
sudo make install

# Install for PDF and html rendenring of manuals
# sudo apt-get install texinfo
# sudo apt-get install texlive texlive-fonts-recommended texlive-latex-recommended texlive-latex-extra

# Install X11 headers and libs
# sudo apt-get install libx11-dev xserver-xorg-dev xorg-dev

export PATH=/usr/local/bin:$PATH
#echo "export PATH=/usr/local/bin:$PATH" >> ~/.bashrc
which R # should show R version 4.2.0

mkdir -p $HOME/R/library





export R_LIBS_USER=$HOME/R/library
alias R='R --no-save'
R

Sys.getenv(HOME)
dir.create(Sys.getenv("R_LIBS_USER"), recursive = TRUE)
.libPaths(Sys.getenv("R_LIBS_USER"))
.libPaths()
options(repos = c(CRAN = "https://cloud.r-project.org"))
install.packages(c("renv", "xml2", "lintr", "roxygen2", "languageserver"), dependencies = TRUE, INSTALL_opts = '--no-lock')
install.packages(c("languageserver"), dependencies = TRUE, INSTALL_opts = '--no-lock')
q()
