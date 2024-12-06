#STATUS:
sudo apt-get remove r-base r-base-core r-recommended
sudo apt-get autoremove
sudo apt-get purge r-base r-base-core r-recommended
sudo rm -rf /usr/local/lib/R
sudo rm -rf /usr/lib/R
sudo rm -rf /etc/R
rm -rf ~/.R
rm -rf ~/R
rm -rf ~/R-4.2.0/
rm -rf ~/.cache/R
sudo apt-get update
sudo apt-get clean
sudo apt-get autoclean
rm -rf ~/R/x86_64-pc-linux-gnu-library
sudo find / -name "*R-*" -type d
