#STATUS:
# git config --global user.email "you@example.com"
# git config --global user.name "Your Name"

#Using system to run the git config command to 
git_config_email <- "git config user.email"
git_config_name <- "git config user.name"
email <- "liusmartinez94@gmail.com"
username <- "Luis"
paste(git_config_email, email)
system(paste(git_config_email, email))
system(paste(git_config_name, username))
