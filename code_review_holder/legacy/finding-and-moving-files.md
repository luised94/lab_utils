#Finding and moving files in current working directory 

>Create destination folder with desired name. Use find and xargs to move files. 
'$mkdir PDFs'

'$find -iname '*.pdf' -print0 | xargs -0 mv -t ./PDFs/'


> Prints file extension if any, sorts and outputs only unique files 
'$find . -type f | perl -ne 'print $1 if m/\.([^.\/]+)$/' | sort -u'

> Using awk
find . -type f | awk -F. '!a[$NF]++{print $NF}'

>Without awk, perl, sed 
find . -type f | rev | cut -d. -f1 | rev  | tr '[:upper:]' '[:lower:]' | sort | uniq | sort -rn

> My attempt, doenst work that well  
find . -type f -exec file {} \; | cut -d: -f 2 | sort -u

find . -type f -exec file {} \; | cut -d: -f 2 | sed 's/^ //; s/ +/ /g' | sort -u

> Pipe any of the results into a text file then use it to create directories.
mkdir fileextensions.txt

>Run the following for loop. Almost destroyed my computer probably by not specifying the current working directory. Unfortunately, command goes through all files, even if they are sorted. 

'for filename in ./**/*; do
  if [[ -f "$filename" ]]; then
      base=${filename%.*}
      ext=${filename#$base.}
    mv "$filename" "${ext}"
  fi
done'

> Repeat for files two directories deep. 
'for filename in ./*; do
  if [[ -f "$filename" ]]; then
      base=${filename%.*}
      ext=${filename#$base.}
    mv "$filename" "${ext}"
  fi
done'

> Delete the remaining empty directories 
'$find -type d -empty'

> Actually still had more levels of directories. Gave same message because it goes through all files
'shopt -s globstar dotglob'
'for filename in ./**/**; do
  if [[ -f "$filename" ]]; then
      base=${filename%.*}
      ext=${filename#$base.}
    mv "$filename" "${ext}"
  fi
done'
'$find -type d -empty -delete'

'$rm -r dir1 dir2 dir3'

> After all this, certain files were not moved and some files were created with no extension
> Moved these manually.


>Pretty sure was what I was looking for. file exits with positive status if it finds something
if [[ -n $(find /var/log/crashes -name "app-*.log" -mmin -5) ]]
then
    service myapp restart
fi
