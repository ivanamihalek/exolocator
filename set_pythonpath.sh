# source this to add the path to PYTHONPATH
myname=${BASH_SOURCE[0]}
fullpath=$(realpath  "$myname")  # get the full path to us
mydir=$(dirname "$fullpath")  # directory only

if [ -z $PYTHONPATH ]
then
      # if PYTHONPATH not set yet, than it's mydir
      export PYTHONPATH=$mydir
else
      # delete all appearance of :$mydir from the existing path
      PYTHONPATH=${PYTHONPATH//:$mydir/}
      #  if $mydir already first in the PYTHONPATH do nothing
      if [[ $PYTHONPATH == "$mydir" || $PYTHONPATH =~ "$mydir:".* ]]
      then # do nothing
          :
      else # prepend :$mydir from the existing path
          export PYTHONPATH=$mydir:$PYTHONPATH
      fi
fi
echo PYTHONPATH is  set to $PYTHONPATH
