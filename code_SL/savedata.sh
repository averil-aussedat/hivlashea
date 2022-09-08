# how to call it : ./savedata.sh yaml/theyamlfile.yaml <newfoldername>

if [ -z "$1" ]
then
  echo No yaml file provided. Abort
else
  if [ -z "$2" ]
  then
    newfolder=run
  else
    newfolder=run_$2
  fi
  echo Saving datas into $newfolder...
  mkdir $newfolder
  cp -r data_input $newfolder
  cp -r data_output $newfolder
  cp -r python_diags $newfolder
  cp $1 $newfolder
  echo ... Done.
fi

