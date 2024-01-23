for p in proc*
do
  echo $p
  cd $p
  rm -r [1-9]*
  cd ..
done
