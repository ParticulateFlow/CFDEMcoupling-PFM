#in controlDict startFrom is "latestTime".
#The first loop is just needed for time directories smaller than 1.

for number in {0..100..1}
do

 if (($number % 100 == 0))
 then
 da=$(bc <<< "($number/100)")
 mkdir "$da"
 mapFields /path/to/source -sourceTime $da -consistent 
 elif (($number % 10 == 0))
 then
 da=$(bc <<< "scale=1;($number/100)")
 mkdir "0$da"
 mapFields /path/to/source -sourceTime $da -consistent
 else
 da=$(bc <<< "scale=2;($number/100)")
 mkdir "0$da"
 mapFields /path/to/source -sourceTime $da -consistent 

 fi
done

for number in {101..2500..1}
do

 if (($number % 100 == 0))
 then
 da=$(bc <<< "($number/100)")
 mkdir "$da"
 mapFields /path/to/source -sourceTime $da -consistent 
 elif (($number % 10 == 0))
 then
 da=$(bc <<< "scale=1;($number/100)")
 mkdir "$da"
 mapFields /path/to/source -sourceTime $da -consistent
 else
 da=$(bc <<< "scale=2;($number/100)")
 mkdir "$da"
 mapFields /path/to/source -sourceTime $da -consistent 

 fi
done

