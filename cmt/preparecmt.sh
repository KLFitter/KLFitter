#!/bin/bash

if [[ ! -h $PWD/../KLFitter ]]; then ln -s $PWD/../library/include $PWD/../KLFitter; fi
if [[ ! -h $PWD/../src      ]]; then ln -s $PWD/../library/src     $PWD/../src;      fi
#if [[ ! -h $PWD/../examples/top_ljets/runKLFitter.cxx ]]; then ln -s $PWD/../examples/top_ljets/runKLFitter.c $PWD/../examples/top_ljets/runKLFitter.cxx; fi

#for i in $PWD/../KLFitter/*.h; do
#  j=${i/$PWD\/..\/KLFitter\//}
#  for k in $PWD/../KLFitter/*.h; do
#    echo $k | xargs perl -p -i -e "s/\"$j\"/\"KLFitter\/$j\"/"
#  done
#  for k in $PWD/../src/*.cxx; do
#    echo $k | xargs perl -p -i -e "s/\"$j\"/\"KLFitter\/$j\"/"
#  done
#done
#
#for i in $PWD/../KLFitter/*.h; do
#  j=${i/$PWD\/..\/KLFitter\//}
#  for k in $PWD/../extras/include/*.h; do
#    echo $k | xargs perl -p -i -e "s/\"$j\"/\"KLFitter\/$j\"/"
#  done
#  for k in $PWD/../extras/src/*.cxx; do
#    echo $k | xargs perl -p -i -e "s/\"$j\"/\"KLFitter\/$j\"/"
#  done
#done
#
#for i in $PWD/../extras/include/*.h; do
#  j=${i/$PWD\/..\/extras\/include\//}
#  for k in $PWD/../extras/include/*.h; do
#    echo $k | xargs perl -p -i -e "s/\"$j\"/\"..\/include\/$j\"/"
#  done
#  for k in $PWD/../extras/src/*.cxx; do
#    echo $k | xargs perl -p -i -e "s/\"$j\"/\"..\/include\/$j\"/"
#  done
#done

for i in $PWD/../KLFitter/*.h; do
  j=${i/$PWD\/..\/KLFitter\//}
  echo $j
  echo ../examples/top_ljets/runKLFitter.c | xargs perl -p -i -e "s/\"$j\"/\"KLFitter\/$j\"/"
done

for i in $PWD/../extras/include/*.h; do
  j=${i/$PWD\/..\/extras\/include\//}
  echo $j
  echo ../examples/top_ljets/runKLFitter.c | xargs perl -p -i -e "s/\"$j\"/\"..\/..\/extras\/include\/$j\"/"
done

exit 0
