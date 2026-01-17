cmake -DCMAKE_INSTALL_PREFIX=../install \
      -DDD4HEP_BUILD_PACKAGES="DDRec;DDDetectors;DDCond;DDAlign;DDDigi;DDG4;DDEve;UtilityApps" \
      -DDD4HEP_USE_LCIO=ON \
      -DDD4HEP_USE_GEANT4=ON ..

make -j8
make install