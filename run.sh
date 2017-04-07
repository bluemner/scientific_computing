#clear && make && ./bin/driver.o
clear && make &&
./bin/ReadWritePGM.o ./data/ship.pgm ./bin/linux_ship_01.pgm &&
./bin/ReadWritePGM.o ./data/street.pgm ./bin/linux_street_01.pgm &&
./bin/ReadWritePGM.o ./data/woman.pgm ./bin/linux_woman_01.pgm 

