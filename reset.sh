rm picdata.bin
rm pictarget.bin
rm PICit
rm *.dat

if (($1 == 1))
then
	./comp_PICit
	./PICit 0
fi
