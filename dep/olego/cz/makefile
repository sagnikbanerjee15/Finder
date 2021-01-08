include ../common.mk
XINC = -I../kent

O = BedFile.o utils.o

%.o: %.cpp
	${CC} ${COPT} ${CFLAGS} ${LOWELAB_DEFS} ${HG_DEFS} ${HG_WARN} ${XINC} -o $@ -c $<

lib:czlib.a

czlib.a: $(O)
	ar rcus $@ $(O)

test:
	cd tests && ${MAKE} test

cleanlocal:
	rm -rf ${O} czlib.a

clean:cleanlocal

#cd tests && ${MAKE} clean
