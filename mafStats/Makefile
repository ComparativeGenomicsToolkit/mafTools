include ../inc/common.mk
SHELL:=/bin/bash
bin = ../bin
inc = ../inc
lib = ../lib
PROGS = mafStats
dependencies = ${inc}/common.h ${inc}/sharedMaf.h ${lib}/common.c ${lib}/sharedMaf.c
objects := ${lib}/common.o ${lib}/sharedMaf.o ../external/CuTest.a src/test.mafStats.o ${sonLibPath}/sonLib.a src/buildVersion.o
testObjects := test/sharedMaf.o test/common.o ../external/CuTest.a src/test.mafStats.o ${sonLibPath}/sonLib.a test/buildVersion.o
sources = src/mafStats.c src/mafStats.h

.PHONY: all clean test buildVersion

all: buildVersion $(foreach f,${PROGS}, ${bin}/$f)
buildVersion: src/buildVersion.c
src/buildVersion.c: ${sources} ${dependencies}
	@python ../lib/createVersionSources.py

../lib/%.o: ../lib/%.c ../inc/%.h
	cd ../lib/ && make

${bin}/mafStats: src/mafStats.c ${dependencies} ${objects}
	mkdir -p $(dir $@)
	${cxx} $< ${objects} -o $@.tmp ${cflags} ${lm}
	mv $@.tmp $@

test/mafStats: src/mafStats.c ${dependencies} ${testObjects}
	mkdir -p $(dir $@)
	${cxx} $< src/allTests.c ${testObjects} -o $@.tmp ${testFlags} ${lm}
	mv $@.tmp $@
%.o: %.c %.h
	${cxx} -c $< -o $@.tmp ${cflags}
	mv $@.tmp $@
test/%.o: ${lib}/%.c ${inc}/%.h
	@echo apple
	mkdir -p $(dir $@)
	${cxx} -c $< -o $@.tmp ${testFlags} ${lm}
	mv $@.tmp $@
test/%.o: src/%.c src/%.h
	@echo orange
	mkdir -p $(dir $@)
	${cxx} -c $< -o $@.tmp ${testFlags}
	mv $@.tmp $@

clean:
	rm -rf $(foreach f,${PROGS}, ${bin}/$f) src/*.o test/ src/buildVersion.c src/buildVersion.h

test: buildVersion test/allTests
	test/allTests && rm -rf ./test/

test/allTests: src/allTests.c ${testObjects}
	mkdir -p $(dir $@)
	${cxx} $^ -o $@.tmp ${testFlags} ${lm}
	mv $@.tmp $@

../external/CuTest.a: ../external/CuTest.c ../external/CuTest.h
	${cxx} -c $< ${cflags}
	ar rc CuTest.a CuTest.o
	ranlib CuTest.a
	rm -f CuTest.o
	mv CuTest.a $@
