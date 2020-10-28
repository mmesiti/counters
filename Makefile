
SRCS=flop_count_template.py.c memory_count_template.py.c transfer_count_template.py.c

CTARGETS=$(patsubst %_template.py.c,%.c,$(SRCS))
PYTARGETS=$(patsubst %_template.py.c,%.py,$(SRCS))
OTESTTARGETS=$(patsubst %_template.py.c,%.test.o,$(SRCS))

TARGETS=$(CTARGETS) $(PYTARGETS) $(OTESTTARGETS)

all: $(TARGETS)

%.py : %_template.py.c memory_base.py.h
	cpp -DPYTHON $< | grep -Ev "^\s*;" | ((which yapf &> /dev/null && yapf) || cat )  > $@

%.c : %_template.py.c memory_base.py.h
	cpp $< | ((which clang-format &> /dev/null && clang-format) || cat )  > $@

%.test.o: %.c
	gcc -o $@ -c -DNF=4 -DT=8 -DX=8 -DY=8 -DZ=8 $<

python: $(PYTARGETS)

c: $(CTARGETS)

testo: $(OTESTTARGETS)

clean:
	rm -f $(TARGETS)

$(TARGETS) :
