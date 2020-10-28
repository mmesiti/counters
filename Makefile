flop_count.py: flop_count_template.py.c
	cpp -DPYTHON flop_count_template.py.c | grep -Ev "^\s*;" > flop_count.py


flop_count.c: flop_count_template.py.c
	cpp -DNF=4 flop_count_template.py.c > flop_count.c
