*.pdf : *.dat Plot.py
	python Plots.py

*.dat : sample_c
	./sample_c

sample_c : sample.c
	gcc -fopenmp sample.c -lm -o sample_c

clean :
	rm *.dat  *_c