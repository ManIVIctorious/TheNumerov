
.PHONY: all
all: Numerov Watson

.PHONY: Numerov
Numerov:
	@cd TheNumerov && make --no-print-directory Numerov

.PHONY: Watson
Watson:
	@cd TheWatson  && make --no-print-directory Watson


.PHONY: clean
clean:
	@cd TheNumerov && make clean
	@cd TheWatson  && make clean
