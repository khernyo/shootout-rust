
.PHONY: all

TARGETS := $(basename $(wildcard *.rs))

all: $(TARGETS)

%: %.rs
	rustc $< -o $@
