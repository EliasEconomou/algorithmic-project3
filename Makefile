
PROGRAMS = $(subst programs/, , $(wildcard programs/*))

# Compile all programs
all: programs tests

programs-%:
	$(MAKE) -C programs/$*

programs: $(addprefix programs-, $(PROGRAMS))

tests:
	$(MAKE) -C tests all

run-programs: $(addprefix run-programs-, $(PROGRAMS))

run-tests:
	$(MAKE) -C tests run

# Clean all programs
clean-programs-%:
	$(MAKE) -C programs/$* clean

clean: $(addprefix clean-programs-, $(PROGRAMS))
	$(MAKE) -C tests clean