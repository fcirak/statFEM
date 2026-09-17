SUBDIRS := del2/apps/reference/laplace

.PHONY: all $(SUBDIRS)

all:
	@echo "Choose a target using VSCode Makefile Tools."

# Each subdir builds itself
$(SUBDIRS):
	$(MAKE) -C $@

clean:
	@echo "Cleaning $(SUBDIRS)..."
	$(MAKE) -C $(SUBDIRS) clean