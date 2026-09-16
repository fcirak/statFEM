SUBDIRS := del2/apps/user/inference1D

.PHONY: all $(SUBDIRS)

all:
	@echo "Choose a target using VSCode Makefile Tools."

# Each subdir builds itself
$(SUBDIRS):
	$(MAKE) -C $@

clean:
	@echo "Cleaning $(SUBDIRS)..."
	$(MAKE) -C $(SUBDIRS) clean