# determine the hostname
HOST=$(shell hostname -s)

# make it generic for continous integration on github
ifeq ($(CIRCLECI),true)
    HOST=circle
endif

# directory of local definitions
LOCALDEFSDIR  = $(STATFEMROOT)/config

# filename of local definitions file
LOCALDEFSFILE = $(LOCALDEFSDIR)/$(addsuffix .$(HOST), statFEM_definitions)

# default file
LOCALDEFSDEFAULT = statFEM_definitions.default

# check existence of definitions file
# "test -f $(LOCALDEFSFILE)" has exit status 1 (failed) and or 0 (successful)
# thus "&& echo 'true'" is only executed, when file is found,
# otherwise call branches to "|| echo 'false'"
ifeq ($(shell test -f $(LOCALDEFSFILE) && echo 'true' || echo 'false'),false)
# put default filename since 
	LOCALDEFSFILE = $(LOCALDEFSDIR)/$(LOCALDEFSDEFAULT)
endif

# inform user
$(info Using definitions file $(LOCALDEFSFILE))

# include local definitions into this file
include $(LOCALDEFSFILE)

