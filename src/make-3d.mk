CC = gcc 
CFLAGS = -I. -O2
LINKER = gcc
LFLAGS = -Wall -I. -O2 -lm -lgsl -lgslcblas
TARGET = 3d-growth-onerun

BINDIR = ../bin
OBJDIR = ../obj

LOCAL_SRC = 3d-growth-onerun.c auxillary-funcs.c basis.c conservation-funcs.c \
            gibbs-thomson.c initial-conditions.c interaction-matrix.c \
            radius-chi-update.c

INCLUDES := headerfile.h


LOCAL_OBJS := $(LOCAL_SRC:%.c=$(OBJDIR)/%.o)
OBJECTS := $(LOCAL_OBJS)
rm = rm -f

$(BINDIR)/$(TARGET): $(OBJECTS)
	@$(LINKER) $(OBJECTS) $(LFLAGS) -o $@
	@echo "Linking complete!"

$(LOCAL_OBJS) : $(OBJDIR)/%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@
	@echo "Compiled "$<" successfully!"


PHONY: clean
clean:
	$(rm) $(OBJECTS)
	@echo "Cleanup complete!"

.PHONY: remove
remove: clean
	$(rm) $(BINDIR)/$(TARGET)
	@echo "Executable removed!"
