RACK_DIR ?= $(HOME)/Rack2-SDK/Rack-SDK

FLAGS +=
CFLAGS +=
CXXFLAGS +=
LDFLAGS +=

SOURCES += $(wildcard src/*.cpp)

DISTRIBUTABLES += res
DISTRIBUTABLES += $(wildcard LICENSE*)
DISTRIBUTABLES += $(wildcard README*)
DISTRIBUTABLES += $(wildcard patches)

include $(RACK_DIR)/plugin.mk

# MinGW binutils 2.30 ld segfaults on the DWARF debug info in these objects;
# disable debug info for Windows cross-builds (appended after plugin.mk so
# -g0 lands after compile.mk's -g and wins).
ifdef ARCH_WIN
  FLAGS += -g0
endif
