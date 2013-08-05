TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

LIBS += -lquadmath

SOURCES += main.cpp \
    functor1d.cpp \
    imc_dfm.cpp \
    diffusion1dp.cpp \
    boundary.cpp \
    tdma.cpp

HEADERS += \
    functor1d.h \
    imc_dfm.h \
    diffusion1dp.h \
    boundary.h \
    tdma.h

