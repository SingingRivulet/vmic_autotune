TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

LIBS += -Bstatic -lsox -lSoundTouch

SOURCES += main.cpp \
    kalman.cpp

HEADERS += \
    freq.hpp \
    note.hpp \
    kalman.h
