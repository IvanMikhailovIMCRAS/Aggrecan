# имя программы:
TARGET=../Aggrecan

# исходники для сборки:
SOURCES = \
        CommonParameters.f90 \
        ErrorList.f90 \
	InputData.f90 OutputData.f90 \
        Walking.f90 \
        OnePoint.f90 \
	Aggrecan.f90 \

OBJECTS=$(SOURCES:%.f90=%.o)

TYPE=master

# про флаги:
# https://www.opennet.ru/docs/RUS/cpp/cpp-10.html

# простая сборка:
all: $(TARGET)

$(OBJECTS): $(SOURCES)

$(TARGET): $(OBJECTS)
	gfortran -O3 -o $(TARGET) -cpp $(SOURCES) 
clean:
	$(RM) $(TARGET) *.mod
	
.PHONY: all clean