CC=g++
OUT=main.exe
OBJS=main.cc

all: $(OUT)

$(OUT): $(OBJS) 
	$(CC) $(OBJS) -o $(OUT) `root-config --cflags --glibs` ../NanoCORE/NANO_CORE.so
	
clean:
	rm main.exe
	#rm *.o

