all:
	#g++ HggUtils.cc -g -c -o HggUtils.o `root-config --cflags`
	#g++ main.cc -g -c -o main.o `root-config --cflags`
	#g++ *.o `root-config --cflags --glibs` NanoTools/NanoCORE/NANO_CORE.so -o looper.exe
	g++ main.cc -o main.exe `root-config --cflags --glibs` ../NanoCORE/NANO_CORE.so
clean:
	rm main.exe
	rm *.o

