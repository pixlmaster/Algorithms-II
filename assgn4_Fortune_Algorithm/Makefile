#Name:Akash Tiwari
#Roll no:17CS10003

a.out: main.o Fortune.o
	g++ main.o Fortune.o
main.o: main.cpp Fortune.o
	g++ -c main.cpp
Fortune.o: Fortune.cpp Fortune_header.hpp
	g++ -c 	Fortune.cpp
clean:
	rm a.out main.o Fortune.o t4.svg
