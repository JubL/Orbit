CFLAGS = --std=c++14 -march=native -Wall -Wextra -O3
OUTPUTNAME = run_orbit


all: functions.o main.o models.o obj.o pbar.o
	g++ $(CFLAGS) functions.o main.o models.o obj.o pbar.o -o $(OUTPUTNAME) -lncurses -lpthread
	rm *.o

cross_section.o: cross_section.cpp
	g++ $(CFLAGS) -c cross_section.cpp

functions.o: functions.cpp
	g++ $(CFLAGS) -c functions.cpp

hello_world.o: hello_world.cpp
	g++ $(CFLAGS) -c hello_world.cpp

main.o: main.cpp
	g++ $(CFLAGS) -c main.cpp

models.o: models.cpp
	g++ $(CFLAGS) -c models.cpp

obj.o: obj.cpp
	g++ $(CFLAGS) -c obj.cpp

pbar.o: pbar.cpp
	g++ $(CFLAGS) -c pbar.cpp

hello: functions.o hello_world.o models.o obj.o
	g++ $(CFLAGS) functions.o hello_world.o models.o obj.o -o run_hello -lpthread -lncurses
	rm *.o

cross: cross_section.o functions.o obj.o
	g++ $(CFLAGS) cross_section.o functions.o obj.o -o run_cross
	rm *.o

data:
	g++ $(CFLAGS) getData.cpp -o run_data

clean:
	rm -f *.o $(OUTPUTNAME) run_cross run_hello run_data
