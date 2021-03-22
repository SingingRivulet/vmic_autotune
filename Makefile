vmic:main.cpp note.hpp freq.hpp
	g++ -std=c++11 main.cpp -o vmic  -Bstatic -lsox -Bstatic -lSoundTouch
clean:
	rm vmic
