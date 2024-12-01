all: standards couplings state modeAmplitude tangent newton agp saveFile fputBreather

standards: src/standards.cpp include/standards.hpp include/eigenClasses.hpp
	g++ -c src/standards.cpp -o lib/standards.o -O3 -IF:\\Libraries\\json\\include -IF:\\Libraries\\eigen -I.\\include

couplings: src/couplings.cpp include/couplings.hpp include/constants.hpp include/eigenClasses.hpp include/standards.hpp
	g++ -c src/couplings.cpp -o lib/couplings.o -O3 -IF:\\Libraries\\json\\include -IF:\\Libraries\\eigen -I.\\include

state: src/state.cpp include/state.hpp include/eigenClasses.hpp include/couplings.hpp include/standards.hpp include/modeAmplitude.hpp
	g++ -c src/state.cpp -o lib/state.o -O3 -IF:\\Libraries\\json\\include -IF:\\Libraries\\eigen -I.\\include

modeAmplitude: src/modeAmplitude.cpp include/modeAmplitude.hpp include/eigenClasses.hpp include/standards.hpp include/couplings.hpp
	g++ -c src/modeAmplitude.cpp -o lib/modeAmplitude.o -O3 -IF:\\Libraries\\json\\include -IF:\\Libraries\\eigen -I.\\include

tangent: src/tangent.cpp include/tangent.hpp include/eigenClasses.hpp include/standards.hpp include/couplings.hpp include/state.hpp
	g++ -c src/tangent.cpp -o lib/tangent.o -O3 -IF:\\Libraries\\json\\include -IF:\\Libraries\\eigen -I.\\include

newton: src/newton.cpp include/newton.hpp include/eigenClasses.hpp include/standards.hpp include/couplings.hpp include/state.hpp include/tangent.hpp include/modeAmplitude.hpp
	g++ -c src/newton.cpp -o lib/newton.o -O2 -IF:\\Libraries\\json\\include -IF:\\Libraries\\eigen -I.\\include

agp: src/agp.cpp include/agp.hpp include/eigenClasses.hpp include/standards.hpp include/state.hpp
	g++ -c src/agp.cpp -o lib/agp.o -O3 -IF:\\Libraries\\json\\include -IF:\\Libraries\\eigen -I.\\include

saveFile: src/saveFile.cpp include/saveFile.hpp include/eigenClasses.hpp include/standards.hpp include/state.hpp include/agp.hpp
	g++ -c src/saveFile.cpp -o lib/saveFile.o -O2 -IF:\\Libraries\\json\\include -IF:\\Libraries\\eigen -I.\\include

fputBreather: src/fputBreather.cpp include/fputFixed.hpp
	g++ -c src/fputBreather.cpp -o lib/fputBreather.o -O3 -IF:\\Libraries\\json\\include -IF:\\Libraries\\eigen -I.\\include
	g++ -o _fputBreather.exe lib/fputBreather.o lib/saveFile.o lib/agp.o lib/newton.o lib/tangent.o lib/modeAmplitude.o lib/state.o lib/couplings.o lib/standards.o
