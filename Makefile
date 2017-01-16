program = black-derman-toy
%.o : %.cpp
	$(CXX) -g -c  $(CPPFLAGS) $(CXXFLAGS) $< -o $@

objects = $(program).o 

bdtmain : $(objects)
	$(CXX)  $(objects) -lnlopt -lm -o $(program)


