/* main.cpp
   OMicroN (optimising microstructures numerically) simulation program
   cpp-file containing main implementation
*/


/* main.cpp is used to creat the instance application which initiates a new simulation
*/

#include "application.h"
#include <iostream>

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <input_file>" << std::endl;
        return 1;
    }

    std::string inputFile = argv[1];
    Application app(inputFile);
    app.run();
    return 0;
}
