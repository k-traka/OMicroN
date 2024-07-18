#include "simulation.h"
#include <cstdlib>


#include <iostream>
#include <fstream>
 #include <cstring>
 
#define EQ_STR(s1, s2) (std::strcmp(s1, s2) == 0)

/** @brief Runs CASIPT simulation
 *
 * @param global Instance of global
 */
void simulation::Run(){
std::cout<<"hey"<<std::endl;
float a = 3.;
float mvTimeStep = a;
std::cout<<"time step is "<<mvTimeStep<<std::endl;
std::cout<<"wrote file "<<std::endl;
}



void simulation::Run(const char* settingsFilePath)
{

    UserSettings* user = new UserSettings(settingsFilePath);
    Simulation::Run(user);
    delete user;
}
