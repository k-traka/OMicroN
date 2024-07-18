#include <algorithm>
#include <cmath>
#include <cstring>
#include <ctime>
#include <fstream>
#include <map>
#include <numeric>
#include <omp.h>
#include <sstream>
#include <utility>
#include <iostream>

#ifndef _GL_LINUX_
#else
#ifdef WIN32
#include <cassert>
#include <windows.h>
#endif
#endif
#ifndef _GL_LINUX_
#pragma hdrstop
#pragma package(smart_init)
#endif


#ifndef M_PI
#define M_PI 3.14159265359
#endif

#ifndef M_SQRT2
#define M_SQRT2 1.4142135623730951
#endif

#ifndef M_SQRT3
#define M_SQRT3 1.7320508075688772
#endif

#ifndef M_SQRT2_2
#define M_SQRT2_2 0.70710678118654755
#endif

#ifndef M_SQRT3_3
#define M_SQRT3_3 0.57735026918962573
#endif


#include "simulation.h"



int main(int argc, char* argv[]){
                UserSettings* user = new UserSettings();

                simulation::Run();

}
