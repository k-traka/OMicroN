/// @file src/simulation.h


#include "vendor/loguru.hpp"
#include "settings.h"


// #include "user.h"

namespace simulation {
void Run();
void Run(UserSettings* user);

// void Run(UserSettings* user, RandomNumberGenerator* rng);
// void Run(UserSettings* user);
// void Run(const char* settingsFilePath);

bool AddLogFile(const char* path, const char* mode = "w", int verbosity = 9);
void StopLogging();
};

