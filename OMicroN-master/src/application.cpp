/* application.cpp
   OMicroN (optimising microstructures numerically) simulation program
   cpp-file containing application class implementation
*/



#include "application.h"
#include <iostream>
#include <filesystem> // for std::filesystem::create_directories and std::filesystem::remove

/**
 * @brief Constructs an Application object with a user input (parameter) file.
 * 
 * @param inputFile Path to the input (parameter) file.
 */

/* Application sets the new simulation based on the input file, creates the instance userSettings, and initalized the log file to write the simulation evolution
 */

Application::Application(const std::string& inputFile) : userSettings(inputFile) {
    // Loguru initialization
    int argc = 1;
    char arg0[] = "app";
    char* argv[] = {arg0, nullptr};
    loguru::init(argc, argv);

    std::string outputFolderPath = userSettings.getOutputFolderPath();
    // Ensure the output directory exists
    std::filesystem::create_directories(outputFolderPath);
    std::string logFilePath = outputFolderPath + "/simulation.log";

    // Remove existing log file if it exists
    if (std::filesystem::exists(logFilePath)) {
        std::filesystem::remove(logFilePath);
    }

    // Set up new log file
    loguru::add_file(logFilePath.c_str(), loguru::Append, loguru::Verbosity_MAX);
}

/**
 * @brief Runs the application.
 */
void Application::run() {
    userSettings.printParameters();
    userSettings.writeToFile("ParametersForSimulation.txt");
    initializeGlobal();
    microstructure->printMicrostructureParameters();
    while (microstructure->IsStillSimulating())
        microstructure->SimulationStep();
}

/**
 * @brief Initializes global variables and objects.
 */
void Application::initializeGlobal() {
    microstructure = std::make_unique<Microstructure>(userSettings);
}
