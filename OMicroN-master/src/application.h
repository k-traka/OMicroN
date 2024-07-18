/* application.h
   OMicroN (optimising microstructures numerically) simulation program
   header file containing application class definitions and implementation
*/


#ifndef Application_H
#define Application_H

#include <string>
#include <memory>
#include "settings.h"
#include "microstructure.h"
#include <loguru.hpp>

/**
 * @class Application
 * @brief Manages the overall application flow.
 *
 * This class handles the initialization of settings, running the simulation based on these settings, 
 * initializing the log file, managing the microstructure class and performing all the simulation steps.
 */
class Application {
public:
    /**
     * @brief Constructs an Application object with a configuration file.
     * 
     * @param configFile Path to the configuration file.
     */
    Application(const std::string& configFile);

    /**
     * @brief Runs the application.
     */
    void run();

private:
    /**
     * @brief Initializes global variables and objects.
     */
    void initializeGlobal();

    UserSettings userSettings; ///< User settings for the application.
    std::unique_ptr<Microstructure> microstructure; ///< Pointer to the microstructure object.
};

#endif // application_H
