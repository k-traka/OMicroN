#include "UserSettings.h"
#include <iostream>
#include <fstream>
#include <sstream>

UserSettings::UserSettings(const std::string& inputFile) {
    setDefaultValues();
    readFromFile(inputFile);
}

void UserSettings::setDefaultValues() {
    param1 = std::make_unique<std::string>("default1");
    param2 = std::make_unique<std::string>("default2");
    param3 = std::make_unique<std::string>("default3");
}

void UserSettings::writeToFile(const std::string& outputFile) {
    std::ofstream out(outputFile);
    if (!out.is_open()) {
        std::cerr << "Unable to open file: " << outputFile << std::endl;
        return;
    }

    out << "param1=" << *param1 << std::endl;
    out << "param2=" << *param2 << std::endl;
    out << "param3=" << *param3 << std::endl;

    out.close();
}

void UserSettings::printParameters() {
    std::cout << "param1 = " << *param1 << std::endl;
    std::cout << "param2 = " << *param2 << std::endl;
    std::cout << "param3 = " << *param3 << std::endl;
}

void UserSettings::readFromFile(const std::string& inputFile) {
    std::ifstream in(inputFile);
    if (!in.is_open()) {
        std::cerr << "Unable to open file: " << inputFile << std::endl;
        return;
    }

    std::string line;
    while (std::getline(in, line)) {
        std::istringstream iss(line);
        std::string key, value;
        if (std::getline(iss, key, '=') && std::getline(iss, value)) {
            if (key == "param1") {
                param1 = std::make_unique<std::string>(value);
            } else if (key == "param2") {
                param2 = std::make_unique<std::string>(value);
            } else if (key == "param3") {
                param3 = std::make_unique<std::string>(value);
            }
        }
    }

    in.close();
}
