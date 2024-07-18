#ifndef UserSettings_H
#define UserSettings_H

#include <string>
#include <unordered_map>
#include <memory>

class UserSettings {
public:
    UserSettings(const std::string& inputFile);
    void writeToFile(const std::string& outputFile);
    void printParameters();

    // Member variables for the parameters
    std::unique_ptr<std::string> param1;
    std::unique_ptr<std::string> param2;
    std::unique_ptr<std::string> param3;

private:
    void readFromFile(const std::string& inputFile);
    void setDefaultValues();
};

#endif // UserSettings_H
