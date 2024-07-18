#ifndef HDF5_UTILS_H
#define HDF5_UTILS_H

#include <vector>
#include <string>
#include <H5Cpp.h>

class H5_IO {
public:
    H5_IO(const std::string& fileName, const std::string& mode);
    ~H5_IO();

    template <typename T>
    std::vector<T> ReadDataset(const std::string& datasetName);

private:
    H5::H5File* file;
};

#endif // HDF5_UTILS_H
