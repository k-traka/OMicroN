#include "hdf5_utils.h"
#include <iostream>

H5_IO::H5_IO(const std::string& fileName, const std::string& mode) {
    if (mode == "r") {
        file = new H5::H5File(fileName, H5F_ACC_RDONLY);
    } else {
        // Handle other modes if needed
    }
}

H5_IO::~H5_IO() {
    if (file) {
        file->close();
        delete file;
    }
}

template <typename T>
std::vector<T> H5_IO::ReadDataset(const std::string& datasetName) {
    std::vector<T> data;

    if (!file) return data;

    try {
        H5::DataSet dataset = file->openDataSet(datasetName);
        H5::DataSpace dataspace = dataset.getSpace();

        int ndims = dataspace.getSimpleExtentNdims();
        hsize_t dims_out[ndims];
        dataspace.getSimpleExtentDims(dims_out, NULL);

        data.resize(dims_out[0]);
        dataset.read(data.data(), H5::PredType::NATIVE_DOUBLE);
    } catch (H5::Exception& e) {
        std::cerr << "Error reading dataset: " << e.getCDetailMsg() << std::endl;
    }

    return data;
}

// Explicit instantiation of the template
template std::vector<double> H5_IO::ReadDataset<double>(const std::string& datasetName);
template std::vector<float> H5_IO::ReadDataset<float>(const std::string& datasetName);
template std::vector<int> H5_IO::ReadDataset<int>(const std::string& datasetName);
