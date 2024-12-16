//
// Created by Aamal Hussain on 16/12/2024.
//

#include "PXH5Dataset.h"
#include <highfive/H5Easy.hpp>


void PXH5Dataset::PrintGroupStructure(const Group &group, const std::string &prefix) {
    std::vector<std::string> keys = group.listObjectNames();

    for (const auto &key : keys) {
        std::cout << prefix << key;
        const ObjectType object_type = group.getObjectType(key);

        switch (object_type) {
            case ObjectType::Group:
                PrintGroupStructure(group.getGroup(key), prefix + " ");
            break;
            case ObjectType::Dataset:
                std::cout << " - Type: Dataset" << std::endl;
            break;
            default:
                std::cout << " - Type: Unknown" << std::endl;
            break;
        }
    }
}

void PXH5Dataset::ReadH5(const fs::path &filepath) {
    try {
        if (!fs::exists(filepath))
            throw std::runtime_error("File not found: " + static_cast<std::string>(filepath));

        const File file(filepath,  File::ReadOnly);
        const Group root_group = file.getGroup("/");


        verts_ = LoadMatrixDataset<MatrixD>(file, "/" + dataset_parameters_.verts);
        faces_ = LoadMatrixDataset<MatrixI>(file, "/" + dataset_parameters_.faces);


    } catch (const std::exception& err) {
        std::cerr << "Failed to load HDF5 file " << filepath << ": " << err.what() << std::endl;
    }
}

void PXH5Dataset::WriteH5(const fs::path &filepath, const MatrixD &out_verts, const MatrixI &out_faces) const {
    try {
        File file(filepath, File::OpenOrCreate);

        Group root_group = file.getGroup("/");
        H5Easy::dump(file, "/" + dataset_parameters_.verts, out_verts);
        H5Easy::dump(file, "/" + dataset_parameters_.faces, out_faces);

    } catch (std::exception &err) {
        std::cerr << "Failed to write to HDF5 file " << filepath << ": " << err.what() << std::endl;
    }
}

template <typename  T>
T PXH5Dataset::LoadMatrixDataset(const File& file, const std::string& path) {
    if (!file.exist(path))
        throw std::runtime_error("Dataset not found: " + path);
    return H5Easy::load<T>(file, path);
}
